function testBuilder_GroundReflection
%%
b = qd_builder('TwoRayGR');
b.scenpar.GR_epsilon = 10;        % Perfect dielectricum

b.simpar.show_progress_bars = false;
b.name = 'Tx1';

d = b.simpar.wavelength * 200;     % LOS path length = 2D distance
x = b.simpar.wavelength * 200.5;   % GR path length (0.5 lambda longer)
h = 0.5*sqrt(x^2-d^2);      % Height

a = qd_arrayant('ula2');    % V-Polarization
a.element_position(:) = 0;
b.rx_array = a;
b.tx_array = a;

b.tx_position = [0;0;h];
b.rx_positions = [d,0,h ; 1e4,0,h]';

% EoD for the GR
theta_r = atan(2*h/d);

% The reflection coefficient
epsilon_r = 10;
Z         = sqrt( epsilon_r - (cos(theta_r)).^2 );
R_par     = (epsilon_r .* sin(theta_r) - Z) ./ (epsilon_r .* sin(theta_r) + Z);
R_per     = ( sin(theta_r) - Z) ./ ( sin(theta_r) + Z);
R         = ( 0.5*(abs(R_par).^2 + abs(R_per).^2) ).';

gen_parameters(b);
c = get_channels(b);

% Test if the same epsilon as in the plpar was set
assertTrue( all( abs( b.gr_epsilon_r - b.scenpar.GR_epsilon ) < 1e-6 ) );

% Are there two paths?
assertTrue( all( cat(1,c.no_path) == 2 ) );

c(1,1).individual_delays = 0;

% Is there 0.5 lambda delay difference ?
assertTrue( c(1,1).delay(2,1) * b.simpar.speed_of_light / b.simpar.wavelength - 0.5 < 1e-6 );

% Both coefficients must be equally strong
A1 = abs( c(1,1).coeff(1,1,1) );
A2 = abs( c(1,1).coeff(1,1,2) );
assertTrue( abs(abs( A1 * R_par ) - A2) < 1e-4 );

% Coefficients must have same phase (180 deg shift due to reflection, 180 deg shift due to length difference)
P1 = angle(c(1,1).coeff(1,1,1,1))*180/pi;
P2 = angle(c(1,1).coeff(1,1,2,1))*180/pi;
assertTrue( abs(P1-P2) < 0.01  );

% Asymptotic path-gain must match the two-ray model
d = c(1,2).rx_position(1);
PGR = 40*log10(d) - 20*log10( h*h );
H = c(1,2).fr(100e6,64);
P = -10*log10( squeeze(mean(abs(H(1,1,:,:)).^2,3)) );
assertTrue(  abs(PGR - P) < 0.5 );

