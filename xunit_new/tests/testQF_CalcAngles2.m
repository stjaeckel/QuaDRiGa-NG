function testQF_CalcAngles2

a = qd_arrayant('testarray',2);
a.center_frequency = 2e9;
single(a);

% a = a.sub_array(1:28);  % Only V
% 
% a.Fb = a.Fa; a.Fa(:) = 0;   % Only H

% Small noise on Fb
%a.Fb = (randn(a.no_el,a.no_az,a.no_elements) + 1j*randn(a.no_el,a.no_az,a.no_elements))*1e-20;

ar = qd_arrayant('omni');
ar.center_frequency = a.center_frequency;
ar.copy_element(1,2:4);
ar.rotate_pattern(45,'x',2);
ar.rotate_pattern(180,'x',3);
ar.rotate_pattern(240,'x',4);

b = qd_builder('Freespace');
b.simpar.center_frequency = a.center_frequency;
b.simpar.show_progress_bars = 0;

ang = [0,0.1,1,2,5];
N = numel(ang);
d3d = b.simpar.wavelength*500;

b.tx_position = [0;0;1];
b.rx_track = qd_track('linear',1,0);
b.rx_track.positions = [ cosd(ang)*d3d, cosd(ang(2:end))*d3d; sind(ang)*d3d, zeros(1,N-1); zeros(1,N), sind(ang(2:end))*d3d ] - repmat([d3d;0;0],1,2*N-1);
b.rx_track.initial_position = [d3d;0;1];
b.tx_array = a;
b.rx_array = ar;

b.gen_parameters;

add_sdc( b, [-d3d;0;1], 0, 'absolute', 'freespace', 100, 0, 1 );
b.pin(:) = 0;

c = b.get_channels;
c.swap_tx_rx;

[ A, E, J, P, R, RX_beam, RX_coeff ] = qf.calc_angles( c.coeff, a, 2, [], [], N+1, 0  );

% Output Variable sizes
assertEqual( size(A), [ar.no_elements,2,2*N-1,2] );
assertEqual( size(E), [ar.no_elements,2,2*N-1,2] );
assertEqual( size(J), [2,ar.no_elements,2,2*N-1,2] );
assertEqual( size(P), [ar.no_elements,2,2*N-1,2] );
assertEqual( size(R), [ar.no_elements,2,2*N-1] );
assertEqual( size(RX_beam),  [a.no_el,a.no_az,2*N-1] );
assertEqual( size(RX_coeff), [a.no_elements,ar.no_elements,2,2*N-1,2] );

assertFalse( any(isnan(RX_beam(:))));
assertFalse( any(isnan(RX_coeff(:))));

% Second path removed
assertTrue( all(abs( reshape(P(:,:,:,2),[],1) ) < 1e-8))

% Check angles
oo = ones(1,ar.no_elements);
x = squeeze(A(:,1,1:N,1))*180/pi-ang(oo,:); assertTrue( all(abs(x(:)) < 0.05) );
x = squeeze(A(:,1,N+1:end,1))*180/pi; assertTrue( all(abs(x(:)) < 0.05) );
x = squeeze(E(:,1,1:N,1))*180/pi; assertTrue( all(abs(x(:)) < 0.05) );
x = squeeze(E(:,1,N+1:end,1))*180/pi-ang(oo,2:end); assertTrue( all(abs(x(:)) < 0.05) );

% Angles of reflected path
x = squeeze(A(:,2,1:N,1))+pi; x = angle(exp(1j*x))*180/pi; assertTrue( all(abs(x(:)) < 0.05) );

% Check path power
x = -10*log10(squeeze(P(:,1,:,1)))+10*log10(b.gain(1)); assertTrue( all(abs(x(:)) < 0.05) );
x = -10*log10(squeeze(P(:,2,:,1)))+10*log10(b.gain(2)); assertTrue( all(abs(x(:)) < 0.05) );

% Check polarizatoion
x = imag( squeeze(J(1,:,1,:,1)) ); assertTrue( all(abs(x(:)) < 1e-2) ); % No phase

Jx = exp(1j*pi*[0,45,180,240]/180).';
Jv = squeeze(real(J(1,:,1,:,1)));
Jh = squeeze(real(J(2,:,1,:,1)));

on = ones(1,2*N-1);
x = Jv-real(Jx(:,on)); assertTrue( all(abs(x(:)) < 2e-2) );
x = Jh-imag(Jx(:,on)); assertTrue( all(abs(x(:)) < 2e-2) );

% Estimating multiple spatial paths
[ A, E, J, P, R, RX_beam, RX_coeff ] = qf.calc_angles( sum(c.coeff,3), a, 3, [], [], 2*N, 0  );

assertFalse( any(isnan(RX_beam(:))));
assertFalse( any(isnan(RX_coeff(:))));

% Check amunt of resolved energy
assertTrue( all(R(:) > 0.95) );

% 3rd and 4th sub-path should have no power
x = P(:,1,:,3); assertTrue( all( x(:) < 1e-40));

% V-Polarizion should detect path powers perfectly
x = -10*log10(squeeze(P(1,1,:,1)))+10*log10(b.gain(1)); assertTrue( all(abs(x(:)) < 0.15) ); % V up
x = -10*log10(squeeze(P(1,1,:,2)))+10*log10(b.gain(2)); assertTrue( all(abs(x(:)) < 0.15) );
x = -10*log10(squeeze(P(3,1,:,1)))+10*log10(b.gain(1)); assertTrue( all(abs(x(:)) < 0.15) ); % V down
x = -10*log10(squeeze(P(3,1,:,2)))+10*log10(b.gain(2)); assertTrue( all(abs(x(:)) < 0.15) );

% Power from H polarization gets transfered to V due to non-exisiting angle resolution on H
x = -10*log10(squeeze(P(2,1,:,1)))+10*log10(b.gain(1)); assertTrue( all(x(:) < -1) );
x = -10*log10(squeeze(P(4,1,:,1)))+10*log10(b.gain(1)); assertTrue( all(x(:) < -1) );

% No second spatial path on H-polarized Tx elements due to 13 dB cutoff limit (0.95)
x = P([2,4],1,:,2); assertTrue( all( x(:) < 1e-40));

% Check angles
x = squeeze(A([1,3],1,1:N,1))*180/pi-ang([1,1],:); assertTrue( all(abs(x(:)) < 0.2) );
x = squeeze(A(:,1,N+1:end,1))*180/pi; assertTrue( all(abs(x(:)) < 0.05) );
x = squeeze(E(:,1,1:N,1))*180/pi; assertTrue( all(abs(x(:)) < 0.05) );
x = squeeze(E([1,3],1,N+1:end,1))*180/pi-ang([1,1],2:end); assertTrue( all(abs(x(:)) < 0.3) );

% Check polarizatoion
x = imag( squeeze(J(1,:,1,:,1)) ); assertTrue( all(abs(x(:)) < 0.1) ); % No phase

Jv = squeeze(real(J(1,:,1,:,1)));
Jh = squeeze(real(J(2,:,1,:,1)));

assertTrue( all(abs(Jv(1,:) - 1) < 1e-4));
assertTrue( all(abs(Jv(3,:) + 1) < 1e-4));
assertTrue( all(all(abs(Jh([1,3],:)) < 1e-4)));


%%
% Test ill-conditioned antenna - V
a = a.sub_array(1:28);  % Only V
b.tx_array = a;
b.rx_track.no_snapshots = 1;

c = b.get_channels;
c.swap_tx_rx;

[ A, E, J, P, R, RX_beam, RX_coeff ] = qf.calc_angles( c.coeff, a, 2, [], [], N+1, 0  );

assertFalse( any(isnan(RX_beam(:))));
assertFalse( any(isnan(RX_coeff(:))));
assertTrue( all(R(:) > 0.95) );

% V-Polarizion should detect path powers perfectly
x = -10*log10(squeeze(P([1,3],1,:,1)))+10*log10(b.gain(1)); assertTrue( all(abs(x(:)) < 0.05) );    
x = -10*log10(squeeze(P([1,3],2,:,1)))+10*log10(b.gain(2)); assertTrue( all(abs(x(:)) < 0.05) );

% Second tx antenne should have 50% of the power du to loss of second polarization
x = -10*log10(squeeze(P(2,1,:,1)))+10*log10(0.5*b.gain(1)); assertTrue( all(abs(x(:)) < 0.05) );  
x = -10*log10(squeeze(P(2,2,:,1)))+10*log10(0.5*b.gain(2)); assertTrue( all(abs(x(:)) < 0.05) );  

% Check angles
assertTrue( all(abs(A(:,1,1,1)) < 1e-2))
assertTrue( all(abs(exp(1j*A(:,2,1,1)) + 1)  < 1e-4))
assertTrue( all(all(abs(E(:,:,1,1)) < 1e-2)))

%%
% Test ill-conditioned antenna - H
a.Fb = a.Fa; a.Fa(:) = 0;   % Only H
b.tx_array = a;
b.rx_track.no_snapshots = 1;

c = b.get_channels;
c.swap_tx_rx;

[ A, E, J, P, R, RX_beam, RX_coeff ] = qf.calc_angles( c.coeff, a, 2, [], [], N+1, 0  );

assertFalse( any(isnan(RX_beam(:))));
assertFalse( any(isnan(RX_coeff(:))));

% Second tx antenne should have 50% of the power du to loss of second polarization
x = -10*log10(squeeze(P(2,1,:,1)))+10*log10(0.5*b.gain(1)); assertTrue( all(abs(x(:)) < 0.05) );  
x = -10*log10(squeeze(P(2,2,:,1)))+10*log10(0.5*b.gain(2)); assertTrue( all(abs(x(:)) < 0.05) );  

% Check angles
assertTrue( all(abs(A([2,4],1,1,1)) < 1e-2))
assertTrue( all(abs(exp(1j*A([2,4],2,1,1)) + 1)  < 1e-4))
assertTrue( all(all(abs(E([2,4],:,1,1)) < 1e-2)))


