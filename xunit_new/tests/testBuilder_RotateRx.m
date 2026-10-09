function testBuilder_RotateRx
%% Set up a LOS scenario
% Transmitter at [0;0;0], rx at [30;0;0];
% Rx has ULA4 with element on y-axis
% Rotate the Rx around the axis and check the angles

% This only works with spherical waves !!!

t = qd_track;
t.initial_position = [3000;0;0];
t.no_snapshots = 17;
t.positions = zeros( size(t.positions));
t.orientation = zeros(3,17);
t.orientation(3,:) = (0:22.5:360)*pi/180;
t.name = 'Rx1';

b = qd_builder('LOSonly');
b.tx_position = [0;0;0];
b.rx_array = qd_arrayant('ula4');
b.rx_array.element_position(2,:) = [-1,-0.5,0,1]*b.simpar.wavelength;
b.tx_array = qd_arrayant('omni');
b.name = 'Tx1';
b.simpar.show_progress_bars = false;
b.simpar.use_absolute_delays = true;
b.rx_positions = [3000;0;0];
b.rx_track = t;

gen_parameters(b);
c = get_channels(b);

ang = unwrap(angle( squeeze(c.coeff)' ));
ang = ang - mean(mean(ang));

tmp = ang(1,:);
assertTrue(  all( abs( tmp(:) ) < 0.01 ) );

tmp = ang(5,:) - [2*pi pi 0 -2*pi];
assertTrue(  all( abs( tmp(:) ) < 0.01 ) );

tmp = ang(9,:);
assertTrue(  all( abs( tmp(:) ) < 0.01 ) );

tmp = ang(13,:) - [-2*pi -pi 0 2*pi];
assertTrue(  all( abs( tmp(:) ) < 0.01 ) );

tmp = ang(17,:);
assertTrue(  all( abs( tmp(:) ) < 0.01 ) );

