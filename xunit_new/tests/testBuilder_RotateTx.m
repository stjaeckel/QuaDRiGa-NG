function testBuilder_RotateTx
%% Set up a LOS scenario
% Transmitter at [0;0;0], rx at [30;0;0];
% Tx has ULA4 with element on y-axis
% Rotate the Tx around the z-axis and check the angles

t = qd_track;
t.initial_position = [3000;0;0];
t.no_snapshots = 1;
t.positions = zeros( size(t.positions));
t.name = 'Rx1';
t.orientation = [0;0;0];

b = qd_builder('LOSonly');
b.tx_position = [0;0;0];
b.rx_array  =qd_arrayant('omni');
b.tx_array =  qd_arrayant('ula4');
b.tx_array.element_position(2,:) = [-1,-0.5,0,1]*b.simpar.wavelength;
b.name = 'Tx1';
b.simpar.show_progress_bars = false;
b.simpar.use_absolute_delays = true;
b.rx_positions = [3000;0;0];
b.rx_track = t;
gen_parameters(b);

a = b.tx_array.copy;
a.set_grid( (-180:10:180)*pi/180 , (-90:10:90)*pi/180 );

rot = 0:22.5:360;
h = zeros(4,numel(rot),'single');
for n=1:numel(rot)
    cc = a.copy;
    cc.rotate_pattern( rot(n) , 'z');
    b.tx_array = cc;
    
    c = b.get_channels;
    h(:,n) = c.coeff(1,:,1);
end

ang = unwrap(angle( squeeze(h)' ));
ang = ang - mean(mean(ang));

tmp = ang(1,:);
assertTrue(  all( abs( tmp(:) ) < 0.01 ) );

tmp = ang(5,:) - [-2*pi -pi 0 2*pi];
assertTrue(  all( abs( tmp(:) ) < 0.01 ) );

tmp = ang(9,:);
assertTrue(  all( abs( tmp(:) ) < 0.01 ) );

tmp = ang(13,:) - [2*pi pi 0 -2*pi];
assertTrue(  all( abs( tmp(:) ) < 0.01 ) );

tmp = ang(17,:);
assertTrue(  all( abs( tmp(:) ) < 0.01 ) );