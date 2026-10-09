function testBuilder_Pol_Rotation
%%
% One feature of the model is that it allows to freely orient the antennas
% at the transmitter and receiver. Here, this feature is tested. Two
% cross-polarized patch antennas were aligned on the optical axis facing
% each other. The surface normal vectors of the transmit and the receive
% patch are aligned with the LOS. The transmitter is rotated from -90° to
% 90° around the optical axis. The real and imaginary parts of the channel
% coefficients are then simulated for each angle. Each real and imaginary
% part is normalized by its maximum.
%
a = qd_arrayant('lhcp-rhcp-dipole');
a.append_array( qd_arrayant('custom',90,90,0) );
a.set_grid( (-180:10:180)*pi/180 , (-90:10:90)*pi/180 );
a.rotate_pattern(180,'z',3);
a.copy_element(3,4);
a.rotate_pattern(90,'x',4);

b = qd_builder('LOSonly');
b.rx_array = a;
b.tx_array = a.copy;
b.name = 'Tx1';
b.simpar.show_progress_bars = false;
b.tx_position = [0;0;0];
b.rx_positions = [11;0;0];
b.rx_track = qd_track.generate('linear',1,0);
b.rx_track.initial_position = b.rx_positions(:,1);
b.rx_track.name = 'Rx1';

gen_parameters(b);

rot = -135:45:135;
h = zeros(4,4,numel(rot),'single');
for n=1 : numel(rot)
    cc = a.copy;
    cc.rotate_pattern( rot(n) , 'x');
    b.tx_array = cc;
    c = b.get_channels;
    h(:,:,n) = c(1,1).coeff(:,:,1,1);
end

hh = h([3,4],[3,4],:);
g = hh;

% Main diagonal phases
tmp = g(1,1,[1,3,4,5,7])  -  g(2,2,[1,3,4,5,7]);
assertTrue(  all( abs( tmp(:) ) < 1e-8 ) );

% Side diagonal phases
tmp =  g(1,2,[1,2,3,5,6,7])  + g(2,1,[1,2,3,5,6,7]);
assertTrue(  all( abs( tmp(:) ) < 1e-8 ) );

g = real(hh)./max(real(reshape(hh,[],1))) + 1j*imag(hh)./max(imag(reshape(hh,[],1)));

% At angle 0, the array is aligned, at +/- 90°, they are inverted
tmp = abs(g(:,:,4)) - [sqrt(2),0;0,sqrt(2)];
assertTrue(  all( abs( tmp(:) ) < 1e-6 ) );

% Inverted powers
tmp = abs(g(:,:,2))  - [0,sqrt(2);sqrt(2),0];
assertTrue(  all( abs( tmp(:) ) < 1e-6 ) );

tmp = abs(g(:,:,6)) - [0,sqrt(2);sqrt(2),0];
assertTrue(  all( abs( tmp(:) ) < 1e-6 ) );

% Crossed powers
pow = [1,1;1,1]*abs(g(1,1,1));

tmp = abs(g(:,:,1))  - pow;
assertTrue(  all( abs( tmp(:) ) < 1e-6 ) );

tmp = abs(g(:,:,3))  - pow;
assertTrue(  all( abs( tmp(:) ) < 1e-6 ) );

tmp = abs(g(:,:,5))  - pow;
assertTrue(  all( abs( tmp(:) ) < 1e-6 ) );

tmp = abs(g(:,:,7))  - pow;
assertTrue(  all( abs( tmp(:) ) < 1e-6 ) );

% Circular Polarization results
hh = h([1,2],[1,2],:);
g = hh;

% Off-Diagonal elements must be 0
assertTrue( all(squeeze(abs(g(1,1,:))) < 1e-7)   );
assertTrue( all(squeeze(abs(g(2,2,:))) < 1e-7)   );

% Off-Diagonal-Power must be identical
tmp = squeeze(abs(g(1,2,:))) -  ones(7,1)*abs(g(1,2,1)) ;
assertTrue(  all( abs( tmp(:) ) < 1e-6 ) );

tmp = squeeze(abs(g(2,1,:))) -  ones(7,1)*abs(g(2,1,1));
assertTrue(  all( abs( tmp(:) ) < 1e-6 ) );

% Angle offsets between two adjacent points must be multiples of pi/2
aa = ( squeeze(angle(g(1,2,:))) - squeeze(angle(g(1,2,:))) )/(pi/2);
assertTrue( all( aa-round(aa) < 1e-8 ) );

