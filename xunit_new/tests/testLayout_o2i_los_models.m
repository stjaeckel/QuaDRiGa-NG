function testLayout_o2i_los_models
% O2I penetration loss of TR 38.901 clause 7.4.3 and LOS probability of InF-DH

set_rand_state( 4 );
f = [ 2, 28 ];                                      % GHz
tw = @( p1, L1, p2, L2 ) 5 - 10*log10( p1*10.^(-L1/10) + p2*10.^(-L2/10) );
L_glass = 2 + 0.2*f;
L_IRR = 25.4 + 0.11*f;
L_concrete = 5 + 4*f;

l = qd_layout;
l.simpar.show_progress_bars = 0;
l.simpar.center_frequency = f * 1e9;
l.no_rx = 500;
l.randomize_rx_positions( 500, 1.5, 1.5, 0 );

%% Low-loss model, uncorrelated, maximum indoor distance 10 m
o2i = l.gen_o2i_loss( '3GPP_38.901', 1, 0, 10 );
[ loss, d2din ] = collect( l, o2i );
assertTrue( all( d2din >= 0 & d2din <= 10 ) );
assertElementsAlmostEqual( loss(:,2) - loss(:,1), ( tw(0.3,L_glass(2),0.7,L_concrete(2)) - ...
    tw(0.3,L_glass(1),0.7,L_concrete(1)) ) * ones(500,1), 'absolute', 1e-10 );
dev = loss(:,1) - tw(0.3,L_glass(1),0.7,L_concrete(1)) - 0.5*d2din;
assertTrue( abs( mean(dev) ) < 0.6 && abs( std(dev) - 4.4 ) < 0.4 );     % Normal, sigma_P = 4.4 dB
assertTrue( abs( mean(d2din) - 10/3 ) < 0.4 );                           % Minimum of two uniforms

%% High-loss model (Release 19 IRR glass)
o2i = l.gen_o2i_loss( '3GPP_38.901', 0, 0, 25 );
loss = collect( l, o2i );
assertElementsAlmostEqual( loss(:,2) - loss(:,1), ( tw(0.7,L_IRR(2),0.3,L_concrete(2)) - ...
    tw(0.7,L_IRR(1),0.3,L_concrete(1)) ) * ones(500,1), 'absolute', 1e-10 );

%% Models of quadriga_lib.o2i_loss by name
o2i = l.gen_o2i_loss( 'car', [], 0 );
[ loss, d2din ] = collect( l, o2i );
assertTrue( all( d2din == 0 ) );
assertTrue( abs( mean( loss(:,1) ) - 9 ) < 0.7 && abs( std( loss(:,1) ) - 5 ) < 0.5 );

%% Spatial consistency with the decorrelation distances of TR 38.901
l.no_rx = 2;
l.rx_position = [ 100, 100.05 ; 50, 50 ; 1.5, 1.5 ];
o2i = l.gen_o2i_loss( '3GPP_38.901', 0.5 );
assertTrue( max( abs( o2i{1}(:) - o2i{2}(:) ) ) < 0.5 );

%% LOS probability of InF-DH (clutter density 0.6)
l = qd_layout;
l.simpar.show_progress_bars = 0;
l.tx_position = [ 0; 0; 8 ];
n = 1000;
l.no_rx = n;
phi = (0:n-1) / n * 2*pi;
l.rx_position = [ 5*cos(phi); 5*sin(phi); 1.5*ones(1,n) ];
l.set_scenario( '3GPP_38.901_InF_DH' );
is_los = false( 1,n );
for r = 1 : n
    is_los(r) = strcmp( l.rx_track(1,r).scenario{1}, '3GPP_38.901_InF_LOS' );
end
p_los = exp( -5 * (-log(1-0.6)/2) * (6-1.5)/(8-1.5) );                 % TR 38.901 Table 7.4.2-1
assertTrue( abs( mean( is_los ) - p_los ) < 0.03 );

end

function [ loss, d2din ] = collect( l, o2i )
% O2I loss [n_rx, n_freq] and indoor distance [n_rx, 1] of the first BS
loss = zeros( l.no_rx, size( o2i{1}, 3 ) );
d2din = zeros( l.no_rx, 1 );
for r = 1 : l.no_rx
    loss(r,:) = reshape( o2i{r}(1,1,:), 1, [] );
    d2din(r) = l.rx_track(1,r).par.o2i_d2din(1,1);
end
end
