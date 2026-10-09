function testBuilder_GroundReflection_Numeric_Stability
% The calculations depent on the phase difference between GR anad LOS
% With single precisoin, the resutls become incorrect ven at close distances
%%

t = qd_track([]);
d = 10.^( 3:7 );
d = reshape( [ d;2*d;5*d ], 1,[] );
t.positions = [d ; zeros(2,numel(d))];

b = qd_builder('TwoRayGR');
b.simpar.center_frequency = 2e9;
b.simpar.show_progress_bars = 0;
b.scenpar.GR_epsilon = 10;
b.rx_track = t;
b.rx_positions = [10,0,2]';
b.tx_position = [0;0;4];
gen_parameters(b);

assertTrue( all( abs( b.EoD- [-atan(2/10) -atan(6/10)]) < 1e-10 ) );
assertTrue( all( abs( b.EoA- [atan(2/10) -atan(6/10)]) < 1e-10 ) );

c = b.get_channels;

H = c.fr(100e6,64);
D  = c.rx_position(1,:);
PG = 10*log10(reshape( mean(abs(H).^2,3) , 1,[] ));
PGR = -40*log10(D) + 20*log10( b.rx_positions(3) * b.tx_position(3) );

assertTrue( all( abs( PG - PGR ) < 0.3 ) )