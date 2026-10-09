function testLayout_Radar_scan
%%

l = qd_layout;
l.simpar.show_progress_bars = 0;
l.no_rx = 4;
l.randomize_rx_positions(100,0,0,0);
l.rx_position = [ 100,0,-100,0 ; 0,100,0,-100 ; 0,0,0,0  ];

t = qd_track('linear',0,0);
t.name = 'Tx';
t.no_snapshots = 4;
t.orientation = [ 0,0,0,0 ; 0,0,0,0 ; 0,2*pi/3,-2*pi/3,0 ];

assertEqual( get_length(t),0 );

x1 = t.interpolate('snapshot',1,[0,36;1,4]);
x2 = t.interpolate('snapshot',1,[0,36;1,4],[],true);
assertEqual( x1,x2 );

assertEqual( t.no_snapshots, 37 );
assertEqual( t.movement_profile, [0 36 ; 1 37] );
l.tx_track = t;

l.set_scenario('LOSonly');

l.tx_array = qd_arrayant('custom',60,60,0);
gain = l.tx_array.calc_gain;

c = l.get_channels;

P = [ c(1,1).coeff(:),c(2,1).coeff(:),c(3,1).coeff(:), c(4,1).coeff(:) ];
P = abs( P ).^2 ./ 10^(0.1*gain); 

assertTrue( all( abs( [ P(1,1) , P(37,1) , P(10,2) , P(19,3) , P(28,4) ] - 1  ) < 1e-7 ) );
assertTrue( all( abs( [ P(4,1) , P(34,1) , P(7,2) , P(13,2) , P(16,3) , P(22,3) , P(25,4), P(31,4) ] - 0.5  ) < 1e-7 ) );
assertTrue( all( [ P(19,1) , P(28,2) , P(1,3) , P(37,3) , P(10,4) ] < 1e-7 ) );


