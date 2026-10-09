function testSOS_generate_2D
%%
set_rand_state( 1 );

dcorr = 9.99;
max_range = 50;

D = (0:99) * max_range/100;
R = exp( -D/dcorr );

x = qd_sos.generate( R,D,100,2,1,13,1,0 );

assertTrue( x.dist_decorr - 10 < 1e-5 )
assertTrue( x.dimensions == 2 )

Ra = x.acf_approx;

assertTrue( all( size(Ra) - [397,397] < 1e-5 ))
assertTrue( abs( Ra(199,199) - 1 ) < 1e-5 )
assertTrue( all( abs(Ra(199,199:298) - R) < 0.15 ) )
assertTrue( all( abs(Ra(199:298,199)'- R) < 0.15 ) )
