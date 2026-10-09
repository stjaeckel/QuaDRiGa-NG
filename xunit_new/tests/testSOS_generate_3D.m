function testSOS_generate_3D
%%
set_rand_state( 1 );

dcorr = 9.99;
max_range = 50;

D = (0:99) * max_range/100;
R = exp( -D/dcorr );

x = qd_sos.generate( R,D,50,3,1,10,1,0 );

assertTrue( x.dist_decorr - 10 < 1e-5 )
assertTrue( x.dimensions == 3 )

Ra = x.acf_approx;

assertTrue( all( size(Ra) - [397,397,3] < 1e-5 ))
assertTrue( all( abs( Ra(199,199,:) - 1 ) < 1e-5) )