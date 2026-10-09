function testSOS_generate_1D
%%
set_rand_state( 1 );

dcorr = 9.99;
max_range = 50;

D = (0:99) * max_range/100;
R = exp( -D/dcorr );

x = qd_sos.generate( R,D,50,1,1,1,1,0 );

assertTrue( x.dist_decorr - 10 < 1e-5 )
assertTrue( x.dimensions == 1 )
assertTrue( x.no_coefficients == 50 )

Ra = x.acf_approx;
assertTrue( numel(Ra) == 397 )
assertTrue( abs( Ra(199) - 1 ) < 1e-5 )
assertTrue( all( Ra(199:298) - R < 0.1 ) )
