function testSOS_default
%%
% Test if the default CDF is initializesd correctly

x = qd_sos;
assertTrue( strcmp(x.name,'Comb256') )
assertTrue( strcmp(x.distribution,'Normal') )
assertTrue( x.dist_decorr == 10 )
assertTrue( x.dimensions == 3 )
assertTrue( x.no_coefficients == 256 )
assertTrue( x.acf( find( x.dist >= 10,1 ) ) - exp(-1) < 0.01 )
