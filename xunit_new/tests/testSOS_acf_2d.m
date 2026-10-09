function testSOS_acf_2d
%%
% Test 2D ACF interpolation

x = qd_sos;
a = x.acf_2d;
nD = numel( x.dist );

assertTrue( a( 2*nD-1, 2*nD-1 ) - 1 <= 1e-5 )
assertTrue( all( abs( a( 2*nD-1 , 2*nD-1:3*nD-3 ) - x.acf(1:end-1) ) <= 1e-5 ) )
assertTrue( all( abs( a( 2*nD-1 , 2*nD-1:-1:nD+1 ) - x.acf(1:end-1) ) <= 1e-5 ) )
assertTrue( all( abs( a( 2*nD-1:3*nD-3,2*nD-1 ).' - x.acf(1:end-1) ) <= 1e-5 ) )
assertTrue( all( abs( a( 2*nD-1:-1:nD+1,2*nD-1 ).' - x.acf(1:end-1) ) <= 1e-5 ) )
