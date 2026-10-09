function testSOS_randn
%%
coord = [ 0,0,0 ; 0,0,0 ; 0,1,0; 0,0,1; 1,0,0 ]';

r = qd_sos.randn( 20, coord );
assertTrue( abs(r(1)-r(2)) < 1e-5 )
assertTrue( abs(r(1)-r(3)) > 1e-5 )
assertTrue( abs(r(1)-r(4)) > 1e-5 )
assertTrue( abs(r(1)-r(5)) > 1e-5 )
assertTrue( abs(r(3)-r(4)) > 1e-5 )
assertTrue( abs(r(3)-r(5)) > 1e-5 )
assertTrue( abs(r(4)-r(5)) > 1e-5 )
