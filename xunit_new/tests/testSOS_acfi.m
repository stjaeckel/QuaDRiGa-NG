function testSOS_acfi
%%
x = qd_sos;
val = x.acfi( x.dist );
assertTrue( all(abs(x.acf-val)<1e-6) );
