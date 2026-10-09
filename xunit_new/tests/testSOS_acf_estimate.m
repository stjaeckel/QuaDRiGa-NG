function testSOS_acf_estimate

%%
x = qd_sos;
[ Re, De, Re_dual ] = x.acf_estimate( 100, 10:10:70 );