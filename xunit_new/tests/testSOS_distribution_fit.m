function testSOS_distribution_fit
%%
% Make sure that the distribution of the gennerted random variables matches the requirement

npos = 10000;
pos = [ rand( 1,npos ); rand( 1,npos ) ; zeros( 1,npos) ]*10000;

bins = -3:0.01:3;
cc = qf.acdf( qd_sos.randn( 10, pos ) , bins);
cr = qf.acdf(        randn( 1, npos ) , bins);
assertTrue( all(abs(cc-cr) < 0.05 ) )


bins = 0:0.00166:1;
cc = qf.acdf( qd_sos.rand( 10, pos ) , bins);
cr = qf.acdf(        rand( 1, npos ) , bins);
assertTrue( all(abs(cc-cr) < 0.05 ) )

cr = qd_sos.randi( 10, pos, 10 );
num = zeros( 1,10 );
for n = 1:10
    num(n) = sum( cr == n );
end
assertTrue( sum(num) == npos )
assertTrue(  all( num > npos / 10 * 0.7 ) )