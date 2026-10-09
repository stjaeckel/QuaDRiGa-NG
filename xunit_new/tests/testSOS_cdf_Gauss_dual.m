function testSOS_cdf_Gauss_dual
%%
% Make sure that the distribution of the gennerted random variables matches the requirement

set_rand_state( 1 );

% Perfect CDF (Gauss-Normal)
bins = -3:0.01:3;
val  = exp(-0.5*bins.^2);
cdf  = cumsum(val) / sum(val);

% Generate random variables 
posA = rand(3,10000) * 10000;         % Random 3D positions
posB = rand(3,10000) * 10000;         % Random 3D positions

v = qd_sos.randn( 10, posA, posB );        % Random spatially correlated variable
c = qf.acdf(v,bins);                % CDF
assertTrue(  all( abs(c-cdf') < 0.02 ) );