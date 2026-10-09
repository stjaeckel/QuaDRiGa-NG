function testSOS_dual_normalization
%%
% Perfect CDF (Gauss-Normal)
bins = -3:0.01:3;
val  = exp(-0.5*bins.^2);
cdf  = cumsum(val) / sum(val);

% Generate random variables 
pos = rand(3,10000) * 10000;         % Random 3D positions
v = qd_sos.randn( 10, pos, pos );        % Random spatially correlated variable
c = qf.acdf(v,bins);                % CDF
assertTrue(  all( abs(c-cdf') < 0.02 ) );