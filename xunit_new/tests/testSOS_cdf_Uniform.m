function testSOS_cdf_Uniform
%%
% Make sure that the distribution of the gennerted random variables matches the requirement

% Perfect CDF (Gauss-Normal)
bins = 0:0.001:1;
cdf  = bins;

% Generate random variables 
pos = randn(3,10000) * 10000;         % Random 3D positions
v = qd_sos.rand( 10, pos );        % Random spatially correlated variable
c = qf.acdf(v,bins);                % CDF

assertTrue(  all( abs(c-cdf') < 0.03 ) );