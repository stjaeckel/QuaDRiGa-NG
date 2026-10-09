function testSOS_map
%%

set_rand_state( 1 );

% Implicitely tests "val"
x = qd_sos;
s = x.map(1:10,1:10,5);
si = quadriga_lib.interp( 1:10, 1:10, s, 1:0.1:10, 1:0.1:10 );

s = x.map(1:100:20000,1:100:20000,5);

% Perfect CDF (Gauss-Normal)
bins = -3:0.01:3;
val  = exp(-0.5*bins.^2);
cdf  = cumsum(val) / sum(val);

% Generate random variables 

c = qf.acdf(s(:),bins);                % CDF
assertTrue(  all( abs(c-cdf') < 0.02 ) );
