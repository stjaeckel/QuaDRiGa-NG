function testSOS_remove_bias

set_rand_state( 1 );
ca = rand( 3,200 ) * 50;
cb = rand( 3,200 ) * 50;

% The correction is exact in double precision; the AVX2 kernel of quadriga_lib.sos_val evaluates
% sine and cosine in single precision, which leaves a residual mean of about 1e-8
tol = 1e-6;

% One device: the mean over the positions becomes zero, both phase columns stay identical
x = qd_sos( 'Comb256', 'Normal', 10 );
assertTrue( abs( mean( val( x, ca ) ) ) > 1e-3 );
remove_bias( x, ca );
assertTrue( abs( mean( val( x, ca ) ) ) < tol );
assertEqual( x.sos_phase(:,1), x.sos_phase(:,2) );

% Two devices on the same map
x = qd_sos( 'Comb256', 'Normal', 10 );
remove_bias( x, ca, cb );
assertTrue( abs( mean( val( x, ca, cb ) ) ) < tol );
assertEqual( x.sos_phase(:,1), x.sos_phase(:,2) );

% Two devices on different maps
x = qd_sos( 'Comb256', 'Normal', 10 );
init( x, 0 );
remove_bias( x, ca, cb );
assertTrue( abs( mean( val( x, ca, cb ) ) ) < tol );

% Different maps without second positions: the first device has a zero mean
x = qd_sos( 'Comb256', 'Normal', 10 );
init( x, 0 );
remove_bias( x, ca );
assertTrue( abs( mean( val( x, ca ) ) ) < tol );

% Object arrays are processed element by element
x = qd_sos( 'Comb256', 'Normal', 10 );
x(1,2) = qd_sos( 'Comb256', 'Normal', 20 );
remove_bias( x, ca );
assertTrue( all( abs( mean( val( x, ca ), 2 ) ) < tol ) );

end
