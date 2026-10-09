function testArrayInterpolation

a = qd_arrayant('dipole');
[Vi,H] = a.interpolate( 0.3*pi/180 , 0.3*pi/180 );
assertTrue( abs( a.Fa(91,181) - Vi ) > 1e-7 );