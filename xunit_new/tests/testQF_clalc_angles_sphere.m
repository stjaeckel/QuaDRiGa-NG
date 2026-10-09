function testQF_clalc_angles_sphere

az = [ 180, 170, -170 ];
el = [ 0,-10,10 ]+5;
pow = [ 1,0.5,0.5 ];

[ as, es, orientation, phi, theta ] = qf.calc_angular_spreads_sphere( az*pi/180, el*pi/180, pow );

orientation = orientation*180/pi;
as = as*180/pi;
es = es*180/pi;
phi = phi*180/pi;
theta = theta*180/pi;

assertTrue( all( abs( phi - [ 0 -14.1 14.1] ) < 1 ) )
assertTrue( all( abs( theta ) < 1 ) )

assertTrue( abs( orientation(1) + 45 ) < 1 );
assertTrue( abs( orientation(2) -5 ) < 1 );
assertTrue( abs( orientation(3) -180 ) < 1 );

assertTrue( all( abs( as - 10 ) < 1 ) )
assertTrue( all( abs( es ) < 1 ) )
