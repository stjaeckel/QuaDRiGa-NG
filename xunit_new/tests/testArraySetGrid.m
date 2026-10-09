function testArraySetGrid

a = qd_arrayant('Dipole');
a(2) =  qd_arrayant('Custom',90,90,1);

set_grid( a, -pi:pi:pi, -pi/2:pi/2:pi/2 );

assertTrue( all( abs( a(1).azimuth_grid  -  (-pi:pi:pi) ) < 1e-5 ) );
assertTrue( all( abs( a(2).azimuth_grid  -  (-pi:pi:pi) ) < 1e-5 ) );
assertTrue( all( abs( a(1).elevation_grid  -  (-pi/2:pi/2:pi/2) ) < 1e-5 ) );
assertEqual(  size(a(1).element_position)  ,  [3 1] );
assertEqual(  size(a(1).Fa)  ,  [3 3] );
assertEqual(  size(a(1).Fb)  ,  [3 3] );
assertEqual(  size(a(1).coupling)  ,  [1 1] );
assertEqual(  a(1).no_az  , 3   );
assertEqual(  a(1).no_el  , 3   );




