function testArayConstructor

a = qd_arrayant('omni');
assertEqual(  a.name  ,  'omni' );
assertTrue(  a.no_elements == 1  );
assertTrue( all( abs( a.elevation_grid  - (-90:90)*pi/180 ) < 1e-5 ) );
assertTrue( all( abs(  a.azimuth_grid  -  (-180:180)*pi/180  ) < 1e-5 ) );
assertTrue( all( a.element_position == [0;0;0] ) );
assertTrue( all( abs(  a.Fa(:)  - ones( 181*361,1) ) < 1e-5 ) );
assertTrue( all( abs(  a.Fb(:)  ) < 1e-5 ) );
assertTrue(  a.coupling  ==  1 );
assertTrue(  a.no_az == 361  );
assertTrue(  a.no_el == 181  );
