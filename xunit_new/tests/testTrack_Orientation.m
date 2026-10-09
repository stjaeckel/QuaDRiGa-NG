function testTrack_Orientation
%%
t = qd_track('linear',1,0);
t.positions(3,2) = 1;
t.calc_orientation;

assertTrue( all( abs( t.orientation(3,:) ) < 0.001 ) )
assertTrue( all( abs( t.orientation(2,:) - atan(1) ) < 0.001 ) )

t.interpolate_positions( 10 )

assertTrue( all( abs( t.orientation(3,:) ) < 0.001 ) )
assertTrue( all( abs( t.orientation(2,:) - atan(1) ) < 0.001 ) )


t = qd_track('linear',0);
t.orientation = [0;-pi/2;0];
t.calc_orientation([],-pi/8);
assertTrue( abs( -pi/2 - pi/8 - t.orientation(2) ) < 1e-9 )
assertTrue( abs( 0 - abs(t.orientation(3)) ) < 1e-9 )

t = qd_track('linear',0);
t.orientation = [0;-pi/2;0];
t.calc_orientation([],[],pi/8);
assertTrue( abs( -pi/2 - t.orientation(2) ) < 1e-9 )
assertTrue( abs( pi/8 - t.orientation(3) ) < 1e-9 )


