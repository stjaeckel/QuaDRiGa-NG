function testTrack_Generate
%%
t = qd_track.generate('linear',0,pi/2);
assertEqual( t.positions , [0,0,0]' );
assertEqual( t.orientation , [0;0;pi/2] );

t = qd_track.generate('linear',1,pi/2);
assertEqual( t.positions , [0,0;0,1;0,0] );

t.calc_orientation;
assertEqual( t.orientation , [0 0;0 0;pi/2 pi/2] );

t = qd_track.generate('circular',pi,pi/2);
tmp = t.positions(:,[1,33,65,97,129]) - [0,-0.5,0,0.5,0;0,-0.5,-1,-0.5,0;0,0,0,0,0];
assertTrue( all( abs( tmp(:) ) <0.001 ) );
assertEqual( t.closed , true );

t.calc_orientation;
assertTrue( all( abs( t.orientation(3,[1,33,65,97,129]) - [-pi,-pi/2,0,pi/2,-pi] ) < 0.05 ) );
assertEqual(  t.orientation(1,:) , zeros(1,129) );
assertEqual(  t.orientation(2,:) , zeros(1,129) );

assertEqual( t.initial_position , [0;0;0] );
