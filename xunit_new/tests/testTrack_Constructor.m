function testTrack_Constructor
% Construct track and check defaults
%%
t = qd_track([]);
assertTrue( t.no_snapshots == 1);
assertEqual(  t.initial_position  ,  [0;0;0] );

t = qd_track('linear');
assertEqual(  t.name  ,  'track' );
assertEqual(  t.initial_position  ,  [0;0;0] );
assertEqual(  t.no_snapshots  , 2  );
