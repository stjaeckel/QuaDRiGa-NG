function testTrack_interpoalte_snap_static
%%

t = qd_track( 'linear', 0 , 0);  
t.positions = rand(3,1);
t.orientation = rand(3,1);

t.set_speed(2);
assertEqual( t.movement_profile,[0 0.5;1 1] )

t.movement_profile = [0 0.2 0.5 ; 1 1 1];

t.scenario = {'A';'B'};

assertFalse( t.closed );

[ dist,ti ] = t.interpolate( 'snapshot', 0.1 ); 

o_snap = ones( 1,ti.no_snapshots );
assertEqual( ti.no_snapshots, 6 );
assertEqual( ti.initial_position, t.initial_position );
assertEqual( ti.positions, t.positions(:,o_snap) );
assertEqual( ti.orientation, t.orientation(:,o_snap) );
assertEqual( ti.no_segments, 1 );
assertEqual( ti.segment_index, 1 );
assertEqual( ti.scenario, t.scenario );
assertElementsAlmostEqual( ti.movement_profile(2,:), [1,3,6],  'absolute', 1e-13 );

assertFalse( ti.closed );
