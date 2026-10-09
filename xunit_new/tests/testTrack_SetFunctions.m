function testTrack_SetFunctions
%%
t = qd_track;
t.no_snapshots = 10;
assertEqual(  size(t.positions)  ,  [3 10] );
t.calc_orientation;
assertEqual(  size(t.orientation)  ,  [3 10] );
t.no_snapshots = 5;
assertEqual(  size(t.positions)  ,  [3 5] );
assertEqual(  size(t.orientation)  ,  [3 5] );

try 
    t.name = [];
    assertTrue( false );
catch err
    assertEqual( err.identifier   ,  'QuaDRiGa:qd_track:wrongInputValue' );
end

try 
    t.name = 'b_l_a';
    assertTrue( false );
catch err
    assertEqual( err.identifier   ,  'QuaDRiGa:qd_track:wrongInputValue' );
end
