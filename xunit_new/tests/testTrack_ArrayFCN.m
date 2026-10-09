function testTrack_ArrayFCN
%%
a = qd_track('linear',200,1);
a.interpolate_positions(2);
a(1,1,2) = qd_track('circular');
a(1,1,3) = a(1,1,1);

set_scenario( a(1,1,1:2), 'Test', 1, 10, 10, 0 )

assertEqual( a(1,1,1).segment_index, a(1,1,3).segment_index  );

tst = true(1,1,3);
tst(1,1,2) = false;
assertEqual( qf.eqo( a(1,1,1), a )  , tst  );

calc_orientation( a );

si = a(1,1,1).segment_index;
[~,di] = get_length( a(1,1,1) );
di=di(si);

correct_overlap( a , 0.2 );

si = a(1,1,3).segment_index;
[~,din] = get_length( a(1,1,3) );
din=din(si);

x = din(2:end) - di(2:end);
x_max = 10.5*0.2*0.67;
x_min = (10-10*0.2*0.66)*0.2*0.66;

assertTrue( all( x > x_min & x < x_max )  );

[len, dst] = get_length( a );

assertEqual( len(1), len(3)  );
assertEqual( dst{1}, dst{3}  );

set_speed( a, 10 )

assertFalse( isempty( a(1,1,1).movement_profile ) );
assertEqual( a(1,1,1).movement_profile, a(1,1,3).movement_profile  );

split_segment( a , 3.9 , 8 , 5 , 0  )
assertEqual( a(1,1,1).segment_index, a(1,1,3).segment_index  );

assertEqual( a(1,1,1).no_segments, 40  );

