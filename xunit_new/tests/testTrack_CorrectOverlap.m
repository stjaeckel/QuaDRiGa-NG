function testTrack_CorrectOverlap
%%
t = qd_track.generate('street',300,pi/2);
t.set_scenario('Bla',1,100,120,5);
assertEqual( t.no_segments , 3 );

par.ds = [1 2 3]*1e-9;
par.pg = (0.1:0.1:30)+90;
t.par = par;

t.split_segment( 10,21,15,5 );

A = t.segment_index;

t.correct_overlap;

B = t.segment_index;

assertTrue( A(1) == 1 );
assertTrue( B(1) == 1 );

for n=2:t.no_segments
    assertTrue( A(n) < B(n) );
end

