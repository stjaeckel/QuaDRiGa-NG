function testTrack_Subtrack_spatial_cconsistency
%%

t = qd_track( 'linear' , 5, pi/2 );    
t.initial_position = [20 ; 30 ; 1.5 ];    
t.interpolate_positions(10);
t.no_segments = t.no_snapshots;
t.scenario = '3GPP_38.901_UMi_LOS';

pos = t.positions_abs;

s = get_subtrack( t );

assertEqual( size(s,2), t.no_snapshots );
assertEqual( size(s,2), t.no_segments );

for n = 1:t.no_snapshots
    assertEqual( pos(:,n),s(1,n).initial_position );
end