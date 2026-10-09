function testBuilder_get_pl
%%

b = qd_builder('3GPP_3D_UMa_LOS');
b.tx_position = [0,0,25]';
b.rx_positions = [0,10,10;0,50,0;200,0,0]';
b.rx_track = qd_track.generate('linear',10);
b.rx_track.initial_position = b.rx_positions(:,1);
b.simpar.show_progress_bars = 0;

b.check_dual_mobility;

pl1 = b.get_pl;
assertEqual( size(pl1),[1,3] );

pl2 = b.get_pl( b.rx_track );
assertEqual( numel(pl2), b.rx_track(1,1).no_snapshots );

assertEqual( pl2(1), pl1(1) )

b.gen_parameters;

[ sf,kf ] = get_sf_profile( b, b.rx_track(1,1) );
assertEqual( size(sf), [1,b.rx_track(1,1).no_snapshots] );
assertEqual( size(kf), [1,b.rx_track(1,1).no_snapshots] );

assertTrue( abs(b.sf(1)-sf(1))<1e-5 )
assertTrue( abs(b.kf(1)-kf(1))<1e-5 )

