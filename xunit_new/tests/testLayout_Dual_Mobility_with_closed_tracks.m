function testLayout_Dual_Mobility_with_closed_tracks
%%

l = qd_layout;
l.simpar.show_progress_bars = 0;
l.simpar.samples_per_meter = 1;

% Tx open, Rx closed

l.tx_track = qd_track( 'linear',9.9 );
l.tx_track(1,1).interpolate_positions( 1 );
l.tx_track(1,1).initial_position(3) = 1.5;
l.rx_track = qd_track('circular',9.9);
l.rx_track(1,1).interpolate_positions( 1 );
l.set_scenario('LOSonly');

l.rx_track(1,1).segment_index(2) = 5;

assertTrue( ~l.tx_track(1,1).closed );
assertTrue( l.rx_track(1,1).closed );
assertTrue( l.dual_mobility );

c = l.get_channels;

assertEqual( c.no_snap, l.rx_track(1,1).no_snapshots )

% Both closed

l.tx_track = qd_track( 'circular',9.9 );
l.tx_track(1,1).interpolate_positions( 1 );
l.tx_track(1,1).initial_position(3) = 1.5;
l.rx_track(1,1).segment_index = 1;

assertTrue( l.tx_track(1,1).closed );
assertTrue( l.rx_track(1,1).closed );
assertTrue( l.dual_mobility );

c = l.get_channels;

assertEqual( c.no_snap, l.rx_track(1,1).no_snapshots-1 );

% Tx closed, Rx open

l.rx_track = qd_track( 'linear',9.9 );
l.rx_track(1,1).interpolate_positions( 1 );
l.set_scenario('LOSonly');

assertTrue( l.tx_track(1,1).closed );
assertTrue( ~l.rx_track(1,1).closed );
assertTrue( l.dual_mobility );

c = l.get_channels;

assertEqual( c.no_snap, l.rx_track(1,1).no_snapshots );

