function testLayout_Builder_Assembly
%%

l = qd_layout;
l.simpar.show_progress_bars = false;
l.no_tx = 3;
l.no_rx = 2;

l.randomize_rx_positions(100,1,1,0);

t = qd_track.generate('linear',10);
t.interpolate_positions(1);
t.initial_position = l.rx_position(:,2);
t.segment_index = [1 5];
l.rx_track(1,2) = t;

l.tx_array(1,2) = qd_arrayant('xpol');
l.rx_array(1,2) = qd_arrayant('xpol');

l.set_scenario( 'LOSonly' );

c = l.get_channels;

assertEqual( size(c),[2,3] );

assertEqual( cat(2,c.no_rxant),[1,2,1,2,1,2] );
assertEqual( cat(2,c.no_txant),[1,1,2,2,1,1] );

