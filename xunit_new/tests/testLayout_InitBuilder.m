function testLayout_InitBuilder
%%
l = qd_layout;
l.no_tx = 2;
l.no_rx = 2;
l.randomize_rx_positions( 100, 1.5, 1.5, 20 );
l.rx_track(1,2).segment_index = [1 300];
l.set_scenario('Freespace');
l.rx_track(1,2).scenario{2,2} = 'LOSonly';

b = l.init_builder;

assertEqual( size(b) , [2 2] );
assertEqual( b(1,1).name , 'Freespace_Tx0001' );
assertEqual( b(1,2).name , 'Freespace_Tx0002' );
assertEqual( b(2,2).name , 'LOSonly_Tx0002' );

assertEqual( b(1,1).no_rx_positions , 3 );
assertEqual( b(1,2).no_rx_positions , 2 );
assertEqual( b(2,1).no_rx_positions , 0 );
assertEqual( b(2,2).no_rx_positions , 1 );

assertEqual( size( b(1,1).rx_track ) , [1 3] );
assertEqual( size( b(1,2).rx_track ) , [1 2] );
assertEqual( size( b(2,2).rx_track ) , [1 1] );

assertEqual(  b(1,1).rx_track(1,1).name  , 'Rx0001');
assertEqual(  b(1,1).rx_track(1,2).name  , 'Rx0002_seg0001');
assertEqual(  b(1,1).rx_track(1,3).name  , 'Rx0002_seg0002');
assertEqual(  b(1,2).rx_track(1,1).name  , 'Rx0001');
assertEqual(  b(1,1).rx_track(1,2).name  , 'Rx0002_seg0001');
assertEqual(  b(2,2).rx_track(1,1).name  , 'Rx0002_seg0002');

