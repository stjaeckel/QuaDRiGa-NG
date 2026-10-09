function testLayout_Constructor
%%
l = qd_layout;
assertEqual(  l.name  ,  'Layout' );
assertEqual(  l.no_tx  ,  1 );
assertEqual(  l.no_rx  ,  1 );
assertEqual(  l.tx_name  ,  {'Tx0001'} );
assertEqual(  l.tx_position  ,  [0;0;25] );
assertEqual(  l.rx_name  ,  {'Rx0001'});
assertEqual(  l.rx_position  ,  [0;0;0] );
assertEqual(  l.pairing  ,  [1;1] );
assertEqual(  l.no_links  ,  1 );

assertTrue( isa(l.tx_array , 'qd_arrayant')  );
assertTrue( isa(l.rx_array , 'qd_arrayant')  );
assertTrue( isa(l.rx_track , 'qd_track')  );
assertTrue( isa(l.tx_track , 'qd_track')  );

assertEqual(  l.rx_track.no_snapshots  ,  1 );

assertEqual(  l.tx_track.no_snapshots  ,  1 );

assertEqual(  numel(l.tx_array)  ,  1 );
assertEqual(  l.tx_array.no_elements  ,  1 );
assertEqual(  l.tx_array.no_az  ,  361 );
assertEqual(  l.tx_array.no_el  ,  181 );

assertEqual(  numel(l.rx_array)  ,  1 );
assertEqual(  l.rx_array.no_elements  ,  1 );
assertEqual(  l.rx_array.no_az  ,  361 );
assertEqual(  l.rx_array.no_el  ,  181 );

