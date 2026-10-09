function testLayout_Generator
%%
l = qd_layout.generate('random',5,400);
assertEqual(  l.no_tx  ,  5 );

l = qd_layout.generate('regular',[],[],qd_arrayant('patch'));
assertEqual(  l.no_tx  ,  7 ); % Default no tx
assertEqual(  l.tx_position(2,3)  , 500 ); % Default ISD
assertEqual(  l.tx_array(7).no_elements  , 3 ); % Default no sectors

l = qd_layout.generate('regular6',37,100);
assertEqual(  l.no_tx  ,  37 );
assertEqual(  l.tx_position(2,3)  , 100 ); 
assertEqual(  l.tx_position(2,11)  , 200 ); 
assertEqual(  l.tx_position(2,25)  , 300 ); 
assertEqual(  l.tx_array(37).no_elements  , 6 ); % 6 sectors

l = qd_layout.generate('regular',19,100);
assertEqual(  l.no_tx  ,  19 );
assertEqual(  l.tx_position(2,3)  , 100 ); % Default ISD

l.no_rx = 4;
assertEqual(  l.no_rx  ,  4 );
assertEqual(  numel(l.rx_name)  ,  4 );
assertEqual(  size(l.rx_position)  ,  [3,4] );
assertEqual(  numel(l.rx_array)  ,  4 );
assertEqual(  numel(l.rx_track)  ,  4 );
assertEqual(  l.no_links  ,  4*19 );
  
l.set_pairing;

l.no_tx = 5;
assertEqual(  size(l.tx_position)  ,  [3,5] );
assertEqual(  numel(l.tx_name)  ,  5 );
assertEqual(  numel(l.tx_array)  ,  5 );
assertEqual(  l.no_links  ,  20 );

l.no_rx = 1;
assertEqual(  l.no_rx  ,  1 );
assertEqual(  numel(l.rx_name)  ,  1);
assertEqual(  size(l.rx_position)  ,  [3,1] );
assertEqual(  numel(l.rx_array)  ,  1 );
assertEqual(  numel(l.rx_track)  ,  1 );
assertEqual(  l.no_links  ,  5 );
