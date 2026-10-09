function testLayout_set_pairing_identical_tx_rx
%%

set_rand_state( 1 );

l = qd_layout;
l.simpar.show_progress_bars = 1;
l.simpar.use_3GPP_baseline = 1;
l.no_tx = 4;
l.no_rx = 3;

l.randomize_rx_positions(100,1.5,1.5,0);
l.tx_position(:,1) = l.rx_position(:,1);
l.tx_position(:,2) = l.rx_position(:,2);
l.tx_position(:,3) = l.rx_position(:,3);

l.set_scenario('QuaDRiGa_UD2D',[],1:3);
l.set_scenario('3GPP_3D_UMi',[],4);

p = l.set_pairing([],[],[],[],0);
assertEqual( p, [1 1 2 2 3 3 4 4 4 ; 2 3 1 3 1 2 1 2 3] );

l.set_pairing('power',-999,[],[],0);
assertEqual( p, [1 1 2 2 3 3 4 4 4 ; 2 3 1 3 1 2 1 2 3] );


