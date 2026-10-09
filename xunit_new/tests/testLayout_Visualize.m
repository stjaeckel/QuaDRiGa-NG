function testLayout_Visualize
%%
l = qd_layout;
l.no_rx = 2;
l.no_tx = 2;
l.tx_position(2,2) = 200;
l.randomize_rx_positions(100,1.5,1.5,50);
l.set_scenario( 'Freespace' );
l.visualize;
close all
