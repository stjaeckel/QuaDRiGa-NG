function testBuilder_Visualize_Clusters
%%

b = qd_builder('3GPP_38.901_UMi_LOS');
b.simpar.show_progress_bars = 0;

b.tx_position  = [0;0;25];
b.rx_positions = [100,0,1.5]';
gen_parameters(b);

visualize_clusters(b);
visualize_clusters(b,1,3);
close all
