function testBuilder_gen_lsf_parameters
%%

b = qd_builder('3GPP_3D_UMa_LOS');
b.simpar.show_progress_bars = 0;
b.rx_positions = [0,10,10;0,50,0]';
b.tx_position = [0;0;0];

b(1,2) = qd_builder('3GPP_3D_UMa_NLOS');
b(1,2).rx_positions = [0,10,10;0,50,0;200,0,0]';
b(1,2).tx_position = [0;0;0];

gen_parameters(b);

assertEqual( b(1,1).scenario,'3GPP_3D_UMa_LOS' );
assertEqual( b(1,2).scenario,'3GPP_3D_UMa_NLOS' );
assertEqual( b(1,1).no_rx_positions,2 );
assertEqual( b(1,2).no_rx_positions,3 );
assertEqual( numel(b(1,1).ds),2 );
assertEqual( numel(b(1,2).ds),3 );

tmp = b(1,1).ds;
b(1,1).gen_parameters;
assertEqual(tmp,b(1,1).ds );