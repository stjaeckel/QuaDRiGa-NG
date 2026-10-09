function testBuilder_get_angles_distances
%%
b = qd_builder('3GPP_3D_UMa_LOS');
b.rx_positions = [0,10,10;50,0,0]';
b.tx_position = [0;0;0];
ang = b.get_angles;
assertEqual( ang, [90,0 ; -90,-180 ; 45 0 ; -45 0] );

dist = b.get_distances;
assertEqual( dist, [10 50] );