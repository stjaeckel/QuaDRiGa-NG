function testBuilder_GR_baseline

b = qd_builder('3GPP_38.901_UMi_LOS_GR');
b.simpar.use_3GPP_baseline = 1;
b.simpar.show_progress_bars = 0;
b.tx_position = [0;0;3];
b.rx_positions = [20 21 ; 0 5;  3 3];

b.init_sos;
b.gen_lsf_parameters;

assertTrue( ~isempty( b.gr_epsilon_r ) );

b.gen_ssf_parameters;

assertEqual( b.check_los, 2 );

c = b.get_channels;

