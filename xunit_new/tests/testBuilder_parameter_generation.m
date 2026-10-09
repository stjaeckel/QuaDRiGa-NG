function testBuilder_parameter_generation

b = qd_builder('3GPP_38.901_UMi_LOS_GR');
b.name = 'b1_t1';
b.simpar.show_progress_bars = 0;
b.simpar.center_frequency = [2,20,30]*1e9;
b.tx_position = [0;0;10];
b.rx_track = qd_track('linear',1);
b.rx_track.initial_position = [10;10;1.5];
b.rx_track.name = 'r1';
b.rx_track(1,2) = qd_track('linear',1);
b.rx_track(1,2).initial_position = [10;-10;1.5];
b.rx_track(1,2).name = 'r2';

b(1,2) = qd_builder('3GPP_38.901_InF_NLOS_DH');
b(1,2).name = 'b2_t2';
b(1,2).simpar(1,1).show_progress_bars = 0;
b(1,2).simpar(1,1).center_frequency = [2,20,30]*1e9;
b(1,2).tx_position = [0;50;25];
b(1,2).rx_positions = [10;10;1.5];

b(1,3) = qd_builder('3GPP_38.901_UMa_NLOS');

b(1,4) = qd_builder('mmMAGIC_UMi_LOS');
b(1,4).name = 'b4_t4';
b(1,4).simpar(1,1).show_progress_bars = 0;
b(1,4).simpar(1,1).use_3GPP_baseline = 1;
b(1,4).simpar(1,1).center_frequency = 14*1e9;
b(1,4).tx_position = [0;50;25];
b(1,4).rx_positions = [10 20 ;10 20 ;1.5 1.4];


%% Init SOS
init_sos( b );

assertFalse( isempty( b(1,1).sos ) );
assertFalse( isempty( b(1,1).gr_sos ) );
assertTrue( isempty( b(1,1).absTOA_sos ) );
assertFalse( isempty( b(1,2).sos ) );
assertFalse( isempty( b(1,2).absTOA_sos ) );
assertFalse( isempty( b(1,3).sos ) );
assertTrue( isempty( b(1,4).absTOA_sos ) );
assertTrue( isempty( b(1,4).path_sos ) );

ph = b(1,1).sos(1,1).sos_phase;

init_sos(b,2);
assertEqual(  ph, b(1,1).sos(1,1).sos_phase );

init_sos(b(1,1));
assertTrue( all(abs(ph(:) - b(1,1).sos(1,1).sos_phase(:) ) > 1e-5) );

init_sos(b(1,3),0);
assertTrue( isempty( b(1,3).sos ) );


%% LSF Parameters
gen_lsf_parameters( b );

assertEqual( size( b(1,1).ds ), [b(1,1).no_freq,b(1,1).no_rx_positions] );
assertTrue( isempty( b(1,1).absTOA_offset ) );
assertFalse( isempty( b(1,1).gr_epsilon_r ) );
assertFalse( isempty( b(1,2).absTOA_offset ) );
assertTrue( isempty( b(1,2).gr_epsilon_r ) );

ds = b(1,1).ds;
gr_epsilon_r = b(1,1).gr_epsilon_r;
absTOA_offset = b(1,2).absTOA_offset;

gen_lsf_parameters( b, 2 );

assertEqual( ds , b(1,1).ds );
assertEqual( gr_epsilon_r , b(1,1).gr_epsilon_r );
assertEqual( absTOA_offset , b(1,2).absTOA_offset );

gen_lsf_parameters( b,0 ); % Same SOS - Same results

assertTrue( isempty( b(1,1).ds ) );
assertTrue( isempty( b(1,1).gr_epsilon_r ) );
assertTrue( isempty( b(1,2).absTOA_offset ) );

gen_lsf_parameters( b,[],0 ); % Same SOS - Same results

assertEqual( ds , b(1,1).ds );
assertEqual( gr_epsilon_r , b(1,1).gr_epsilon_r );
assertEqual( absTOA_offset , b(1,2).absTOA_offset );

%% SSF Parameters
gen_ssf_parameters( b );

gr_epsilon_r = b(1,1).gr_epsilon_r;

kf = b(1,1).kf;
sf = b(1,1).sf;


%% Estimation of LSF parameters from SSF parameters

gen_lsf_parameters(b,0);
gen_lsf_from_ssf( b ,0 );

assertTrue( all( abs( gr_epsilon_r(:) - b(1,1).gr_epsilon_r(:) ) < 1e-12 ) );   % GR Reflectivity is correct?
assertTrue( all( abs( kf(:) - b(1,1).kf(:) ) < 1e-12 ) );
assertTrue( all( abs( sf(:) - b(1,1).sf(:) ) < 1e-12 ) );

gen_lsf_from_ssf( b ,1 );

assertTrue( all( abs( gr_epsilon_r(:) - b(1,1).gr_epsilon_r(:) ) < 1e-12 ) );   % GR Reflectivity is correct?
assertTrue( all( abs( kf(:) - b(1,1).kf(:) ) < 1e-12 ) );
assertTrue( all( abs( sf(:) - b(1,1).sf(:) ) < 1e-12 ) );

%% FBS and LBS positions
gen_fbs_lbs( b );

taus1 = b(1,1).taus;
taus2 = b(1,2).taus;

gen_ssf_from_scatterers( b );

assertTrue( all(abs(taus1(:) - b(1,1).taus(:)) < 1e-16) )
assertTrue( all(abs(taus2(:) - b(1,2).taus(:)) < 1e-16) )

