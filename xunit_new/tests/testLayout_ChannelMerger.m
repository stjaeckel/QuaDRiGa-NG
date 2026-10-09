function testLayout_ChannelMerger
%%
a = qd_arrayant('dipole');

l = qd_layout;
l.simpar.center_frequency = 2.185e9;
l.simpar.show_progress_bars = 0;
l.simpar.use_absolute_delays = 1;
l.tx_position = [0;0;125];
l.tx_array = a;
l.rx_array = a;

% Linear track with one segment
t = qd_track('linear',10,0);
t.interpolate_positions(2);
t.scenario = 'LOSonly';
t.segment_index = [1,10];
t.initial_position = [100,-100,0]';
l.rx_track = t;

% Test 1 : LOSonly - WINNER_UMa_C2_LOS
t.scenario = {'LOSonly','WINNER_UMa_C2_LOS'};
cb = l.init_builder;
gen_parameters( cb );
c = get_channels( cb );
d = merge(c,0.2);

% All LOS coeffs should be > 0
assertTrue( all( abs( reshape( d.coeff(1,1,1,:) , [],1) ) > 0 ) );

% All NLOS coeffs at the beginning sould be 0
assertTrue( all( reshape( d.coeff(1,1,2:end,1) , [],1) == 0 ) );

% All NLOS coeffs at the end sould be > 0
assertTrue( all( abs(reshape( d.coeff(1,1,2:end,end) , [],1)) >= 1e-10 ) );

% Test 2: WINNER_UMa_C2_LOS - LOSonly
t.scenario = {'WINNER_UMa_C2_LOS','LOSonly'};
cb = l.init_builder;
gen_parameters( cb );
c = get_channels( cb );
d = merge(c,0.8);

% All LOS coeffs should be > 0
assertTrue( all( abs( reshape( d.coeff(1,1,1,:) , [],1) ) > 0 ) );

% All NLOS coeffs at the beginning sould be > 0
assertTrue( all( abs(reshape( d.coeff(1,1,2:end,1) , [],1)) >= 1e-10 ) );

% All NLOS coeffs at the end sould be 0
assertTrue( all( reshape( d.coeff(1,1,2:end,end) , [],1) == 0 ) );


% Test 3 : LOSonly - LOSonly
% Create two tracks and then merge them
t.scenario = 'LOSonly';
cb = l.init_builder;
gen_parameters( cb );
c = get_channels( cb );
d = merge(c);

% Remove second segment from track and assign single track to CB
% Then run the calcualtion again
t.no_segments = 1;
cb(1,1).rx_track(1,1) = t;
c = cb.get_channels;

% "c(1)" and "d" should be identical 
diff_coeff = reshape( c(1,1).coeff - d.coeff , [],1);
assertTrue( all( abs( diff_coeff ) < 1e-14 ) );

diff_delay = reshape( c(1,1).delay - d.delay , [],1);
assertTrue( all( abs( diff_delay ) < 1e-14 ) );
