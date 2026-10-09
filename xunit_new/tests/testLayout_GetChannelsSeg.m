function testLayout_GetChannelsSeg

%% Position-based channel generation

l = qd_layout;
l.simpar.samples_per_meter = 3;
l.simpar.show_progress_bars = 0;
l.no_tx = 3;
l.no_rx = 2;

l.tx_track(1,2) = qd_track('linear',5);
l.tx_track(1,2).name = 'mobile-tx';
l.tx_track(1,2).initial_position = [5,5,1.5]';

l.tx_position(1,1) = 20;

l.randomize_rx_positions(100,1.5,1.5,0);
l.rx_track(1,2) = qd_track('linear',5,0);
l.rx_track(1,2).name = 'mobile-rx';
l.rx_track(1,2).initial_position = [-20,5,1.5]';

add_segment( l.rx_track, [-17,5,1.5]' );

l.set_scenario('3GPP_38.901_UMa_LOS',[],[1,3]);
l.set_scenario('QuaDRiGa_UD2D_LOS',[],2);
l.rx_track(1,2).scenario{2,2} = 'QuaDRiGa_UD2D_NLOS';
l.rx_track(1,2).scenario{3,1} = '3GPP_38.901_UMa_NLOS';

% Check if the number of snapshots is set correctly
assertEqual( l.tx_track(1,2).no_snapshots, 2 );
assertEqual( l.rx_track(1,2).no_snapshots, 3 );

% Calculate track checksum
hsh = checksum( l.rx_track ) + checksum( l.tx_track ) + checksum( l.simpar );

% Get all channels
c1 = l.get_channels; 

% Check if tracks didnt change
assertEqual( l.track_checksum, hsh );

% Builders should be initialized
assertTrue( ~isempty( l.h_qd_builder_init ) );

% No update-rate is used, so channel interpolation should be disabled
assertFalse( l.use_channel_interpolation );

% Test if automatic interpolation of snapshots is correct (only interpolate mobile channels)
assertEqual( [ c1.no_snap ], [1 16 16 16 1 16] );

% Creating the channels should not change the snapshot definitions
assertEqual( l.tx_track(1,2).no_snapshots, 2 );
assertEqual( l.rx_track(1,2).no_snapshots, 3 );

% Call get_channels a second time should reuse the exisiting builders and greate identical channels
c2 = l.get_channels;

% Check if tracks didnt change
assertEqual( l.track_checksum, hsh );

% Compare c2 and c2
assertEqual( size(c1), size(c2) );
for ic = 1 : numel( c1 )
    [ i1,i2 ] = qf.qind2sub( size( c1 ), ic );
    assertTrue( all( abs( c1(i1,i2).coeff(:) - c2(i1,i2).coeff(:) ) < 1e-12 ) )
    assertTrue( all( abs( c1(i1,i2).delay(:) - c2(i1,i2).delay(:) ) < 1e-12 ) )
end

% Generate segment-by-segment channels
c3 = l.get_channels_seg( 2,2,1 ); 
c4 = l.get_channels_seg( 2,2,2 ); 
c5 = l.get_channels_seg( 3,1 ); 

% Check if tracks didnt change
assertEqual( l.track_checksum, hsh );

assertEqual( c3.name, 'mobile-tx_mobile-rx-S01' );
assertTrue( all(all(abs( c1(2,2).tx_position(:,1:c3.no_snap) - c3.tx_position ) < 1e-12)) );
assertTrue( all(all(abs( c2(2,2).rx_position(:,1:c3.no_snap) - c3.rx_position ) < 1e-12)) );
tmp = c1(2,2).coeff(:,:,1:c3.no_path,1:c3.no_snap) - c3.coeff;  % Identical order of paths
assertTrue( all( abs( tmp(:) ) < 1e-12 ) );

ii = c3.no_snap+1:c1(2,2).no_snap;
assertTrue( all(all(abs( c1(2,2).tx_position(:,ii) - c4.tx_position ) < 1e-12)) );
assertTrue( all(all(abs( c2(2,2).rx_position(:,ii) - c4.rx_position ) < 1e-12)) );
tmp1 = fr( c1(2,2), 10e6, 5, ii );    % Different order after merging
tmp2 = fr( c4, 10e6, 5 );
assertTrue( all( abs( tmp1(:)-tmp2(:) ) < 1e-8 ) );

assertEqual( c5.name, 'Tx0003_Rx0001' );
assertTrue( all( abs( c1(1,3).coeff(:) - c5.coeff(:) ) < 1e-12 ) );

%% Time-based channel generation without interpolation

% Update layout
l.no_tx = 2;
l.tx_track(1,2) = qd_track('linear',2.5);
l.tx_track(1,2).name = 'mobile-tx';
l.tx_track(1,2).initial_position = [5,5,1.5]';

l.rx_track(1,1).scenario = l.rx_track(1,1).scenario(1:2,:);
l.rx_track(1,2).scenario = l.rx_track(1,2).scenario(1:2,:);

l.tx_track(1,2).set_speed(5);
l.rx_track(1,2).set_speed(10);
l.update_rate = 0.05;
l.simpar(1,1).sample_density = 2.1;

l.simpar(1,1).center_frequency = [2,20]*1e9;

% Builders should NOT be initialized
assertTrue( isempty( l.h_qd_builder_init ) );
assertTrue( isempty( l.builder_index ) );

hsh = checksum( l.rx_track ) + checksum( l.tx_track ) + checksum( l.simpar );

% Generate some single channels
c81 = l.get_channels_seg( 1,2,[],[],0.8 );    	% All segments and frequencies
c82 = l.get_channels_seg( 1,2,1,[],0.8 );      	% First segment, all frequencies
c83 = l.get_channels_seg( 1,2,2,[],0.8 );      	% Second segment, all frequencies
c84 = l.get_channels_seg( 1,2,1,2,0.8 );        % First segment, first frequency
c85 = l.get_channels_seg( 1,2,2,1,0.8 );        % Second segment, second frequency

% Check if tracks didnt change
assertEqual( l.track_checksum, hsh );

% Check channel names
assertEqual( c81(1,1).name ,'F01-Tx0001_mobile-rx' );
assertEqual( c81(1,2).name ,'F02-Tx0001_mobile-rx' );
assertEqual( c82(1,1).name ,'F01-Tx0001_mobile-rx-S01' );
assertEqual( c82(1,2).name ,'F02-Tx0001_mobile-rx-S01' );
assertEqual( c83(1,1).name ,'F01-Tx0001_mobile-rx-S02' );
assertEqual( c83(1,2).name ,'F02-Tx0001_mobile-rx-S02' );
assertEqual( c84(1,1).name ,'F02-Tx0001_mobile-rx-S01' );
assertEqual( c85(1,1).name ,'F01-Tx0001_mobile-rx-S02' );

% The update rate is too low, so channel interpolation should be disabled
assertFalse( l.use_channel_interpolation );

% Get all channels
c6 = l.get_channels(0,0.8);
assertEqual( size(c6),[2 2 2] );

% Test if automatic interpolation of snapshots is correct
assertEqual( [ c6.no_snap ], ones(1,8)*11 );

% Compare coefficients c6 vs c8
assertTrue( all( abs( c6(2,1,2).coeff(:) - c81(1,2).coeff(:) ) < 1e-12 ) )
assertTrue( all( abs( c6(2,1,1).coeff(:) - c81(1,1).coeff(:) ) < 1e-12 ) )

tmp = c6(2,1,2).coeff(:,:,1:c82(1,1).no_path,1:c82(1,1).no_snap) - c82(1,2).coeff;  % Identical order of paths
assertTrue( all( abs( tmp(:) ) < 1e-12 ) );

tmp = c6(2,1,2).coeff(:,:,1:c82(1,1).no_path,1:c82(1,1).no_snap) - c84(1,1).coeff;  % Identical order of paths
assertTrue( all( abs( tmp(:) ) < 1e-12 ) );

ii = c82(1,1).no_snap+1 : c6(2,1,1).no_snap;
tmp1 = fr( c6(2,1,2), 10e6, 5, ii );    % Different order after merging
tmp2 = fr( c83(1,2), 10e6, 5 );
assertTrue( all( abs( tmp1(:)-tmp2(:) ) < 1e-8 ) );

tmp1 = fr( c6(2,1,1), 10e6, 5, ii );    % Different order after merging
tmp2 = fr( c85, 10e6, 5 );
assertTrue( all( abs( tmp1(:)-tmp2(:) ) < 1e-8 ) );

% Get the static channel
c7 = l.get_channels_seg( 1,1 ); 
assertEqual( c7(1,2).name ,'F02-Tx0001_Rx0001'  )
assertEqual( c7(1,2).no_snap , 11 )
assertTrue( all( abs( c6(1,1,2).coeff(:) - c7(1,2).coeff(:) ) < 1e-12 ) )
assertTrue( all( abs( c6(1,1,1).coeff(:) - c7(1,1).coeff(:) ) < 1e-12 ) )
