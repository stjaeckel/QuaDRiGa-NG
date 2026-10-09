function testLayout_GetChannelsSeg_StaticRx
%% Time-based channel generation with dual mobility, interpolation, and multi-frequency

l = qd_layout;
l.simpar.center_frequency = [0.5,0.4]*1e9;
l.simpar.sample_density = 2.1;
l.simpar.show_progress_bars = 0;

l.tx_track = qd_track('linear',5,0);
l.tx_track.set_speed(0.5);
l.tx_position = [0,0,1]';
l.tx_track.name = 'Tx';

l.rx_track = qd_track('linear',0,0);
l.rx_track.positions = zeros(3,5);
l.rx_track.orientation(3,:) = (0:4)*pi/2;
l.rx_track.movement_profile = [0,10;1,5];
l.rx_track.segment_index = [1,3];
l.rx_position = [-10;10;1.5];
l.rx_track.name = 'Rx';
l.rx_track.scenario = {'Freespace','TwoRayGR'};

l.rx_array = qd_arrayant('patch');

l.update_rate = 0.02;

c9 = l.get_channels;

assertEqual( c9(1,1,1).no_path, 2 );
assertEqual( c9(1,1,2).no_path, 2 );

assertEqual( c9(1,1).center_frequency, l.simpar.center_frequency(1) );
assertEqual( c9(1,2).center_frequency, l.simpar.center_frequency(2) );

% The update rate is higher than the sample density, so channel interpolation should be used
assertTrue( l.use_channel_interpolation );

% The tx-track must have the same number of snapshots
[len,dist] = get_length( l.tx_track );
assertTrue( abs(len - 5) < 0.01 );
assertTrue( abs( 1/mean(diff(dist)) - l.simpar.samples_per_meter ) < 0.1 );

% The rx_track is interpolated to match the sample density
[len,dist] = get_length( l.rx_track );
assertTrue( abs(len) < 0.01 );

ca1     = l.get_channels_seg(1,1,1);
ca2     = l.get_channels_seg(1,1,2);
ca12    = l.get_channels_seg(1,1,[1,2],2);

assertEqual( ca1(1,1).center_frequency, l.simpar(1,1).center_frequency(1) );
assertEqual( ca1(1,2).center_frequency, l.simpar(1,1).center_frequency(2) );
assertEqual( ca12.center_frequency, l.simpar(1,1).center_frequency(2) );

% Compare number of snapshots
assertEqual( c9(1,1).no_snap, ca1(1,1).no_snap + ca2(1,1).no_snap );
assertEqual( c9(1,2).no_snap, ca1(1,2).no_snap + ca2(1,1).no_snap );
assertEqual( c9(1,2).no_snap, ca12.no_snap );

% Compare RX positions
assertTrue( all(all( abs( c9(1,1).rx_position - cat(2, ca1(1,1).rx_position, ca2(1,1).rx_position ) ) < 1e-13 )) );
assertTrue( all(all( abs( c9(1,1).rx_position - ca12(1,1).rx_position  ) < 1e-13 )) );

% Compare TX positions
assertEqual( size(c9(1,1).tx_position,2), c9(1,1).no_snap );
assertTrue( all(all( abs( c9(1,1).tx_position - cat(2, ca1(1,1).tx_position, ca2(1,1).tx_position ) ) < 1e-13 )) );
assertTrue( all(all( abs( c9(1,2).tx_position - ca12(1,1).tx_position  ) < 1e-13 )) );

% Compare coefficients
fr91  = permute( fr( c9(1,1),    100e6, 5 ),[3,4,1,2] );
fr92  = permute( fr( c9(1,2),    100e6, 5 ),[3,4,1,2] );
fr1   = permute( fr( ca1(1,1),   100e6, 5 ),[3,4,1,2] );
fr2   = permute( fr( ca2(1,1),   100e6, 5 ),[3,4,1,2] );
fr12  = permute( fr( ca12,  100e6, 5 ),[3,4,1,2] );

assertTrue( all(all(  abs(fr91 - cat(2, fr1,fr2 )) < 1e-13 )) );
assertTrue( all(all(  abs(fr92 - fr12) < 1e-13 )) );

assertEqual( ca12(1,1).par.update_rate, l.update_rate );
