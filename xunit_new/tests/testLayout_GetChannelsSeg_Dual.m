function testLayout_GetChannelsSeg_Dual
%% Time-based channel generation with dual mobility and interpolation

l = qd_layout;
l.simpar.center_frequency = 0.5*1e9;
l.simpar.sample_density = 2.1;
l.simpar.show_progress_bars = 0;

l.rx_track = qd_track('linear',10,pi/2);
l.rx_track.set_speed(1);
l.rx_position = [-10;10;1.5];
l.rx_track.name = 'Rx';
l.rx_track.add_segment([-10;16;1.5]);

l.tx_track = qd_track('linear',5,0);
l.tx_track.set_speed(0.5);
l.tx_position = [0,0,1]';
l.tx_track.name = 'Tx';

l.set_scenario('3GPP_38.901_UMa_LOS');

l.update_rate = 0.02;

c9 = l.get_channels;

% The update rate is higher than the sample density, so channel interpolation should be used
assertTrue( l.use_channel_interpolation );

% The rx_track is interpolated to match the sample density
[len,dist] = get_length( l.rx_track );
assertTrue( abs(len - 10) < 0.01 );
assertTrue( abs( 1/mean(diff(dist)) - l.simpar.samples_per_meter ) < 0.1 );

% The (shorter) tx-track must have the same number of snapshots
[len,dist] = get_length( l.tx_track );
assertTrue( abs(len - 5) < 0.01 );
assertTrue( abs( 1/mean(diff(dist)) - 2*l.simpar.samples_per_meter ) < 0.1 );

ca1     = l.get_channels_seg(1,1,1);
ca2     = l.get_channels_seg(1,1,2);
ca12    = l.get_channels_seg(1,1,[1,2]);

% Compare number of snapshots
assertEqual( c9.no_snap, ca1.no_snap + ca2.no_snap )
assertEqual( c9.no_snap, ca12.no_snap )

% Compare RX positions
assertTrue( all(all( abs( c9.rx_position - cat(2, ca1.rx_position, ca2.rx_position ) ) < 1e-13 )) )
assertTrue( all(all( abs( c9.rx_position - ca12.rx_position  ) < 1e-13 )) )

% Compare TX positions
assertEqual( size(c9.tx_position,2), c9.no_snap )
assertTrue( all(all( abs( c9.tx_position - cat(2, ca1.tx_position, ca2.tx_position ) ) < 1e-13 )) )
assertTrue( all(all( abs( c9.tx_position - ca12.tx_position  ) < 1e-13 )) )

% Compare coefficients
fr9   = permute( fr( c9,    100e6, 5 ),[3,4,1,2] );
fr1   = permute( fr( ca1,   100e6, 5 ),[3,4,1,2] );
fr2   = permute( fr( ca2,   100e6, 5 ),[3,4,1,2] );
fr12  = permute( fr( ca12,  100e6, 5 ),[3,4,1,2] );

assertTrue( all(all(  abs(fr9 - cat(2, fr1,fr2 )) < 1e-8 )) );
assertTrue( all(all(  abs(fr9 - fr12) < 1e-8 )) );

assertEqual( ca12.par.update_rate, l.update_rate )
