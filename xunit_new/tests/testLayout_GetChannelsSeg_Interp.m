function testLayout_GetChannelsSeg_Interp
%% Time-based channel generation with closed track and interpolation

l = qd_layout;
l.simpar.center_frequency = 0.9*1e9;
l.simpar.sample_density = 1.1;
l.simpar.show_progress_bars = 0;
l.tx_position = [10;10;20];
l.rx_track = qd_track('circular',pi*2,0);
l.rx_track.initial_position(3,1) = 1.5;
l.rx_track.add_segment([-1,1,1.5]');
l.rx_track.add_segment([-1,-1,1.5]');
l.rx_track.set_speed(0.5);
l.set_scenario('3GPP_38.901_UMa_LOS');
l.update_rate = 0.2;

c9 = l.get_channels;

% The update rate is higher than the sample density, so channel interpolation should be used
assertTrue( l.use_channel_interpolation );

% The rx_track is interpolated to match the sample density
[len,dist] = get_length( l.rx_track );
assertTrue( abs(len - 2*pi) < 0.01 );
assertTrue( abs( 1/mean(diff(dist)) - l.simpar.samples_per_meter ) < 0.1 );

ca1     = l.get_channels_seg(1,1,1);
ca2     = l.get_channels_seg(1,1,2);
ca3     = l.get_channels_seg(1,1,3);
ca12    = l.get_channels_seg(1,1,[1,2]);
ca23    = l.get_channels_seg(1,1,[2,3]);
ca123   = l.get_channels_seg(1,1,[1,2,3]);

% Compare number of snapshots
assertEqual( c9.no_snap, ca1.no_snap + ca2.no_snap + ca3.no_snap )
assertEqual( c9.no_snap, ca12.no_snap + ca3.no_snap )
assertEqual( c9.no_snap, ca1.no_snap + ca23.no_snap )
assertEqual( c9.no_snap, ca123.no_snap )

% Compare RX positions
assertTrue( all(all( abs( c9.rx_position - cat(2, ca1.rx_position, ca2.rx_position, ca3.rx_position ) ) < 1e-13 )) )
assertTrue( all(all( abs( c9.rx_position - cat(2, ca12.rx_position, ca3.rx_position ) ) < 1e-13 )) )
assertTrue( all(all( abs( c9.rx_position - cat(2, ca1.rx_position, ca23.rx_position ) ) < 1e-13 )) )
assertTrue( all(all( abs( c9.rx_position - ca123.rx_position ) < 1e-13 )) )

% Compare coefficients
fr9   = permute( fr( c9,    100e6, 5 ),[3,4,1,2] );
fr1   = permute( fr( ca1,   100e6, 5 ),[3,4,1,2] );
fr2   = permute( fr( ca2,   100e6, 5 ),[3,4,1,2] );
fr3   = permute( fr( ca3,   100e6, 5 ),[3,4,1,2] );
fr12  = permute( fr( ca12,  100e6, 5 ),[3,4,1,2] );
fr23  = permute( fr( ca23,  100e6, 5 ),[3,4,1,2] );
fr123 = permute( fr( ca123, 100e6, 5 ),[3,4,1,2] );

assertTrue( all(all(  abs(fr9 - cat(2, fr1,fr2,fr3 )) < 1e-8 )) );
assertTrue( all(all(  abs(fr9 - cat(2, fr12,fr3 )) < 1e-8 )) );
assertTrue( all(all(  abs(fr9 - cat(2, fr1,fr23 )) < 1e-8 )) );
assertTrue( all(all(  abs(fr9 - fr123 ) < 1e-8 )) );

assertEqual( ca23.par.update_rate, l.update_rate )
