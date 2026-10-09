function testLayout_GetChannelsSeg_Baseline
%% Time-based channel generation with closed track and interpolation

l = qd_layout;
l.simpar.use_3GPP_baseline = 1;
l.simpar.center_frequency = 0.9*1e9;
l.simpar.sample_density = 1.1;
l.simpar.show_progress_bars = 0;
l.no_tx = 2;
l.tx_position = [10,30;10,10;20,20];
l.no_rx = 3;
l.randomize_rx_positions(100,1,2,1);
l.set_scenario('3GPP_38.901_InF_NLOS_DH');
l.rx_track(1,2).scenario{2,1} = '3GPP_38.901_InF_LOS';

% Generate channels using segment-by-segment method
ca1 = l.get_channels_seg(1,1);
ca2 = l.get_channels_seg(2,2);

% Generate all channels
[c9,b] = l.get_channels;

% Generate channels using segment-by-segment method
cb1 = l.get_channels_seg(1,1);
cb2 = l.get_channels_seg(2,2);
cb3 = l.get_channels_seg(1,3);

% Compare names
assertEqual( c9(1,1).name, ca1.name );
assertEqual( c9(2,2).name, ca2.name );

% Compare positions
assertTrue( all( c9(1,1).rx_position(:) - ca1.rx_position(:) < 1e-13 ) );

% Compare coefficients
fr91  = permute( fr( c9(1,1),    100e6, 5 ),[3,4,1,2] );
fr92  = permute( fr( c9(2,2),    100e6, 5 ),[3,4,1,2] );
fr93  = permute( fr( c9(3,1),    100e6, 5 ),[3,4,1,2] );
fra1  = permute( fr( ca1,   100e6, 5 ),[3,4,1,2] );
fra2  = permute( fr( ca2,   100e6, 5 ),[3,4,1,2] );
frb1  = permute( fr( cb1,   100e6, 5 ),[3,4,1,2] );
frb2  = permute( fr( cb2,   100e6, 5 ),[3,4,1,2] );
frb3  = permute( fr( cb3,   100e6, 5 ),[3,4,1,2] );

assertTrue( all( abs(fr91(:) - fra1(:) ) < 1e-13 ) );
assertTrue( all( abs(fr91(:) - frb1(:) ) < 1e-13 ) );
assertTrue( all( abs(fr92(:) - fra2(:) ) < 1e-13 ) );
assertTrue( all( abs(fr92(:) - frb2(:) ) < 1e-13 ) );
assertTrue( all( abs(fr93(:) - frb3(:) ) < 1e-13 ) );

% Check if absTOA_offset is initialized
assertTrue( ~isempty( b(1,1).absTOA_offset ) )

% Change baseline setting --> should triggier parameter reset
l.simpar(1,1).use_3GPP_baseline = 0;

cc1  = l.get_channels_seg(1,1);
frc1 = permute( fr( cc1,   100e6, 5 ),[3,4,1,2] );
assertFalse( all( abs(fr91(:) - frc1(:) ) < 1e-13 ) );

