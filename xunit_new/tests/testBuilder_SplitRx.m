function testBuilder_SplitRx

%% Single-Mobility

b = qd_builder('3GPP_38.901_UMi_LOS_GR');
b.simpar.show_progress_bars = 0;
b.name = 'scenario_txname';
b.simpar.center_frequency = [2,20]*1e9;
b.tx_position = [0;0;10];
b.rx_track = qd_track('linear',1);
b.rx_track.initial_position = [10;10;1.5];
b.rx_track(1,2) = qd_track('linear',1);
b.rx_track(1,2).initial_position = [10;-10;1.5];
b.rx_track(1,2).name = 'track2';

b.check_dual_mobility;
assertFalse( b.dual_mobility );

b.init_sos;

% Option 1
b1 = split_rx( b );
gen_parameters( b1 );
b1 = split_multi_freq( b1 );
c = get_channels( b1 );

% Option 2
gen_parameters( b );
b2 = split_multi_freq( b );
b2 = split_rx( b2 );
c(2,:) = get_channels( b2 );

% Option 3
b3 = split_rx( b );
b3 = split_multi_freq( b3 );
c(3,:) = get_channels( b3 );

% Option 4
b4 = split_multi_freq( b );
c(4,:) = get_channels( b4 );

for n = 2:4
    for m = 1:4
        assertTrue( all( abs( c(1,m).coeff(:) - c(n,m).coeff(:) ) < 1e-12 ) )
        assertTrue( all( abs( c(1,m).delay(:) - c(n,m).delay(:) ) < 1e-12 ) )
        assertTrue( all( abs( c(1,m).par.fbs_pos(:) - c(n,m).par.fbs_pos(:) ) < 1e-11 ) ) % mm-precision
        assertTrue( all( abs( c(1,m).par.lbs_pos(:) - c(n,m).par.lbs_pos(:) ) < 1e-11 ) ) % mm-precision
    end
end


%% Dual-Mobility

b = qd_builder('3GPP_38.901_InF_NLOS_DH');
b.simpar.show_progress_bars = 0;
b.name = 'scenario_txname';
b.simpar.center_frequency = [2,20]*1e9;

b.tx_track = qd_track('linear',1);
b.tx_track.initial_position = [0;0;1.5];

b.rx_track = qd_track('linear',1);
b.rx_track.initial_position = [10;10;1.5];
b.rx_track(1,2) = qd_track('linear',1);
b.rx_track(1,2).initial_position = [10;-10;1.5];
b.rx_track(1,2).name = 'track2';

b.check_dual_mobility;
assertTrue( b.dual_mobility );

b.init_sos;

% Option 1
b1 = split_rx( b );
gen_parameters( b1 );
b1 = split_multi_freq( b1 );
c = get_channels( b1 );

% Option 2
gen_parameters( b );
b2 = split_multi_freq( b );
b2 = split_rx( b2 );
c(2,:) = get_channels( b2 );

% Compare LSF parameters
assertTrue( all( abs( cat(1,b1.ds) - cat(1,b2.ds) ) < 1e-12 ) )
assertTrue( all( abs( cat(1,b1.sf) - cat(1,b2.sf) ) < 1e-12 ) )
assertTrue( all( abs( cat(1,b1.asD) - cat(1,b2.asD) ) < 1e-12 ) )
assertTrue( all( abs( cat(1,b1.asA) - cat(1,b2.asA) ) < 1e-12 ) )
assertTrue( all( abs( cat(1,b1.esD) - cat(1,b2.esD) ) < 1e-12 ) )
assertTrue( all( abs( cat(1,b1.esA) - cat(1,b2.esA) ) < 1e-12 ) )

% Compare SSF Parameters
assertTrue( all(all( abs( cat(1,b1.taus) - cat(1,b2.taus) ) < 1e-12 )) )
assertTrue( all(all( abs( cat(1,b1.pow) - cat(1,b2.pow) ) < 1e-12 )) )
assertTrue( all(all( abs( cat(1,b1.AoD) - cat(1,b2.AoD) ) < 1e-12 )) )
assertTrue( all(all( abs( cat(1,b1.AoA) - cat(1,b2.AoA) ) < 1e-12 )) )
assertTrue( all(all( abs( cat(1,b1.EoD) - cat(1,b2.EoD) ) < 1e-12 )) )
assertTrue( all(all( abs( cat(1,b1.EoA) - cat(1,b2.EoA) ) < 1e-12 )) )

% Option 3
b3 = split_rx( b );
b3 = split_multi_freq( b3 );
c(3,:) = get_channels( b3 );

% Option 4
b4 = split_multi_freq( b );
c(4,:) = get_channels( b4 );

for n = 2:4
    for m = 1:4
        assertTrue( all( abs( c(1,m).coeff(:) - c(n,m).coeff(:) ) < 1e-12 ) )
        assertTrue( all( abs( c(1,m).delay(:) - c(n,m).delay(:) ) < 1e-12 ) )
        assertTrue( all( abs( c(1,m).par.fbs_pos(:) - c(n,m).par.fbs_pos(:) ) < 1e-11 ) ) % mm-precision
        assertTrue( all( abs( c(1,m).par.lbs_pos(:) - c(n,m).par.lbs_pos(:) ) < 1e-11 ) ) % mm-precision
    end
end

