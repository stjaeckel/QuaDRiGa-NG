function testBuilder_spatial_consistency
%%

set_rand_state(1);

b = qd_builder('3GPP_38.901_UMi_LOS');
b.simpar.use_absolute_delays;
b.name = 'Tx1';
b.simpar.show_progress_bars = false;
b.tx_position = [0;0;25];

% b.simpar.use_random_initial_phase = 1;
% b.plpar = [];
% b.scenpar.NumClusters = 4;
% b.scenpar.SF_sigma = 0;

b.rx_positions = [100,0,1.5 ; 100,0,1.5 ; 0,100,1.5]';

gen_parameters(b,[],0);

assertTrue(  all( abs( b.taus(1,:) - b.taus(2,:) ) < 1e-12 ) )
assertFalse( all( abs( b.taus(1,:) - b.taus(3,:) ) < 1e-12 ) )

assertTrue(  all( abs( b.pow(1,:) - b.pow(2,:) ) < 1e-12 ) )
assertFalse( all( abs( b.pow(1,:) - b.pow(3,:) ) < 1e-12 ) )

assertTrue(  all( abs( b.AoD(1,:) - b.AoD(2,:) ) < 1e-12 ) )
assertFalse( all( abs( b.AoD(1,:) - b.AoD(3,:) ) < 1e-12 ) )

assertTrue(  all( abs( b.AoA(1,:) - b.AoA(2,:) ) < 1e-12 ) )
assertFalse( all( abs( b.AoA(1,:) - b.AoA(3,:) ) < 1e-12 ) )

assertTrue(  all( abs( b.EoD(1,:) - b.EoD(2,:) ) < 1e-12 ) )
assertFalse( all( abs( b.EoD(1,:) - b.EoD(3,:) ) < 1e-12 ) )

assertTrue(  all( abs( b.EoA(1,:) - b.EoA(2,:) ) < 1e-12 ) )
assertFalse( all( abs( b.EoA(1,:) - b.EoA(3,:) ) < 1e-12 ) )

assertTrue(  all(all(abs( b.xprmat(:,:,1) - b.xprmat(:,:,2) ) < 1e-12)) );
assertFalse(  all(all(abs( b.xprmat(:,:,1) - b.xprmat(:,:,3) ) < 1e-12)) );

assertTrue(  all( abs( b.pin(1,:) - b.pin(2,:) ) < 1e-12 ) )
assertFalse( all( abs( b.pin(1,:) - b.pin(3,:) ) < 1e-12 ) )

c = b.get_channels;

assertTrue( all( abs( c(1,1).coeff(:) - c(1,2).coeff(:) ) < 1e-12 ) )