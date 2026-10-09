function testBuilder_absTOA_model

% Test setting the 3 model parameters
b = qd_builder('3GPP_38.901_InF_NLOS_DH');
b.simpar.show_progress_bars = 0;

assertTrue( isempty( b.absTOA_sos ));

b.scenpar.SC_lambda = 0;
b.scenpar.PerClusterDS = 5;
b.scenpar.absTOA_mu = -7+rand;
b.scenpar.absTOA_sigma = rand;
b.scenpar.absTOA_lambda = rand*10;

b.write_conf_file('bla.conf');

b1 = qd_builder('bla');
assertTrue( abs( b1.scenpar.absTOA_mu - b.scenpar.absTOA_mu ) < 1e-4 )
assertTrue( abs( b1.scenpar.absTOA_sigma - b.scenpar.absTOA_sigma ) < 1e-4 )
assertTrue( abs( b1.scenpar.absTOA_lambda - b.scenpar.absTOA_lambda ) < 1e-4 )

delete('bla.conf')

b.tx_position = [ [0;0;5]*[1 1 1],[10;0;1]];
b.rx_positions = [ 10,0,1 ; 10,0,1 ; 0,10,1 ; 0,0,5  ]';

b.check_dual_mobility;

b.init_sos;
assertEqual( b.absTOA_sos.dist_decorr, b.scenpar.absTOA_lambda );
assertTrue(  all( abs(b.absTOA_sos.sos_phase(:,1) - b.absTOA_sos.sos_phase(:,2)) < 1e-16 ) ); % Channel reciprocity

b.gen_parameters;

ii = b.taus(:,2:end) >= (b.absTOA_offset'*ones(1,b.NumClusters-1));
assertTrue( all(ii(:)))

assertEqual( b.NumClusters, b.scenpar.NumClusters + 4 ); % Split only 2 clusters

assertTrue(  abs( b.absTOA_offset(1) - b.absTOA_offset(2) ) < 1e-15 )
assertTrue(  abs( b.absTOA_offset(1) - b.absTOA_offset(4) ) < 1e-15 )
assertTrue(  abs( b.absTOA_offset(1) - b.absTOA_offset(3) ) > 1e-15 )

b = qd_builder('3GPP_38.901_InF_NLOS_DH');
b.simpar.show_progress_bars = 0;
b.simpar.use_3GPP_baseline = 1;
b.tx_position = [0;0;5];
b.rx_positions = [ 10,0,1 ; 10,0,1 ; 0,10,1 ]';

b.scenpar.PerClusterDS = 5;

b.gen_parameters;

assertTrue( ~isempty( b.absTOA_offset ));
assertEqual( b.NumClusters, b.scenpar.NumClusters + 4 ); % Split only 2 clusters



