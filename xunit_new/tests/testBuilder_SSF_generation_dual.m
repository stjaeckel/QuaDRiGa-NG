function testBuilder_SSF_generation_dual
%%

set_rand_state(2);

b = qd_builder('3GPP_38.901_UMi_NLOS'); % Must have spatial consistency enabled

% Full-Dual mobility is only achieved if parameters are symmetric
scp = b.scenpar;
scp.AS_D_lambda = scp.AS_A_lambda;
scp.AS_D_mu = scp.AS_A_mu;
scp.AS_D_gamma = scp.AS_A_gamma;
scp.AS_D_sigma = scp.AS_A_sigma;
scp.AS_D_delta = scp.AS_A_delta;
scp.ES_D_lambda = scp.ES_A_lambda;
scp.ES_D_mu = scp.ES_A_mu;
scp.ES_D_gamma = scp.ES_A_gamma;
scp.ES_D_sigma = scp.ES_A_sigma;
scp.ES_D_delta = scp.ES_A_delta;
scp.ES_D_omega = scp.ES_A_omega;
scp.ES_D_mu_A = 0;
scp.asA_ds = 0;
scp.asA_sf = 0;
scp.esD_ds = 0;
scp.esA_asA = 0;
scp.esD_asD = 0;
scp.esA_asD = 0;
scp.PerClusterES_A = 0;
scp.PerClusterES_D = 0;
scp.PerClusterAS_A = 0;
scp.PerClusterAS_D = 0;

b.scenpar = scp;

b.lsp_xcorr = eye(8);
b.simpar.show_progress_bars = 0;
b.simpar.use_absolute_delays = 1;
b.rx_positions = [0,10,10;0,500,0.1;0,0,1]';
b.tx_position = b.rx_positions(:,[2 1 3]);

gen_parameters(b);

assertEqual( b.dual_mobility, true );

assertTrue( all( abs( b.taus(1,:) - b.taus(2,:) ) < 1e-12 ) );      % Same delays for swapped TX-Rx
assertTrue( all( abs( b.pow(1,:) - b.pow(2,:) ) < 1e-12 ) );        % Same power for swapped Tx-Rx
assertTrue( all( abs( b.AoA(1,:) - b.AoD(2,:) ) < 1e-7 ) );        % Swapped AoA and AoD for swapped Tx-Rx
assertTrue( all( abs( b.AoD(1,:) - b.AoA(2,:) ) < 1e-7 ) );        % Swapped AoD and AoA for swapped Tx-Rx
assertTrue( all( abs( b.EoA(1,:) - b.EoD(2,:) ) < 1e-7 ) );        % Swapped EoA and EoD for swapped Tx-Rx
assertTrue( all( abs( b.EoD(1,:) - b.EoA(2,:) ) < 1e-7 ) );        % Swapped AoD and EoA for swapped Tx-Rx

