function testBuilder_gen_subpaths
%%

% legacy, PerClusterDS not Freq-dependent
b = qd_builder('3GPP_38.901_UMi_LOS',0); 
b.simpar.center_frequency = [6e9 20e9];
b.simpar.show_progress_bars = 0;
b.scenpar.SC_lambda = 0;
b.rx_positions = [20,0,0 ; 200.1,0,0]'; 
b.tx_position = [0;0;25];
b.gen_parameters; 

% legacy, PerClusterDS Freq-dependent
b = qd_builder('3GPP_38.901_UMa_NLOS',0); 
b.simpar.center_frequency = [6e9 20e9];
b.simpar.show_progress_bars = 0;
b.scenpar.SC_lambda = 0;
b.rx_positions = [20,0,0 ; 200.1,0,0]'; 
b.tx_position = [0;0;25];
b.gen_parameters; 

% mmMAGIC, PerClusterDS Freq-dependent
b = qd_builder('3GPP_38.901_UMa_NLOS',0); 
b.simpar.center_frequency = [6e9 20e9];
b.simpar.show_progress_bars = 0;
b.scenpar.SubpathMethod = 'mmMAGIC';
b.scenpar.NumClusters = 6;
b.scenpar.NumSubPaths = 7;
b.scenpar.SC_lambda = 0;
b.rx_positions = [20,0,0 ; 200.1,0,0]'; 
b.tx_position = [0;0;25];
b.gen_parameters; 

b.scenpar.NumClusters = 6;

% mmMAGIC, PerClusterDS not Freq-dependent
b = qd_builder('mmMAGIC_UMi_NLOS',0); 
b.simpar.center_frequency = [6e9 20e9];
b.simpar.show_progress_bars = 0;
b.scenpar.SC_lambda = 0;
b.rx_positions = [20,0,0 ; 200.1,0,0]'; 
b.tx_position = [0;0;25];
b.gen_parameters; 

