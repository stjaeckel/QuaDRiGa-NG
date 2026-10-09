function testBuilder_gen_lsf_ssf

b = qd_builder('3GPP_38.901_UMi_LOS');
b.simpar.use_3GPP_baseline = 0;
b.simpar.show_progress_bars = 0;
b.simpar.center_frequency(2) = 10e9;
b.scenpar.SC_lambda = 0;

b.rx_positions = (rand(3,5)-0.5)*100;
b.rx_positions(3,:) = b.rx_positions(3,:)+51;
b.tx_position = [0;0;25];

b.gen_parameters(1);    % Gen LSF only

assertTrue( ~isempty( b.sos ) );
assertTrue( ~isempty( b.ds ) );
assertTrue( isempty( b.NumClusters ) );
assertTrue( isempty( b.taus ) );

b.kf(:) = 10;       % Fix KF

b.gen_parameters(2); % Gen SSF

assertTrue( ~isempty( b.taus ) );
assertTrue( isempty( b.fbs_pos ) );
assertTrue( isempty( b.lbs_pos ) );

% Check if KF comes from the SSF parameters
t = b.pow(:,1,:) ./ sum(b.pow(:,2:end,:),2);
assertTrue( all( abs(t(:) - 10) < 1e-13 ));

assertEqual( b.NumClusters, b.scenpar.NumClusters+4 );  % Split 2 strongest clusters into sub-clusters

b.gen_parameters(3); % Sub-clusters and scatterers

assertTrue( ~isempty( b.fbs_pos ) );
assertTrue( ~isempty( b.lbs_pos ) );

assertTrue( b.NumClusters > b.scenpar.NumClusters );

b.gen_parameters(0);

assertTrue( isempty( b.sos ) );
assertTrue( isempty( b.ds ) );
assertTrue( isempty( b.NumClusters ) );
assertTrue( isempty( b.taus ) );
assertTrue( isempty( b.fbs_pos ) );
assertTrue( isempty( b.lbs_pos ) );

b.gen_parameters(4);

assertTrue( b.NumClusters > b.scenpar.NumClusters );
assertTrue( ~isempty( b.sos ) );
assertTrue( ~isempty( b.ds ) );
assertTrue( ~isempty( b.NumClusters ) );
assertTrue( ~isempty( b.taus ) );
assertTrue( ~isempty( b.fbs_pos ) );
assertTrue( ~isempty( b.lbs_pos ) );

t = b.pow(:,1,:) ./ sum(b.pow(:,2:end,:),2);
assertTrue( all( abs(t(:) - 10) > 1e-13 ));

