function testBuilder_ZeroPower
%%
b = qd_builder('3GPP_38.901_UMi_LOS_GR');
b.simpar.show_progress_bars = 0;
b.tx_position = [0;0;25];
b.rx_positions = [100,0,15;102,0,15;104,0,15;105,0,15]';
gen_parameters(b);

ii = randperm(b.NumClusters-2)+2;
b.pow(1,ii(1:12)) = 0;
b.pow(1,end-1:end) = 0;
b.pow(2,2) = 0;
b.pow(3,1) = 0;

c = get_channels( b );

assertEqual( c(1,1).no_path, sum(b.pow(1,:)~=0) );
assertEqual( c(1,2).no_path, b.NumClusters );
assertFalse( isnan( c(1,2).coeff(:,:,2,:) ) );
assertEqual( c(1,2).coeff(:,:,2,:), 0 );
assertEqual( c(1,3).no_path, b.NumClusters );
assertFalse( isnan( c(1,3).coeff(:,:,1,:) ) );
assertEqual( c(1,3).coeff(:,:,1,:), 0 );
assertEqual( c(1,4).no_path, b.NumClusters );

% Baseline model
b = qd_builder('3GPP_38.901_UMi_LOS');
b.simpar.use_3GPP_baseline = 1;
b.simpar.show_progress_bars = 0;
b.tx_position = [0;0;25];
b.rx_positions = [100,0,15;102,0,15;104,0,15]';
gen_parameters(b);

ii = randperm(b.NumClusters-1)+1;
b.pow(1,ii(1:12)) = 0;
b.pow(1,end-1:end) = 0;
b.pow(2,1) = 0;

c = get_channels( b );

assertEqual( c(1,1).no_path, sum(b.pow(1,:)~=0) );
assertEqual( c(1,2).no_path, b.NumClusters );
assertEqual( c(1,2).no_path, b.NumClusters );
assertFalse( isnan( c(1,2).coeff(:,:,1,:) ) );
assertEqual( c(1,2).coeff(:,:,1,:), 0 );
assertEqual( c(1,3).no_path, b.NumClusters );



