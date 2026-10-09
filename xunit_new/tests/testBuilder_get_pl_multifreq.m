function testBuilder_get_pl_multifreq
%%

conf_files = { 'Freespace', ...             % logdist
    'BERLIN_UMa_LOS' , ...                  % logdist_simple
    'MIMOSA_10-45_LOS' , ...                % constant
    'WINNER_UMa_C2_LOS' ,...                % winner_los
    'WINNER_SMa_C1_NLOS',...                % winner_nlos
    'WINNER_Indoor_A1_NLOS' ,...            % winner_pathloss
	'3GPP_38.901_UMi_LOS' ,...              % dual-slope
    '3GPP_38.901_UMi_LOS_GR',...            % tripple-slope
    '3GPP_3D_UMa_NLOS',...                  % 3gpp_3d_uma_nlos
    '3GPP_38.901_UMi_NLOS'};                % 3gpp_3d_umi_nlos

evaltrack = qd_track( 'linear', 100, 0 );
evaltrack.initial_position = [100;0;1.5];

rx_pos = [100 200 ; 0 0 ; 1.5 1.5 ];

for n = 1 : numel( conf_files )
    b = qd_builder( conf_files{n});
    
    b.tx_position = [0;0;4];
    b.rx_positions = rx_pos;
    b.simpar.center_frequency = [ 0.5e9 1.9e9 2.6e9 28e9];
    
    b.check_dual_mobility;
    
    [ l1,SF_sigma1] = b.get_pl;
    [ l2,SF_sigma2] = b.get_pl( evaltrack,[], b.tx_position(:,1)  );
    
    assertEqual( size(l1),[4,2]);
    assertEqual( size(SF_sigma1),[4,2] );
    assertEqual( size(SF_sigma2),[4,2] );
    assertEqual( l1,l2 );
end
