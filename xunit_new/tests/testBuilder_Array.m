function testBuilder_Array
%%
b = qd_builder('3GPP_38.901_UMi_LOS_GR');
b.simpar.show_progress_bars = 0;
b(1,2) = qd_builder('LOSonly');
b(2,1) = qd_builder('TwoRayGR');
b(2,2) = qd_builder('Ul');

for n = 1:2
    for m = 1:2
        b(n,m).tx_position = [0;0;25];
        b(n,m).rx_positions = [100,0,15]';
    end
end

gen_parameters(b);
c = get_channels( b );

assertEqual( size(c),[1,4])
assertEqual( c(1,1).no_path , b(1,1).NumClusters );
assertEqual( c(1,2).no_path , b(2,1).NumClusters );
assertEqual( c(1,3).no_path , b(1,2).NumClusters );
assertEqual( c(1,4).no_path , b(2,2).NumClusters );