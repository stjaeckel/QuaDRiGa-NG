function testBuilder_mmMAGIC_subpaths
%%

b = qd_builder('mmMAGIC_UMi_NLOS');
b.simpar.show_progress_bars = 0;
b.simpar.center_frequency(2) = 20e9;
b.tx_array = qd_arrayant;
b.rx_array = b.tx_array;

assertEqual( b.scenpar.SubpathMethod, 'mmMAGIC' );

b.tx_position = [0;0;10];
b.rx_positions = [20 20 ; 0 0 ; 1.5 1.5 ];

b.gen_parameters;

assertEqual( numel( b.NumSubPaths ), b.NumClusters );
assertTrue( all( abs( b.taus(1,:) - b.taus(2,:) ) < 1e-12 ) ); % Spatial consistency

for n = 1 : b.scenpar.NumClusters-1
    ii = b.scenpar.NumSubPaths*(n-1)+2 : b.scenpar.NumSubPaths*n+1;
    a = b.AoA(1,ii);
    a = a - angle(mean(exp(1j*a)));
    a = angle(exp(1j*a));
    a = a-mean(a);
    assertTrue( abs( a(1)+a(end) ) < 1e-11 );
end



