function testQF_calc_delay_spread

% Two paths with equal power: DS is half the delay difference, mean delay is in the middle
[ ds, mean_delay ] = qf.calc_delay_spread( [ 0, 100e-9 ], [ 1, 1 ] );
assertElementsAlmostEqual( ds, 50e-9, 'absolute', 1e-14 );
assertElementsAlmostEqual( mean_delay, 50e-9, 'absolute', 1e-14 );

% Multiple CIRs (rows), reference calculation
taus = [ 0, 10, 30, 70 ; 5, 6, 7, 8 ; 0, 0, 100, 100 ] * 1e-9;
pow  = [ 1, 0.5, 0.25, 0.125 ; 1, 1, 1, 1 ; 3, 1, 1, 3 ];
[ ds, mean_delay ] = qf.calc_delay_spread( taus, pow );
assertEqual( size(ds), [3,1] );
assertEqual( size(mean_delay), [3,1] );
pn = pow ./ repmat( sum(pow,2), 1, 4 );
mean_ref = sum( pn.*taus, 2 );
ds_ref = sqrt( sum( pn.*taus.^2, 2 ) - mean_ref.^2 );
assertElementsAlmostEqual( ds, ds_ref, 'absolute', 1e-14 );
assertElementsAlmostEqual( mean_delay, mean_ref, 'absolute', 1e-14 );

% Same powers for all CIRs
ds = qf.calc_delay_spread( taus, pow(1,:) );
pn = pow(1,:) ./ sum(pow(1,:));
pn = pn( [1,1,1],: );
ds_ref = sqrt( sum( pn.*taus.^2, 2 ) - sum( pn.*taus, 2 ).^2 );
assertElementsAlmostEqual( ds, ds_ref, 'absolute', 1e-14 );

% Threshold: the weak path at 1000 ns is removed with a 20 dB threshold
taus = [ 0, 100e-9, 1000e-9 ];
pow  = [ 1, 1, 0.001 ];
ds = qf.calc_delay_spread( taus, pow );
assertTrue( ds > 54e-9 );
ds = qf.calc_delay_spread( taus, pow, 20 );
assertElementsAlmostEqual( ds, 50e-9, 'absolute', 1e-14 );

% Granularity: paths at 100 and 110 ns are grouped into the same 50 ns bin
taus = [ 0, 100e-9, 110e-9 ];
pow  = [ 2, 1, 1 ];
ds = qf.calc_delay_spread( taus, pow, [], 50e-9 );
assertElementsAlmostEqual( ds, 50e-9, 'absolute', 1e-14 );

end
