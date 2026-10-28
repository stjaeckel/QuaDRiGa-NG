function testBuilder_lsf_model
% Large-scale parameters of the '38.901' and the 'legacy' model

set_rand_state( 3 );
rx = [ 30 100 100 400 1500 3000 ; 20 50 50 -300 800 -200 ; 1.5 1.5 22.5 1.5 4.5 1.5 ];

%% Correlation in the horizontal plane: equal LSPs for positions above each other
b = new_builder( '3GPP_38.901_UMa_NLOS', rx, '38.901' );
assertTrue( b.lsf_generator.horizontal );
xi = 10*log10( b.sf );
assertElementsAlmostEqual( xi(2), xi(3), 'absolute', 1e-12 );   % SF has no height terms
assertElementsAlmostEqual( log10( b.asA(2) ), log10( b.asA(3) ), 'absolute', 1e-12 );

b = new_builder( '3GPP_38.901_UMa_NLOS', rx, 'legacy' );
assertFalse( b.lsf_generator.horizontal );
xi = 10*log10( b.sf );
assertTrue( abs( xi(2) - xi(3) ) > 1e-6 );

%% Shadow-fading std of the path-loss model, the SF value is not mixed with the other LSPs
b = new_builder( '3GPP_38.901_RMa_LOS', rx, '38.901' );
[ ~, sf_sigma ] = b.get_pl;
assertEqual( unique( sf_sigma ), [4, 6] );                        % Before and after the breakpoint
g = b.lsf_generator;
rnd = quadriga_lib.lsf_draw( g.lambda, g.phase, b.tx_position.', rx.', [], false, g.preset, g.horizontal );
assertElementsAlmostEqual( 10*log10( b.sf ), rnd(:,3).' .* sf_sigma, 'absolute', 1e-10 );

sf = b.sf;
M = b.lsp_xcorr;
M(1,3) = -0.6; M(3,1) = -0.6;
b.lsp_xcorr = M;
gen_lsf_parameters( b );
assertElementsAlmostEqual( b.sf, sf, 'relative', 1e-12 );         % Cholesky: SF first

b = new_builder( '3GPP_38.901_RMa_LOS', rx, 'legacy' );
xi = 10*log10( b.sf );
b.lsp_xcorr = M;
gen_lsf_parameters( b );
assertTrue( max( abs( 10*log10( b.sf ) - xi ) ) > 1e-3 );         % Symmetric square root

%% Limits of the angular spreads
b = qd_builder( '3GPP_38.901_UMa_NLOS' );
b.scenpar.AS_A_mu = 2.1;
b.scenpar.ES_A_mu = 1.8;
b = new_builder( b, rx, '38.901' );
assertTrue( all( b.asA <= 104 ) && any( b.asA == 104 ) );
assertTrue( all( b.esA <= 52 ) && any( b.esA == 52 ) );
b.simpar.lsf_model = 'legacy';
gen_lsf_parameters( b );
assertTrue( any( b.asA > 104 ) );

%% SF profile along a track: equal to the builder value at the segment start
b = new_builder( '3GPP_38.901_RMa_LOS', rx(:,1), '38.901' );
t = qd_track( 'linear', 4000, 0 );
t.initial_position = rx(:,1);
t.interpolate_positions( 0.05 );
sf = b.get_sf_profile( t );
assertElementsAlmostEqual( sf(1), b.sf, 'relative', 1e-12 );
[ ~, sf_sigma ] = b.get_pl( t );
assertEqual( unique( sf_sigma ), [4, 6] );                        % The track crosses the breakpoint
g = b.lsf_generator;
pos = t.positions + t.initial_position * ones( 1, t.no_snapshots );
rnd = quadriga_lib.lsf_draw( g.lambda, g.phase, b.tx_position.', pos.', [], false, g.preset, g.horizontal );
assertElementsAlmostEqual( 10*log10( sf ), 4 * rnd(:,3).', 'absolute', 1e-10 );  % Std of the segment start

%% Spatial consistency disabled for the ground permittivity and the absolute TOA
b = qd_builder( '3GPP_38.901_InF_NLOS_DH' );
b.simpar.show_progress_bars = 0;
b.simpar.autocorrelation_function = 'Disable';
b.tx_position = [ 0; 0; 5 ];
b.rx_positions = rx(:,1:3);
init_sos( b );
assertEqual( b.lsf_generator.preset, 'Comb256' );
assertEqual( b.lsf_generator.lambda.SC_lambda, 0 );
assertEqual( b.lsf_generator.lambda.absTOA_lambda, 0 );

%% No Rx positions
b = qd_builder( '3GPP_38.901_UMa_NLOS' );
b.simpar.show_progress_bars = 0;
b.tx_position = [ 0; 0; 25 ];
init_sos( b );
gen_lsf_parameters( b );
assertTrue( isempty( b.ds ) );

end

function b = new_builder( scen, rx, model )
if ischar( scen )
    b = qd_builder( scen );
else
    b = scen;
end
b.simpar.show_progress_bars = 0;
b.simpar.lsf_model = model;
b.simpar.center_frequency = 3.5e9;
b.tx_position = [ 0; 0; 25 ];
b.rx_positions = rx;
b.dual_mobility = 0;
init_sos( b );
gen_lsf_parameters( b );
end
