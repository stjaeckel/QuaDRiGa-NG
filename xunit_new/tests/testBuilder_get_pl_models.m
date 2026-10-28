function testBuilder_get_pl_models
% Path-loss models: WINNER LOS beyond the breakpoint, O2I loss, position input, multi-frequency

%% WINNER LOS beyond the breakpoint: C2 and per-link heights
b = qd_builder( 'WINNER_UMa_C2_LOS' );
b.simpar.show_progress_bars = 0;
b.simpar.center_frequency = 2.6e9;
b.tx_position = [ 0; 0; 25 ];
b.rx_positions = [ 3000, 0, 1.5; 0, 5000, 1.0; -4000, 0, 2.5 ].';
p = b.plpar;
hb = 25;
hm = max( b.rx_positions(3,:), 1.5 );
d2 = sqrt( sum( b.rx_positions(1:2,:).^2, 1 ) );
pl_ref = p.A2*log10(d2) + p.B2 + p.C2*log10(2.6) + p.D2*log10(hb) + p.E2*log10(hm) + p.F2*hm;
assertElementsAlmostEqual( b.get_pl, pl_ref, 'absolute', 1e-10 );

%% The O2I loss is added once, also for the 'nlos' model
b = qd_builder( '3GPP_38.901_UMa_NLOS' );
b.simpar.show_progress_bars = 0;
b.simpar.center_frequency = [ 3.5e9, 28e9 ];
b.tx_position = [ 0; 0; 25 ];
b.rx_positions = [ 15; 5; 1.5 ];
t = qd_track( [] );
t.positions = [ 0:10:200; zeros(2,21) ];
t.initial_position = [ 15; 5; 1.5 ];
pl_out = b.get_pl( t );
t.par = struct( 'o2i_loss', reshape( [ 15, 20 ], 1, 1, 2 ), 'o2i_d2din', 12 );
[ pl_in, ~, pl_init ] = b.get_pl( t );
assertElementsAlmostEqual( pl_in - pl_out, [ 15; 20 ] * ones(1,21), 'absolute', 1e-10 );
assertElementsAlmostEqual( pl_init, pl_in(:,1), 'absolute', 1e-10 );

% Legacy model: the outdoor distance is shorter than the geometric distance
b.simpar.lsf_model = 'legacy';
pl_legacy = b.get_pl( t );
assertTrue( all( pl_legacy(:) < pl_in(:) ) );

%% Position input: 'nlos' model at the given positions
b.simpar.lsf_model = '38.901';
b.rx_positions = [ 2000 2500 3000 ; 0 0 0 ; 1.5 1.5 1.5 ];
pos = [ 20 30 40 50 ; 0 0 0 0 ; 1.5 1.5 1.5 1.5 ];
pl_pos = b.get_pl( pos );
b.rx_positions = pos;
assertEqual( pl_pos, b.get_pl );

%% Satellite model at two frequencies
b = qd_builder( 'QuaDRiGa_NTN_Rural_LOS' );
b.simpar.show_progress_bars = 0;
b.tx_position = [ 300e3; 0; 600e3 ];
b.rx_positions = [ 0 100 ; 0 0 ; 1.5 1.5 ];
b.simpar.center_frequency = [ 2e9, 20e9 ];
pl = b.get_pl;
b.simpar.center_frequency = 2e9;
pl1 = b.get_pl;
b.simpar.center_frequency = 20e9;
pl2 = b.get_pl;
assertElementsAlmostEqual( pl, [ pl1; pl2 ], 'absolute', 1e-10 );

%% No Rx positions
b.rx_positions = zeros( 3,0 );
assertEqual( size( b.get_pl ), [ 1, 0 ] );

end
