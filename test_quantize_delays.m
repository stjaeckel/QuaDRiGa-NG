%TEST_QUANTIZE_DELAYS  Compare original quantize_delays with quadriga_lib wrapper
%
%   Runs both implementations on identical synthetic channels and checks
%   that the outputs match within floating-point tolerance.

clear; clc;
fprintf('========== quantize_delays : original vs. quadriga_lib ==========\n\n');

rng(42);                              % Reproducible random data
cnt = [0 0];                          % [pass, fail]

% ---- Test parameters -------------------------------------------------------
tap_spacing  = 5e-9;                  % 5 ns  (200 MHz sample rate)
tol_coeff    = 2e-5;                  % Tolerance for coefficient comparison
tol_delay    = 1e-13;                 % Tolerance for delay comparison (quantised)

% ============================================================================
%% Test 1 - Basic operation, default parameters (fix_taps = 0)
% ============================================================================
fprintf('--- Test 1: fix_taps = 0, default parameters ---\n');

no_rx   = 2;
no_tx   = 2;
no_path = 8;
no_snap = 10;

coeff = complex( randn(no_rx,no_tx,no_path,no_snap,'single'), ...
                 randn(no_rx,no_tx,no_path,no_snap,'single') );
delay = single( sort( abs(randn(no_rx,no_tx,no_path,no_snap,'single'))*50e-9, 3 ) );

h = qd_channel( coeff, delay );
h.name = 'test1';
h.individual_delays = true;
h.center_frequency  = 3.5e9;
h.tx_position = [0;0;25];
h.rx_position = [100;0;1.5];

h1 = quantize_delays( copy(h), tap_spacing, [], [], [], 0, 0 );
h2 = quantize_delays_new( copy(h), tap_spacing, [], [], [], 0, 0 );

ok = isequal(size(h1.coeff),size(h2.coeff));
if ok; fprintf('  PASS  coeff size\n'); cnt(1)=cnt(1)+1;
else;  fprintf('  FAIL  coeff size  (%s vs %s)\n',mat2str(size(h1.coeff)),mat2str(size(h2.coeff))); cnt(2)=cnt(2)+1; end

ok = isequal(size(h1.delay),size(h2.delay));
if ok; fprintf('  PASS  delay size\n'); cnt(1)=cnt(1)+1;
else;  fprintf('  FAIL  delay size  (%s vs %s)\n',mat2str(size(h1.delay)),mat2str(size(h2.delay))); cnt(2)=cnt(2)+1; end

err_c = max(abs( single(h1.coeff(:)) - single(h2.coeff(:)) ));
ok = err_c < tol_coeff;
if ok; fprintf('  PASS  coeff match\n'); cnt(1)=cnt(1)+1;
else;  fprintf('  FAIL  coeff match  (max err = %g)\n',err_c); cnt(2)=cnt(2)+1; end

err_d = max(abs( double(h1.delay(:)) - double(h2.delay(:)) ));
ok = err_d < tol_delay;
if ok; fprintf('  PASS  delay match\n'); cnt(1)=cnt(1)+1;
else;  fprintf('  FAIL  delay match  (max err = %g)\n',err_d); cnt(2)=cnt(2)+1; end

ok = strcmp(h1.name, h2.name);
if ok; fprintf('  PASS  name\n'); cnt(1)=cnt(1)+1;
else;  fprintf('  FAIL  name\n'); cnt(2)=cnt(2)+1; end

ok = h1.center_frequency == h2.center_frequency;
if ok; fprintf('  PASS  fc\n'); cnt(1)=cnt(1)+1;
else;  fprintf('  FAIL  fc\n'); cnt(2)=cnt(2)+1; end

ok = isequal(h1.tx_position, h2.tx_position);
if ok; fprintf('  PASS  tx_pos\n'); cnt(1)=cnt(1)+1;
else;  fprintf('  FAIL  tx_pos\n'); cnt(2)=cnt(2)+1; end

ok = isequal(h1.rx_position, h2.rx_position);
if ok; fprintf('  PASS  rx_pos\n'); cnt(1)=cnt(1)+1;
else;  fprintf('  FAIL  rx_pos\n'); cnt(2)=cnt(2)+1; end
fprintf('\n');

% ============================================================================
%% Test 2 - fix_taps = 1 (single grid for all)
% ============================================================================
fprintf('--- Test 2: fix_taps = 1 ---\n');

h1 = quantize_delays( copy(h), tap_spacing, [], [], [], 1, 0 );
h2 = quantize_delays_new( copy(h), tap_spacing, [], [], [], 1, 0 );

ok = isequal(size(h1.coeff),size(h2.coeff));
if ok; fprintf('  PASS  coeff size\n'); cnt(1)=cnt(1)+1;
else;  fprintf('  FAIL  coeff size  (%s vs %s)\n',mat2str(size(h1.coeff)),mat2str(size(h2.coeff))); cnt(2)=cnt(2)+1; end

err_c = max(abs( single(h1.coeff(:)) - single(h2.coeff(:)) ));
ok = err_c < tol_coeff;
if ok; fprintf('  PASS  coeff match\n'); cnt(1)=cnt(1)+1;
else;  fprintf('  FAIL  coeff match  (max err = %g)\n',err_c); cnt(2)=cnt(2)+1; end

err_d = max(abs( double(h1.delay(:)) - double(h2.delay(:)) ));
ok = err_d < tol_delay;
if ok; fprintf('  PASS  delay match\n'); cnt(1)=cnt(1)+1;
else;  fprintf('  FAIL  delay match  (max err = %g)\n',err_d); cnt(2)=cnt(2)+1; end
fprintf('\n');

% ============================================================================
%% Test 3 - fix_taps = 2 (per snapshot)
% ============================================================================
fprintf('--- Test 3: fix_taps = 2 ---\n');

h1 = quantize_delays( copy(h), tap_spacing, [], [], [], 2, 0 );
h2 = quantize_delays_new( copy(h), tap_spacing, [], [], [], 2, 0 );

ok = isequal(size(h1.coeff),size(h2.coeff));
if ok; fprintf('  PASS  coeff size\n'); cnt(1)=cnt(1)+1;
else;  fprintf('  FAIL  coeff size  (%s vs %s)\n',mat2str(size(h1.coeff)),mat2str(size(h2.coeff))); cnt(2)=cnt(2)+1; end

err_c = max(abs( single(h1.coeff(:)) - single(h2.coeff(:)) ));
ok = err_c < tol_coeff;
if ok; fprintf('  PASS  coeff match\n'); cnt(1)=cnt(1)+1;
else;  fprintf('  FAIL  coeff match  (max err = %g)\n',err_c); cnt(2)=cnt(2)+1; end

err_d = max(abs( double(h1.delay(:)) - double(h2.delay(:)) ));
ok = err_d < tol_delay;
if ok; fprintf('  PASS  delay match\n'); cnt(1)=cnt(1)+1;
else;  fprintf('  FAIL  delay match  (max err = %g)\n',err_d); cnt(2)=cnt(2)+1; end
fprintf('\n');

% ============================================================================
%% Test 4 - fix_taps = 3 (per tx-rx pair)
% ============================================================================
fprintf('--- Test 4: fix_taps = 3 ---\n');

h1 = quantize_delays( copy(h), tap_spacing, [], [], [], 3, 0 );
h2 = quantize_delays_new( copy(h), tap_spacing, [], [], [], 3, 0 );

ok = isequal(size(h1.coeff),size(h2.coeff));
if ok; fprintf('  PASS  coeff size\n'); cnt(1)=cnt(1)+1;
else;  fprintf('  FAIL  coeff size  (%s vs %s)\n',mat2str(size(h1.coeff)),mat2str(size(h2.coeff))); cnt(2)=cnt(2)+1; end

err_c = max(abs( single(h1.coeff(:)) - single(h2.coeff(:)) ));
ok = err_c < tol_coeff;
if ok; fprintf('  PASS  coeff match\n'); cnt(1)=cnt(1)+1;
else;  fprintf('  FAIL  coeff match  (max err = %g)\n',err_c); cnt(2)=cnt(2)+1; end

err_d = max(abs( double(h1.delay(:)) - double(h2.delay(:)) ));
ok = err_d < tol_delay;
if ok; fprintf('  PASS  delay match\n'); cnt(1)=cnt(1)+1;
else;  fprintf('  FAIL  delay match  (max err = %g)\n',err_d); cnt(2)=cnt(2)+1; end
fprintf('\n');

% ============================================================================
%% Test 5 - Limited taps (max_no_taps = 6)
% ============================================================================
fprintf('--- Test 5: max_no_taps = 6, fix_taps = 0 ---\n');

h1 = quantize_delays( copy(h), tap_spacing, 6, [], [], 0, 0 );
h2 = quantize_delays_new( copy(h), tap_spacing, 6, [], [], 0, 0 );

ok = isequal(size(h1.coeff),size(h2.coeff));
if ok; fprintf('  PASS  coeff size\n'); cnt(1)=cnt(1)+1;
else;  fprintf('  FAIL  coeff size  (%s vs %s)\n',mat2str(size(h1.coeff)),mat2str(size(h2.coeff))); cnt(2)=cnt(2)+1; end

err_c = max(abs( single(h1.coeff(:)) - single(h2.coeff(:)) ));
ok = err_c < tol_coeff;
if ok; fprintf('  PASS  coeff match\n'); cnt(1)=cnt(1)+1;
else;  fprintf('  FAIL  coeff match  (max err = %g)\n',err_c); cnt(2)=cnt(2)+1; end

err_d = max(abs( double(h1.delay(:)) - double(h2.delay(:)) ));
ok = err_d < tol_delay;
if ok; fprintf('  PASS  delay match\n'); cnt(1)=cnt(1)+1;
else;  fprintf('  FAIL  delay match  (max err = %g)\n',err_d); cnt(2)=cnt(2)+1; end
fprintf('\n');

% ============================================================================
%% Test 6 - Shared delays (individual_delays = false)
% ============================================================================
fprintf('--- Test 6: shared delays (individual_delays = false) ---\n');

delay_shared = single( sort( abs(randn(1,1,no_path,no_snap,'single'))*50e-9, 3 ) );
delay_full   = repmat( delay_shared, [no_rx, no_tx, 1, 1] );

h_shared = qd_channel( coeff, delay_full );
h_shared.name = 'test6_shared';
h_shared.individual_delays = false;
h_shared.center_frequency  = 3.5e9;
h_shared.tx_position = [0;0;25];
h_shared.rx_position = [100;0;1.5];

h1 = quantize_delays( copy(h_shared), tap_spacing, [], [], [], 0, 0 );
h2 = quantize_delays_new( copy(h_shared), tap_spacing, [], [], [], 0, 0 );

ok = isequal(size(h1.coeff),size(h2.coeff));
if ok; fprintf('  PASS  coeff size\n'); cnt(1)=cnt(1)+1;
else;  fprintf('  FAIL  coeff size  (%s vs %s)\n',mat2str(size(h1.coeff)),mat2str(size(h2.coeff))); cnt(2)=cnt(2)+1; end

err_c = max(abs( single(h1.coeff(:)) - single(h2.coeff(:)) ));
ok = err_c < tol_coeff;
if ok; fprintf('  PASS  coeff match\n'); cnt(1)=cnt(1)+1;
else;  fprintf('  FAIL  coeff match  (max err = %g)\n',err_c); cnt(2)=cnt(2)+1; end

err_d = max(abs( double(h1.delay(:)) - double(h2.delay(:)) ));
ok = err_d < tol_delay;
if ok; fprintf('  PASS  delay match\n'); cnt(1)=cnt(1)+1;
else;  fprintf('  FAIL  delay match  (max err = %g)\n',err_d); cnt(2)=cnt(2)+1; end
fprintf('\n');

% ============================================================================
%% Test 7 - Antenna sub-selection (i_rxant, i_txant)
% ============================================================================
fprintf('--- Test 7: antenna sub-selection ---\n');

no_rx_big = 4;  no_tx_big = 4;
coeff_big = complex( randn(no_rx_big,no_tx_big,no_path,no_snap,'single'), ...
                     randn(no_rx_big,no_tx_big,no_path,no_snap,'single') );
delay_big = single( sort( abs(randn(no_rx_big,no_tx_big,no_path,no_snap,'single'))*50e-9, 3 ) );

h_big = qd_channel( coeff_big, delay_big );
h_big.name = 'test7_subsel';
h_big.individual_delays = true;
h_big.center_frequency  = 3.5e9;
h_big.tx_position = [0;0;25];
h_big.rx_position = [100;0;1.5];

i_rx = [1 3];  i_tx = [2 4];

h1 = quantize_delays( copy(h_big), tap_spacing, [], i_rx, i_tx, 0, 0 );
h2 = quantize_delays_new( copy(h_big), tap_spacing, [], i_rx, i_tx, 0, 0 );

ok = isequal(size(h1.coeff),size(h2.coeff));
if ok; fprintf('  PASS  coeff size\n'); cnt(1)=cnt(1)+1;
else;  fprintf('  FAIL  coeff size  (%s vs %s)\n',mat2str(size(h1.coeff)),mat2str(size(h2.coeff))); cnt(2)=cnt(2)+1; end

err_c = max(abs( single(h1.coeff(:)) - single(h2.coeff(:)) ));
ok = err_c < tol_coeff;
if ok; fprintf('  PASS  coeff match\n'); cnt(1)=cnt(1)+1;
else;  fprintf('  FAIL  coeff match  (max err = %g)\n',err_c); cnt(2)=cnt(2)+1; end

err_d = max(abs( double(h1.delay(:)) - double(h2.delay(:)) ));
ok = err_d < tol_delay;
if ok; fprintf('  PASS  delay match\n'); cnt(1)=cnt(1)+1;
else;  fprintf('  FAIL  delay match  (max err = %g)\n',err_d); cnt(2)=cnt(2)+1; end
fprintf('\n');

% ============================================================================
%% Test 8 - Already-quantised input (delays on grid)
% ============================================================================
fprintf('--- Test 8: already-quantised input ---\n');

delay_quant = single( round( abs(randn(no_rx,no_tx,no_path,no_snap,'single'))*50e-9 ...
    / tap_spacing ) * tap_spacing );
delay_quant = sort( delay_quant, 3 );

h_q = qd_channel( coeff, delay_quant );
h_q.name = 'test8_prequant';
h_q.individual_delays = true;
h_q.center_frequency  = 3.5e9;
h_q.tx_position = [0;0;25];
h_q.rx_position = [100;0;1.5];

h1 = quantize_delays( copy(h_q), tap_spacing, [], [], [], 0, 0 );
h2 = quantize_delays_new( copy(h_q), tap_spacing, [], [], [], 0, 0 );

ok = isequal(size(h1.coeff),size(h2.coeff));
if ok; fprintf('  PASS  coeff size\n'); cnt(1)=cnt(1)+1;
else;  fprintf('  FAIL  coeff size  (%s vs %s)\n',mat2str(size(h1.coeff)),mat2str(size(h2.coeff))); cnt(2)=cnt(2)+1; end

err_c = max(abs( single(h1.coeff(:)) - single(h2.coeff(:)) ));
ok = err_c < tol_coeff;
if ok; fprintf('  PASS  coeff match\n'); cnt(1)=cnt(1)+1;
else;  fprintf('  FAIL  coeff match  (max err = %g)\n',err_c); cnt(2)=cnt(2)+1; end

err_d = max(abs( double(h1.delay(:)) - double(h2.delay(:)) ));
ok = err_d < tol_delay;
if ok; fprintf('  PASS  delay match\n'); cnt(1)=cnt(1)+1;
else;  fprintf('  FAIL  delay match  (max err = %g)\n',err_d); cnt(2)=cnt(2)+1; end
fprintf('\n');

% ============================================================================
%% Test 9 - Power conservation check
% ============================================================================
fprintf('--- Test 9: power conservation ---\n');

P_in  = sum( abs(coeff(:)).^2 );

h_out1 = quantize_delays( copy(h), tap_spacing, [], [], [], 0, 0 );
h_out2 = quantize_delays_new( copy(h), tap_spacing, [], [], [], 0, 0 );

P_out1 = sum( abs( single(h_out1.coeff(:)) ).^2 );
P_out2 = sum( abs( single(h_out2.coeff(:)) ).^2 );

rel_err1 = abs(P_out1 - P_in) / P_in;
ok = rel_err1 < 0.02;
if ok; fprintf('  PASS  power orig\n'); cnt(1)=cnt(1)+1;
else;  fprintf('  FAIL  power orig  (rel err = %g)\n',rel_err1); cnt(2)=cnt(2)+1; end

rel_err2 = abs(P_out2 - P_in) / P_in;
ok = rel_err2 < 0.02;
if ok; fprintf('  PASS  power new\n'); cnt(1)=cnt(1)+1;
else;  fprintf('  FAIL  power new  (rel err = %g)\n',rel_err2); cnt(2)=cnt(2)+1; end

ok = abs(P_out1 - P_out2)/P_in < tol_coeff;
if ok; fprintf('  PASS  power match\n'); cnt(1)=cnt(1)+1;
else;  fprintf('  FAIL  power match  (rel diff = %g)\n',abs(P_out1-P_out2)/P_in); cnt(2)=cnt(2)+1; end
fprintf('\n');

% ============================================================================
%% Test 10 - Different tap spacing (10 ns)
% ============================================================================
fprintf('--- Test 10: tap_spacing = 10 ns ---\n');

ts2 = 10e-9;
h1 = quantize_delays( copy(h), ts2, [], [], [], 0, 0 );
h2 = quantize_delays_new( copy(h), ts2, [], [], [], 0, 0 );

ok = isequal(size(h1.coeff),size(h2.coeff));
if ok; fprintf('  PASS  coeff size\n'); cnt(1)=cnt(1)+1;
else;  fprintf('  FAIL  coeff size  (%s vs %s)\n',mat2str(size(h1.coeff)),mat2str(size(h2.coeff))); cnt(2)=cnt(2)+1; end

err_c = max(abs( single(h1.coeff(:)) - single(h2.coeff(:)) ));
ok = err_c < tol_coeff;
if ok; fprintf('  PASS  coeff match\n'); cnt(1)=cnt(1)+1;
else;  fprintf('  FAIL  coeff match  (max err = %g)\n',err_c); cnt(2)=cnt(2)+1; end

err_d = max(abs( double(h1.delay(:)) - double(h2.delay(:)) ));
ok = err_d < tol_delay;
if ok; fprintf('  PASS  delay match\n'); cnt(1)=cnt(1)+1;
else;  fprintf('  FAIL  delay match  (max err = %g)\n',err_d); cnt(2)=cnt(2)+1; end
fprintf('\n');

% ============================================================================
%% Summary
% ============================================================================
fprintf('=================================================================\n');
fprintf('  Total: %d passed, %d failed\n', cnt(1), cnt(2));
if cnt(2) == 0
    fprintf('  ALL TESTS PASSED.\n');
else
    fprintf('  *** SOME TESTS FAILED ***\n');
end
fprintf('=================================================================\n');
