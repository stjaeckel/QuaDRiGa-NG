function make_lsf_parity_ref
%MAKE_LSF_PARITY_REF Writes the reference values of the large-scale model parity test
%
% Description:
%   Evaluates the large-scale model of QuaDRiGa (path loss, shadow-fading scaling, the 8 LSPs,
%   absolute TOA and ground permittivity) for all scenario files with one fixed set of phases and
%   stores the results in 'lsf_parity_ref.mat'. The file is read by 'testBuilder_lsf_parity'.
%
%   For each scenario:
%   - drop:   24 Rx positions at 3.5 GHz
%   - drop3:  3 Rx positions at 3.5 and 28 GHz
%   - trk(1): outdoor track with 20 snapshots starting at a drop position
%   - trk(2): indoor track (O2I) close to the BS, the minimum outdoor distance of 10 m applies
%   - trk(3): indoor track (O2I) far from the BS
%
%   The masks mark the values that the reference is expected to match:
%   - lsp_valid: no LSP has a negative standard deviation (distance, height or elevation terms)
%   - pl_valid: 'winner_los' links before the breakpoint (per-link and mean-height breakpoint)
%   - abs_toa / gr are empty if they are drawn without a spatially consistent generator
%   - sf_base_valid: scenario has no frequency offset, distance, height or elevation SF terms
%   Indoor tracks are compared at the first snapshot and at the initial position only. The O2I
%   loss of 'nlos' files is set to 0.

phase = quadriga_lib.sos_init_phases( 256, 11, 20261028 );

list = qd_builder.supported_scenarios( 0 );
set_rand_state( 7 );

ref = struct( 'scenario', list );
for n = 1 : numel( list )
    b = qd_builder( list{n} );
    b.simpar.show_progress_bars = 0;
    is_sat = isstruct( b.plpar ) && isfield( b.plpar, 'model' ) && strcmp( b.plpar.model, 'satellite' );
    if is_sat
        tx = [ 300e3; 0; 600e3 ];
    else
        tx = [ 10; -20; 25 ];
    end

    % Drop with 24 positions at one frequency
    d2d = 10.^( 1 + rand( 1, 24 ) * 2.3 );
    phi = rand( 1, 24 ) * 2 * pi;
    rx = [ tx(1) + d2d.*cos(phi); tx(2) + d2d.*sin(phi); 1.5 + 3*floor( rand(1,24)*8 ) ];
    if is_sat
        rx(1:2,:) = rx(1:2,:) - tx(1:2) * ones(1,24);
    end
    ref(n).drop = eval_drop( b, phase, tx, rx, 3.5e9 );

    % Drop with 3 positions at two frequencies
    ref(n).drop3 = eval_drop( b, phase, tx, rx(:,[2,9,17]), [3.5e9, 28e9] );

    % Tracks
    if is_sat
        start = [ rx(:,1), [ 8; 3; 1.5 ], [ 300; 150; 4.5 ] ];
    else
        start = [ rx(:,1), tx + [ 8; 3; -23.5 ], tx + [ 300; 150; -20.5 ] ];
    end
    b = setup_builder( b, phase, tx, rx, 3.5e9 );
    is_nlos = isstruct( b.plpar ) && isfield( b.plpar, 'model' ) && strcmp( b.plpar.model, 'nlos' );
    for k = 1 : 3
        t = qd_track([]);
        t.positions = [ (0:19)*0.5*cosd(30); (0:19)*0.5*sind(30); zeros(1,20) ];
        t.initial_position = start(:,k);
        r.d2d_in = 0;
        if k > 1
            r.d2d_in = 7.3;
            d3in = indoor_d3d( tx, start(:,k), r.d2d_in );
            t.par = struct( 'o2i_loss', 12.3 * ~is_nlos, 'o2i_d3din', d3in );
        end
        [ r.pl, r.scale_sf, r.pl_init, r.scale_sf_init ] = b.get_pl( t );
        [ r.sf, r.kf ] = b.get_sf_profile( t );
        r.initial_position = t.initial_position;
        r.positions = t.positions;
        r.par = t.par;
        pos = [ t.positions + t.initial_position * ones(1,20), t.initial_position ];
        d3in = 0;
        if k > 1; d3in = t.par.o2i_d3din; end
        r.pl_valid = winner_los_valid( b.plpar, 3.5e9, tx, pos, d3in );
        r.lsp_valid = lsp_sigma_valid( b, tx, pos(:,1:20) );
        ref(n).trk(k) = r;
    end
end

save( fullfile( fileparts( mfilename('fullpath') ), 'lsf_parity_ref.mat' ), 'ref', 'phase', '-v7' );
end

function r = eval_drop( b, phase, tx, rx, fc )
b = setup_builder( b, phase, tx, rx, fc );
gen_lsf_parameters( b );
r.tx = tx;
r.rx = rx;
r.fc = fc;
r.lsp = cat( 3, b.ds, b.kf, b.sf, b.asD, b.asA, b.esD, b.esA, b.xpr );  % [n_freq, n_rx, 8]
r.abs_toa = b.absTOA_offset;
if isempty( b.absTOA_sos ); r.abs_toa = []; end
r.gr = b.gr_epsilon_r;
if isempty( b.gr_sos ); r.gr = []; end
if isstruct( b.plpar ) && isfield( b.plpar, 'model' ) && strcmp( b.plpar.model, 'satellite' )
    % The satellite model of 'get_pl' supports one frequency per call
    r.pl = zeros( numel(fc), size(rx,2) );
    r.scale_sf = r.pl;
    for k = 1 : numel( fc )
        b.simpar.center_frequency = fc(k);
        [ r.pl(k,:), r.scale_sf(k,:) ] = b.get_pl;
    end
    b.simpar.center_frequency = fc;
else
    [ r.pl, r.scale_sf ] = b.get_pl;
end
r.sf_base = b.sf .^ r.scale_sf;                 % SF of the 3GPP baseline mode (linear)
s = b.scenpar;
r.sf_base_valid = ( s.SF_omega == 0 || s.SF_delta == 0 ) && s.SF_kappa == 0 && s.SF_tau == 0 && s.SF_beta == 0;
r.lsp_valid = lsp_sigma_valid( b, tx, rx );
r.pl_valid = winner_los_valid( b.plpar, fc, tx, rx, 0 );
end

function b = setup_builder( b, phase, tx, rx, fc )
b.simpar.center_frequency = fc;
b.tx_position = tx;
b.rx_positions = rx;
b.dual_mobility = 0;
gen_lsf_parameters( b, 0 );
s = b.scenpar;
lambda = [ s.DS_lambda, s.KF_lambda, s.SF_lambda, s.AS_D_lambda, s.AS_A_lambda, ...
    s.ES_D_lambda, s.ES_A_lambda, s.XPR_lambda ];
mirror = [ 1 2 3 5 4 7 6 8 ];
lsf_sos = qd_sos([]);
for g = 1 : 8
    lsf_sos(1,g) = qd_sos( 'Comb256', 'Normal', lambda(g) );
    lsf_sos(1,g).sos_phase = [ phase(:,g), phase(:,mirror(g)) ];
end
b.sos = lsf_sos;
b.absTOA_sos = [];
if s.absTOA_mu > -30 && s.absTOA_lambda > 0
    b.absTOA_sos = qd_sos( 'Comb256', 'Normal', s.absTOA_lambda );
    b.absTOA_sos.sos_phase = phase(:,[9,9]);
end
b.gr_sos = [];
if s.SC_lambda > 0
    b.gr_sos = qd_sos( 'Comb256', 'Uniform', s.SC_lambda );
    b.gr_sos.sos_phase = phase(:,[10,10]);
end
end

function valid = lsp_sigma_valid( b, tx, rx )
% Links where no LSP standard deviation is negative
lv = b.lsp_vals;
d2d = sqrt( sum( ( rx(1:2,:) - tx(1:2) * ones(1,size(rx,2)) ).^2, 1 ) );
d2d = max( d2d, 1e-4 );
h = max( tx(3), 1e-4 );
a = atan( h ./ d2d );
valid = true( 1, size(rx,2) );
for f = 1 : size( lv, 3 )
    sig = lv(:,2,f) * ones(size(d2d)) + lv(:,7,f) * log10(d2d) + lv(:,8,f) * log10(h) + lv(:,9,f) * log10(a);
    valid = valid & all( sig >= 0, 1 );
end
end

function valid = winner_los_valid( plpar, fc, tx, rx, d3in )
% 'winner_los' links before the breakpoint (mean-height and per-link breakpoint)
valid = true( 1, size(rx,2) );
if ~isstruct( plpar ) || ~isfield( plpar, 'model' ) || ~strcmp( plpar.model, 'winner_los' )
    return
end
p = plpar;
f = fc(:) / 1e9;
o = ones( 1, size(rx,2) );
d3 = sqrt( sum( ( rx - tx * o ).^2, 1 ) );
d2 = sqrt( sum( ( rx(1:2,:) - tx(1:2) * o ).^2, 1 ) ) .* ( 1 - d3in ./ d3 );
d2 = max( d2, 0.1 );
hb = max( tx(3), 1.5 ) * o;
hm = max( rx(3,:), 1.5 );
for k = 1 : numel( f )
    G = p.B1 + p.C1*log10(f(k)) + p.D1*log10(mean(hb)) + p.E1*log10(mean(hm));
    H = p.B2 + p.C2*log10(f(k)) + p.D2*log10(mean(hb)) + p.E2*log10(mean(hm));
    bp_mean = 10.^( (H-G)./( p.A1-p.A2 ) );
    G = p.B1 + p.C1*log10(f(k)) + p.D1*log10(hb) + p.E1*log10(hm);
    H = p.B2 + p.C2*log10(f(k)) + p.D2*log10(hb) + p.E2*log10(hm);
    bp_link = 10.^( (H-G)./( p.A1-p.A2 ) );
    valid = valid & d2 <= min( bp_mean, bp_link );
end
end

function d3in = indoor_d3d( tx, rx, d2d_in )
% Indoor 3D distance as computed by qd_layout.gen_o2i_loss, minimum outdoor distance 10 m
d3 = sqrt( sum( ( rx - tx ).^2 ) );
d2 = sqrt( sum( ( rx(1:2) - tx(1:2) ).^2 ) );
d2d_in = min( d2d_in, d2 );
d3in = min( d3 * d2d_in / d2, d3 );
if d3 - d3in < 10
    d3in = d3 - 10;
end
d3in = max( d3in, 0 );
end
