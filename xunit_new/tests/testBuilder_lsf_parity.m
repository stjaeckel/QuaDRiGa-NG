function testBuilder_lsf_parity
% Path loss and large-scale parameters reproduce the stored reference values of all scenarios
%
% The reference values are written by 'data/make_lsf_parity_ref.m'. Each scenario is evaluated
% with the stored phases on a drop of 24 positions (3.5 GHz), a drop of 3 positions (3.5 and
% 28 GHz), an outdoor track and two indoor tracks.

d = load( fullfile( fileparts( mfilename('fullpath') ), 'data', 'lsf_parity_ref.mat' ) );
tol = 1e-12;

for n = 1 : numel( d.ref )
    ref = d.ref(n);
    b = qd_builder( ref.scenario );
    b.simpar.show_progress_bars = 0;

    for drop = { ref.drop, ref.drop3 }
        r = drop{1};
        msg = [ ref.scenario, ', ', num2str( size(r.rx,2) ), ' positions' ];
        setup_builder( b, d.phase, r.tx, r.rx, r.fc );
        gen_lsf_parameters( b );

        lsp = cat( 3, b.ds, b.kf, b.sf, b.asD, b.asA, b.esD, b.esA, b.xpr );
        assert_close( lsp(:,r.lsp_valid,:), r.lsp(:,r.lsp_valid,:), tol, [msg,', LSPs'] );
        if ~isempty( r.abs_toa )
            assert_close( b.absTOA_offset, r.abs_toa, tol, [msg,', abs. TOA'] );
        end
        if ~isempty( r.gr )
            assert_close( b.gr_epsilon_r, r.gr, tol, [msg,', ground permittivity'] );
        end

        [ pl, scale_sf ] = get_pl_per_freq( b );
        assert_close( pl(:,r.pl_valid), r.pl(:,r.pl_valid), tol, [msg,', path loss'] );
        if r.sf_base_valid
            sf_base = b.sf .^ scale_sf;
            assert_close( sf_base(:,r.pl_valid), r.sf_base(:,r.pl_valid), tol, [msg,', baseline SF'] );
        end
    end

    setup_builder( b, d.phase, ref.drop.tx, ref.drop.rx, ref.drop.fc );
    for k = 1 : numel( ref.trk )
        r = ref.trk(k);
        msg = [ ref.scenario, ', track ', num2str(k) ];
        t = qd_track([]);
        t.positions = r.positions;
        t.initial_position = r.initial_position;
        t.par = r.par;

        [ pl, ~, pl_init ] = b.get_pl( t );
        [ sf, kf ] = b.get_sf_profile( t );
        if k == 1                   % Outdoor track: all snapshots
            is = 1 : size( r.positions, 2 );
        else                        % Indoor track: first snapshot only
            is = 1;
        end
        ip = is( r.pl_valid(is) );
        assert_close( pl(:,ip), r.pl(:,ip), tol, [msg,', path loss'] );
        if r.pl_valid(end)
            assert_close( pl_init, r.pl_init, tol, [msg,', path loss at initial position'] );
        end
        il = is( r.lsp_valid(is) );
        assert_close( sf(:,il), r.sf(:,il), tol, [msg,', SF profile'] );
        assert_close( kf(:,il), r.kf(:,il), tol, [msg,', KF profile'] );
    end
end
end

function setup_builder( b, phase, tx, rx, fc )
% Builder with the given positions and the stored phases of the large-scale generators
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

function [ pl, scale_sf ] = get_pl_per_freq( b )
% Path loss at the builder positions, one call per frequency for the satellite model
fc = b.simpar.center_frequency;
if isstruct( b.plpar ) && isfield( b.plpar, 'model' ) && strcmp( b.plpar.model, 'satellite' )
    pl = zeros( numel(fc), b.no_rx_positions );
    scale_sf = pl;
    for k = 1 : numel( fc )
        b.simpar.center_frequency = fc(k);
        [ pl(k,:), scale_sf(k,:) ] = b.get_pl;
    end
    b.simpar.center_frequency = fc;
else
    [ pl, scale_sf ] = b.get_pl;
end
end

function assert_close( a, b, tol, msg )
assertEqual( size(a), size(b), msg );
ok = a == b | abs( a - b ) <= tol * abs( b );
assertTrue( all( ok(:) ), msg );
end
