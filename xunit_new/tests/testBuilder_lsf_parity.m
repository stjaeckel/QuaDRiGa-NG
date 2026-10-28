function testBuilder_lsf_parity
% Path loss and large-scale parameters reproduce the stored reference values of all scenarios
%
% The reference values were computed with the original QuaDRiGa implementation of the large-scale
% model. Each scenario is evaluated with the stored phases on a drop of 24 positions (3.5 GHz), a
% drop of 3 positions (3.5 and 28 GHz), an outdoor track and two indoor tracks, with
% qd_simulation_parameters.lsf_model set to 'legacy'. The masks in the data exclude links with a
% negative LSP standard deviation and 'winner_los' links beyond the breakpoint. Indoor tracks are
% compared at the first snapshot, where the indoor distance of the original model applies. The
% shadow fading of the 3GPP baseline mode is compared for scenarios without frequency-offset,
% distance, height and elevation terms of the SF standard deviation ('sf_base_valid').

d = load( fullfile( fileparts( mfilename('fullpath') ), 'data', 'lsf_parity_ref.mat' ) );
tol = 1e-12;

for n = 1 : numel( d.ref )
    ref = d.ref(n);
    b = qd_builder( ref.scenario );
    b.simpar.show_progress_bars = 0;
    b.simpar.lsf_model = 'legacy';

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

        pl = b.get_pl;
        assert_close( pl(:,r.pl_valid), r.pl(:,r.pl_valid), tol, [msg,', path loss'] );
        if r.sf_base_valid
            b.simpar.use_3GPP_baseline = 1;
            gen_lsf_parameters( b );
            b.simpar.use_3GPP_baseline = 0;
            assert_close( b.sf(:,r.pl_valid), r.sf_base(:,r.pl_valid), tol, [msg,', baseline SF'] );
        end
    end

    setup_builder( b, d.phase, ref.drop.tx, ref.drop.rx, ref.drop.fc );
    for k = 1 : numel( ref.trk )
        r = ref.trk(k);
        msg = [ ref.scenario, ', track ', num2str(k) ];
        t = qd_track([]);
        t.positions = r.positions;
        t.initial_position = r.initial_position;
        if ~isempty( r.par )
            t.par = struct( 'o2i_loss', r.par.o2i_loss, 'o2i_d2din', r.d2d_in );
        end

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
names = { 'DS_lambda','KF_lambda','SF_lambda','AS_D_lambda','AS_A_lambda','ES_D_lambda',...
    'ES_A_lambda','XPR_lambda','absTOA_lambda','SC_lambda' };
lambda = struct();
for g = 1 : numel( names )
    lambda.( names{g} ) = b.scenpar.( names{g} );
end
b.lsf_generator = struct( 'phase', cat( 3, phase, phase(:,[1 2 3 5 4 7 6 8 9 10 11]) ), ...
    'preset', 'Comb256', 'horizontal', false, 'lambda', lambda );
end

function assert_close( a, b, tol, msg )
assertEqual( size(a), size(b), msg );
ok = a == b | abs( a - b ) <= tol * abs( b );
assertTrue( all( ok(:) ), msg );
end
