function testBuilder_add_clusters

%nnn=85
%RandStream.setGlobalStream(RandStream('mt19937ar','seed',nnn));

% Create layout
l = qd_layout;
l.simpar.autocorrelation_function = 'Disable';
l.simpar.center_frequency(2) = 3.7e9;
l.simpar.show_progress_bars = 0;
l.tx_position(:,2) = [100;0;25];
l.no_rx = 5;
l.randomize_rx_positions(200,1.5,1.5,0);
l.set_scenario('3GPP_38.901_UMi',[],1,0,0,0,0);
l.set_scenario('3GPP_38.901_UMi_LOS_GR',[],2,0,0,0,0);

%%
% Initialize builders
b = init_builder( l );

% Split into one builder per RX
bs = split_rx( b );

% Call add_paths on an empty array
add_paths( bs );

% Initialize SOS
init_sos( b );

% Add a ground-reflection component to all builders
add_paths( b, 'TwoRayGR' );

% Check results - results should only have a GR path
sic = size( b );
gain = {}; sf = {}; gr_epsilon_r = {};
for n = 1 : sum(sic)
    [ i1,i2 ] = qf.qind2sub( sic, n );
    gain{n} = b(i1,i2).gain;
    sf{n} = b(i1,i2).sf;
    gr_epsilon_r{n} = b(i1,i2).gr_epsilon_r;
    if b(i1,i2).no_rx_positions > 0
        assertEqual( b(i1,i2).NumClusters, 2 ); % LOS + GR
        assertEqual( b(i1,i2).NumSubPaths, [1 1] );
        assertEqual( size(b(i1,i2).subpath_coupling),...
            [4,(b(i1,i2).scenpar.NumClusters - 1 - b(i1,i2).scenpar.GR_enabled) * b(i1,i2).scenpar.NumSubPaths + 2,2] );
        assertTrue( ~isempty( b(i1,i2).gr_epsilon_r )  )
        assertTrue( all(isinf(b(i1,i2).kf(:)))  )
        assertEqual( size(b(i1,i2).gain), [b(i1,i2).no_rx_positions,2,2] );
        assertEqual( size(b(i1,i2).taus), [b(i1,i2).no_rx_positions,2] );
        assertEqual( size(b(i1,i2).fbs_pos), [3,2,b(i1,i2).no_rx_positions,b(i1,i2).no_freq] );
    end
end

% Create a NLOS-only builder, modify some parameters
ba = qd_builder('3GPP_38.901_InF_NLOS_DL');
ba.scenpar.NumClusters = 3;
ba.scenpar.PerClusterDS = 1;
ba.scenpar.PerClusterDS_gamma = 1;

% Add NLOS paths to exisiting builder
bld = add_paths( b, ba );

% Check values of the added data
assertEqual( numel(ba), 1 );
assertEqual( ba.simpar.autocorrelation_function, 'Disable' );
assertEqual( ba.no_freq, 2 );
assertEqual( ba.no_rx_positions, l.no_rx * l.no_tx );
assertEqual( ba.NumClusters, 7 ); % LOS + 2*NLOS*3 (NLOS-Split)
assertEqual( ba.NumSubPaths, [1 10 6 4 10 6 4] );
assertEqual( size(ba.taus), [ba.no_rx_positions,7,2] );     % Different delays over frequency
assertTrue( isempty( ba.gr_epsilon_r )  );

% "bld" and "ba" should be identical
bx = split_rx(ba);
assertEqual( bld(1,5).gain, bx(1,5).gain );
assertEqual( bld(1,3).taus, bx(1,3).taus );
assertEqual( bld(1,7).xprmat, bx(1,7).xprmat );
assertEqual( bld(1,2).sos(1,1).sos_phase, bx(1,2).sos(1,1).sos_phase );

% Check results - LOS and GR path should not change, NLOS should be added
cnt = 0;
lbs_pos = {};
for n = 1 : sum(sic)
    [ i1,i2 ] = qf.qind2sub( sic, n );
    lbs_pos{n} = b(i1,i2).lbs_pos;
    if b(i1,i2).no_rx_positions > 0
        assertEqual( b(i1,i2).NumClusters, 8 ); % LOS + GR + 2*NLOS
        assertEqual( b(i1,i2).NumSubPaths, [1 1 10 6 4 10 6 4] );
        assertEqual( size(b(i1,i2).subpath_coupling),...
            [4,(b(i1,i2).scenpar.NumClusters - 1 - b(i1,i2).scenpar.GR_enabled) * b(i1,i2).scenpar.NumSubPaths + 2,2] );
        assertTrue( ~isempty( b(i1,i2).gr_epsilon_r )  )
        assertTrue( all(~isinf(b(i1,i2).kf(:)))  )
        assertEqual( size(b(i1,i2).gain), [b(i1,i2).no_rx_positions,8,2] );
        assertEqual( size(b(i1,i2).taus), [b(i1,i2).no_rx_positions,8,2] );
        assertEqual( b(i1,i2).gain(:,1:2,:), gain{n} );             % GR and LOS stay the same
        assertEqual( b(i1,i2).gain(:,3:4,:), ba.gain(cnt+1:cnt+b(i1,i2).no_rx_positions,[2,3],:) ); % NLOS Paths
        assertEqual( b(i1,i2).gr_epsilon_r, gr_epsilon_r{n} );      % epsilon_r stays the same
        assertTrue( all( b(i1,i2).sf(:) > sf{n}(:) ) );             % SF must increase due to added paths
        assertEqual( size(b(i1,i2).fbs_pos), [3,42,b(i1,i2).no_rx_positions,b(i1,i2).no_freq] );
        cnt = cnt + b(i1,i2).no_rx_positions;
    end
    sf{n} = b(i1,i2).sf;
end

% Adding a ground-reflection component a second time sould replace the exisiting value, but not
% change NLOS data
add_paths( b, 'TwoRayGR' );

% Check results
cnt = 0;
for n = 1 : sum(sic)
    [ i1,i2 ] = qf.qind2sub( sic, n );
    if b(i1,i2).no_rx_positions > 0
        assertEqual( b(i1,i2).NumClusters, 8 ); % LOS + GR + 2*NLOS
        assertEqual( b(i1,i2).NumSubPaths, [1 1 10 6 4 10 6 4]  );
        assertEqual( size(b(i1,i2).subpath_coupling),...
            [4,(b(i1,i2).scenpar.NumClusters - 1 - b(i1,i2).scenpar.GR_enabled) * b(i1,i2).scenpar.NumSubPaths + 2,2] );
        assertTrue( ~isempty( b(i1,i2).gr_epsilon_r )  )
        assertTrue( all(~isinf(b(i1,i2).kf(:)))  )
        assertEqual( size(b(i1,i2).gain), [b(i1,i2).no_rx_positions,8,2] );
        assertEqual( size(b(i1,i2).taus), [b(i1,i2).no_rx_positions,8,2] );
        x = b(i1,i2).gain(:,1:2,:) - gain{n};
        assertTrue( all( x(:)~=0 ) );             % GR and LOS differ due to different epsilon_r
        assertEqual( b(i1,i2).gain(:,3:4,:), ba.gain(cnt+1:cnt+b(i1,i2).no_rx_positions,[2,3],:) ); % NLOS Paths
        x = b(i1,i2).gr_epsilon_r - gr_epsilon_r{n};
        assertTrue( all( x(:)~=0 ) );       % epsilon_r stays the same
        assertTrue( all(abs( b(i1,i2).sf(:) - sf{n}(:) ) < 1e-12));   % SF must stay the same
        assertEqual(  b(i1,i2).lbs_pos, lbs_pos{n} );
        cnt = cnt + b(i1,i2).no_rx_positions;
    end
end


% Adding a freespace component should replace the LOS and GR components, but not change NLOS data
add_paths( b, 'Freespace' );

% Check results
cnt = 0;
for n = 1 : sum(sic)
    [ i1,i2 ] = qf.qind2sub( sic, n );
    if b(i1,i2).no_rx_positions > 0
        assertEqual( b(i1,i2).NumClusters, 8 ); % LOS + GR + 2*NLOS
        assertEqual( b(i1,i2).NumSubPaths, [1 1 10 6 4 10 6 4] );
        assertEqual( size(b(i1,i2).subpath_coupling),...
            [4,(b(i1,i2).scenpar.NumClusters - 1 - b(i1,i2).scenpar.GR_enabled) * b(i1,i2).scenpar.NumSubPaths + 2,2] );
        assertTrue( ~isempty( b(i1,i2).gr_epsilon_r )  )
        assertTrue( all(~isinf(b(i1,i2).kf(:)))  )
        assertEqual( size(b(i1,i2).gain), [b(i1,i2).no_rx_positions,8,2] );
        assertEqual( size(b(i1,i2).taus), [b(i1,i2).no_rx_positions,8,2] );
        x = b(i1,i2).gain(:,2,:);
        assertTrue( all( x(:)==0 ) );             % GR must be zero
        assertEqual( b(i1,i2).gain(:,3:4,:), ba.gain(cnt+1:cnt+b(i1,i2).no_rx_positions,[2,3],:) ); % NLOS Paths
        assertTrue( all(abs( b(i1,i2).sf(:) - sf{n}(:) ) < 1e-12));   % SF must stay the same
        assertEqual(  b(i1,i2).lbs_pos, lbs_pos{n} );
        cnt = cnt + b(i1,i2).no_rx_positions;
    end
end

% Calling "add_paths" without a parameter definition should add all paths from the original builders
bx = add_paths( b );

% Check if the same SOS random generator were used
assertEqual( bx(1,1).sos(1,1).sos_phase, b(1,1).sos(1,1).sos_phase );

% Check results
cnt = 0;
for n = 1 : sum(sic)
    [ i1,i2 ] = qf.qind2sub( sic, n );
    if b(i1,i2).no_rx_positions > 0
        NumSubPaths = (b(i1,i2).scenpar.NumClusters-1-b(i1,i2).scenpar.GR_enabled)*b(i1,i2).scenpar.NumSubPaths + 2 + 40;
        try
            assertEqual( sum(b(i1,i2).NumSubPaths), NumSubPaths );
        catch
            error('FixMe')
        end
        assertEqual( size(b(i1,i2).subpath_coupling),[4,NumSubPaths,2] );
        assertTrue( ~isempty( b(i1,i2).gr_epsilon_r )  )
        assertEqual( size(b(i1,i2).gain), [b(i1,i2).no_rx_positions,b(i1,i2).NumClusters,2] );
        assertEqual( size(b(i1,i2).taus), [b(i1,i2).no_rx_positions,b(i1,i2).NumClusters,2] );
        x = b(i1,i2).gain(:,2,:);
        if b(i1,i2).scenpar.GR_enabled == 0
            assertTrue( all( x(:) < 1e-19 ) );             % GR must be zero
        end
        assertEqual( b(i1,i2).gain(:,3:4,:), ba.gain(cnt+1:cnt+b(i1,i2).no_rx_positions,[2,3],:) ); % NLOS Paths
        assertEqual( b(i1,i2).lbs_pos(:,1:42,:,:), lbs_pos{n} );
        cnt = cnt + b(i1,i2).no_rx_positions;
    end
end

% Add clusters from exisiting builder
add_paths( b, b );

% Check results
cnt = 0;
for n = 1 : sum(sic)
    [ i1,i2 ] = qf.qind2sub( sic, n );
    if b(i1,i2).no_rx_positions > 0
        NumSubPaths = 2*(b(i1,i2).scenpar.NumClusters-1-b(i1,i2).scenpar.GR_enabled)*b(i1,i2).scenpar.NumSubPaths + 2 + 40;
        assertEqual( sum(b(i1,i2).NumSubPaths), NumSubPaths );
        assertEqual( size(b(i1,i2).subpath_coupling),[4,NumSubPaths,2] );
        assertTrue( ~isempty( b(i1,i2).gr_epsilon_r )  )
        assertEqual( size(b(i1,i2).gain), [b(i1,i2).no_rx_positions,b(i1,i2).NumClusters,2] );
        assertEqual( size(b(i1,i2).taus), [b(i1,i2).no_rx_positions,b(i1,i2).NumClusters,2] );
        x = b(i1,i2).gain(:,2,:);
        if b(i1,i2).scenpar.GR_enabled == 0
            assertTrue( all( x(:) < 1e-19 ) );             % GR must be zero
        end
        assertEqual( b(i1,i2).gain(:,3:4,:), ba.gain(cnt+1:cnt+b(i1,i2).no_rx_positions,[2,3],:) ); % NLOS Paths
        assertEqual( b(i1,i2).lbs_pos(:,1:42,:,:), lbs_pos{n} );
        cnt = cnt + b(i1,i2).no_rx_positions;
    end
end

b = split_multi_freq( b );
b = split_rx( b );
c = get_channels( b );

for n = 1 : numel(c)
    assertEqual( b(1,n).NumClusters, c(1,n).no_path );
end

