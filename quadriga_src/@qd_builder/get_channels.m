function h_channel = get_channels( h_builder, vb_dots, only_coeff )
%GET_CHANNELS Calculate the channel coefficients
%
% Calling object:
%   Object array
%
% Output:
%   h_channel
%   A vector of qd_channel objects
%
%
% QuaDRiGa Copyright (C) 2011-2025
% Fraunhofer-Gesellschaft zur Foerderung der angewandten Forschung e.V. acting on behalf of its
% Fraunhofer Heinrich Hertz Institute, Einsteinufer 37, 10587 Berlin, Germany
% All rights reserved.
%
% e-mail: quadriga@hhi.fraunhofer.de
%
% This file is part of QuaDRiGa.
%
% The Quadriga software is provided by Fraunhofer on behalf of the copyright holders and
% contributors "AS IS" and WITHOUT ANY EXPRESS OR IMPLIED WARRANTIES, including but not limited to
% the implied warranties of merchantability and fitness for a particular purpose.
%
% You can redistribute it and/or modify QuaDRiGa under the terms of the Software License for
% The QuaDRiGa Channel Model. You should have received a copy of the Software License for The
% QuaDRiGa Channel Model along with QuaDRiGa. If not, see <http://quadriga-channel-model.de/>.

if ~exist( 'vb_dots','var' ) || isempty( vb_dots )
    vb_dots = [];
end

if ~exist( 'only_coeff','var' ) || isempty( only_coeff )
    only_coeff = 0;
elseif numel(h_builder) > 1
    error('QuaDRiGa:qd_builder:get_channels','Coefficient only output does not support dq_builder arrays.');
end

% Array indexing is needed for Octave
verbose = h_builder(1,1,1,1).simpar(1,1).show_progress_bars;
if verbose && isempty( vb_dots )
    fprintf('Channels     [');
    vb_dots = 50;
    tStart = clock;
    show_progress = true;
else
    show_progress = false;
end
m0=0;

if numel(h_builder) > 1

    % Equally distribute the dots in the progress bar
    sic = size( h_builder );
    vb_dots = zeros( 1,numel(h_builder) );
    for i_cb = 1 : numel(h_builder)
        [ i1,i2,i3,i4 ] = qf.qind2sub( sic, i_cb );
        if verbose
            vb_dots(i_cb) = h_builder(i1,i2,i3,i4).no_rx_positions;
        else
            % Workaround for Octave
            if numel( sic ) == 4
                h_builder(i1,i2,i3,i4).simpar(1,1).show_progress_bars = false;
            elseif numel( sic ) == 3
                h_builder(i1,i2,i3).simpar(1,1).show_progress_bars = false;
            else % 2 and 1
                h_builder(i1,i2).simpar(1,1).show_progress_bars = false;
            end
        end
    end
    if verbose
        vb_dots = init_progress_dots(vb_dots);
    end
    
    % Call each builder in the builder array and concatenate the output channels
    cnt = 1;
    h_channel = qd_channel;
    for i_cb = 1 : numel(h_builder)
        [ i1,i2,i3,i4 ] = qf.qind2sub( sic, i_cb );
        if h_builder( i1,i2,i3,i4 ).no_rx_positions > 0
            tmp = h_builder( i1,i2,i3,i4 ).get_channels( vb_dots(i_cb) );
            h_channel( 1, cnt : cnt+size(tmp,2)-1 ) = tmp;
            cnt = cnt + size(tmp,2);
        end
    end

else
    % Fix for Octave (conversion from object-array to single object)
    h_builder = h_builder(1,1);

    % Check if we have a single-frequency builder
    if numel( h_builder.simpar(1,1).center_frequency ) > 1
        error('QuaDRiGa:qd_builder:get_channels','get_channels only works for single-frequency simulations.');
    end
    center_frequency = h_builder.simpar(1,1).center_frequency;
    wave_no = 2*pi/h_builder.simpar(1,1).wavelength;

    % Check if SSF parameters have been generated already
    if isempty( h_builder.taus )
        error('QuaDRiGa:qd_builder:get_channels','Small-scale fading parameters have not been generated yet.');
    end
   
    % Check if the builder is a dual-mobility builder and that the inputs are correctly formatted
    dual_mobility = h_builder.dual_mobility;
    if dual_mobility == -1
        h_builder.check_dual_mobility;
    end

    % Set initial parameters
    n_mobiles = h_builder.no_rx_positions;

    % These variables are often needed. Pre-computing them saves a lot of time
    use_absolute_delays = h_builder.simpar(1,1).use_absolute_delays;
    if only_coeff
        if any( h_builder.NumSubPaths > 1 )
            error('QuaDRiGa:qd_builder:get_channels','Sub-paths are not supported in coefficients-only mode.');
        end
        use_3GPP_baseline = false;
        use_ground_reflection = false;
    else
        use_3GPP_baseline = h_builder.simpar(1,1).use_3GPP_baseline; % logical
        use_ground_reflection = h_builder.check_los > 1.5; % logical
        if use_3GPP_baseline && use_ground_reflection
            % For 3GPP-Baseline, GR is just another path. No need for GR-drifting
            use_ground_reflection = false;
        end
    end

    % If Laplacian PAS is used, the intra-cluster angles are increased by a factor of sqrt(2). To
    % compensate, the intra-cluster powers must be adjusted. This is done by weighting the path
    % amplitudes, depending on the number of subpaths per cluster. The weights are set here.
    if ~only_coeff && strcmp( h_builder.scenpar.SubpathMethod, 'Laplacian' )
        use_laplacian_pas = true;
        laplacian_weights = {1, [1.18 0.78], ...
            [0.60 1.05 1.24], ...
            [0.86 0.40 1.71 0.42], ...
            [1.05 0.45 0.85 0.85 1.50], ...
            [0.65 0.59 1.87 0.39 1.17 0.46], ...
            [0.77 0.75 1.06 0.61 1.02 0.53 1.74], ...
            [0.74 0.79 1.07 0.51 0.94 0.57 1.52 1.38], ...
            [0.96 0.62 0.92 0.75 1.08 0.51 1.71 1.24 0.63], ...
            [0.90 0.71 0.96 0.78 1.04 0.51 1.52 1.14 0.57 1.37], ...
            [0.92 0.66 0.89 0.79 1.03 0.52 1.47 1.13 0.60 1.36 1.15], ...
            [0.84 0.71 0.91 0.69 1.03 0.54 1.53 1.17 0.50 1.37 1.13 1.01], ...
            [0.79 0.67 1.01 0.61 0.96 0.49 1.61 1.24 0.54 1.33 1.08 1.13 0.86], ...
            [0.98 0.71 0.91 0.84 1.15 0.53 1.45 1.12 0.71 1.43 1.25 1.01 0.83 0.47], ...
            [0.99 0.70 0.89 0.83 1.09 0.57 1.32 1.11 0.64 1.37 1.16 1.01 0.80 0.53 1.41], ...
            [0.96 0.75 0.93 0.82 1.02 0.59 1.35 1.07 0.70 1.30 1.16 0.99 0.82 0.45 1.36 1.18], ...
            [0.95 0.70 0.92 0.80 1.04 0.58 1.28 1.07 0.71 1.28 1.10 1.01 0.79 0.50 1.32 1.17 1.25], ...
            [0.89 0.83 0.97 0.82 1.02 0.67 1.30 1.16 0.69 1.17 1.05 1.10 0.89 0.58 1.26 1.25 1.29 0.53], ...
            [0.91 0.79 1.01 0.85 0.99 0.66 1.29 1.14 0.70 1.20 1.04 1.03 0.92 0.56 1.23 1.21 1.27 0.53 1.15], ...
            [0.89 0.79 0.98 0.83 0.97 0.71 1.27 1.11 0.67 1.15 1.08 1.05 0.89 0.50 1.21 1.24 1.27 0.59 1.13 1.15]};
    else
        use_laplacian_pas = false;
    end

    % Create new channel object
    if only_coeff
        h_channel = zeros(h_builder.rx_array(1,1).no_elements, ...
            h_builder.tx_array(1,1).no_elements, ...
            h_builder.NumClusters, h_builder.no_rx_positions );
    else
        h_channel = qd_channel;
    end

    % Reference to the current antenna objects
    tx_reference = qd_arrayant([]);
    rx_reference = qd_arrayant([]);

    % The loop for each user position
    for i_mobile = 1 : n_mobiles

        % Update progress bar
        if verbose
            m1=ceil(i_mobile/n_mobiles*vb_dots);
            if m1>m0
                for m2=1:m1-m0
                    fprintf('o');
                end
                m0=m1;
            end
        end

        % Check if the positions are correct
        if any( abs( h_builder.rx_track(1,i_mobile).initial_position - h_builder.rx_positions(:,i_mobile) ) > 1e-5 )
            if i_mobile > 1 && ~qf.eqo( h_builder.rx_track(1,i_mobile), h_builder.rx_track(1,1) )
                warning('QuaDRiGa:qd_builder:get_channels',...
                    'Rx position in track does not match the initial position, using initial position.');
            end
            h_builder.rx_track(1,i_mobile).initial_position = h_builder.rx_positions(:,i_mobile);
        end
        if any( abs( h_builder.tx_track(1,i_mobile).initial_position - h_builder.tx_position(:,i_mobile) ) > 1e-5 )
            if i_mobile > 1 && ~qf.eqo( h_builder.tx_track(1,i_mobile), h_builder.tx_track(1,1) )
                warning('QuaDRiGa:qd_builder:get_channels',...
                    'Tx position in track does not match the initial position, using initial position.');
            end
            h_builder.tx_track(1,i_mobile).initial_position = h_builder.tx_position(:,i_mobile);
        end

        % Extract TX antenna data
        if ~qf.eqo( h_builder.tx_array(1,i_mobile), tx_reference )
            tx_reference = h_builder.tx_array(1,i_mobile);
            tx_ant = struct( ...
                'e_theta_re', real( h_builder.tx_array(1,i_mobile).Fa ), ...
                'e_theta_im', imag( h_builder.tx_array(1,i_mobile).Fa ), ...
                'e_phi_re', real( h_builder.tx_array(1,i_mobile).Fb ), ...
                'e_phi_im', imag( h_builder.tx_array(1,i_mobile).Fb ), ...
                'azimuth_grid', h_builder.tx_array(1,i_mobile).azimuth_grid, ...
                'elevation_grid', h_builder.tx_array(1,i_mobile).elevation_grid, ...
                'element_pos', h_builder.tx_array(1,i_mobile).element_position, ...
                'coupling_re', real( h_builder.tx_array(1,i_mobile).coupling ), ...
                'coupling_im', imag( h_builder.tx_array(1,i_mobile).coupling ) );
        end

        % Extract RX antenna data
        if ~qf.eqo( h_builder.rx_array(1,i_mobile), rx_reference )
            rx_reference = h_builder.rx_array(1,i_mobile);
            rx_ant = struct( ...
                'e_theta_re', real( h_builder.rx_array(1,i_mobile).Fa ), ...
                'e_theta_im', imag( h_builder.rx_array(1,i_mobile).Fa ), ...
                'e_phi_re', real( h_builder.rx_array(1,i_mobile).Fb ), ...
                'e_phi_im', imag( h_builder.rx_array(1,i_mobile).Fb ), ...
                'azimuth_grid', h_builder.rx_array(1,i_mobile).azimuth_grid, ...
                'elevation_grid', h_builder.rx_array(1,i_mobile).elevation_grid, ...
                'element_pos', h_builder.rx_array(1,i_mobile).element_position, ...
                'coupling_re', real( h_builder.rx_array(1,i_mobile).coupling ), ...
                'coupling_im', imag( h_builder.rx_array(1,i_mobile).coupling ) );
        end

        % Get the list of zero-power paths - we do not return paths with zero-power
        iClst       = h_builder.pow(i_mobile,:) ~= 0;               % Index-list of active clusters
        iClst(1)    = true;                                         % LOS is always present
        if use_ground_reflection
            iClst(2) = true;
        end
        iPath       = clst_expand( iClst, h_builder.NumSubPaths );  % Index-list of active (sub-)paths
        n_clusters  = sum( iClst );                                 % Number of clusters
        n_paths     = sum( iPath );                                 % Number od (sub-)paths
        n_subpaths  = h_builder.NumSubPaths( iClst );               % Number of paths per clusters (vector)
        n_snapshots = h_builder.rx_track(1,i_mobile).no_snapshots;  % Number of snapshots along the track segment
        n_tx_ports  = size( tx_ant.coupling_re, 2 );                % Number of transmit antenna ports
        n_rx_ports  = size( rx_ant.coupling_re, 2 );                % Number of transmit antenna ports
        n_links     = n_tx_ports * n_rx_ports;                      % Number of MIMO links in the output
        initial_pos = h_builder.rx_track(1,i_mobile).segment_index( min( [h_builder.rx_track(1,i_mobile).no_segments,2] ));

        % Extract orientation data
        rx_orientation = h_builder.rx_track(1,i_mobile).orientation;
        tx_orientation = h_builder.tx_track(1,i_mobile).orientation;
        if size( tx_orientation,2 ) == 1
            tx_orientation = repmat( tx_orientation, 1, n_snapshots);
        end

        % Extract polarization coupling Matrix in interleaved form
        M = zeros( 8, n_paths );
        M([1,3,5,7],:) = real( h_builder.xprmat(:,iPath,i_mobile) );
        M([2,4,6,8],:) = imag( h_builder.xprmat(:,iPath,i_mobile) );

        % Set the initial path gain to 1, the correct values are applied later
        path_gain = ones( 1, n_paths );

        % Extract the random initial phases
        pin = double( h_builder.pin( i_mobile, iPath ) );

        % Placeholder for the coefficient calculation
        coeff = zeros( n_links, n_clusters, n_snapshots );   % Coefficients
        delay = zeros( n_links, n_clusters, n_snapshots );   % Delays
        ppat  = zeros( n_links, n_clusters, n_snapshots );   % Radiated power

        % Extract TX and RX position
        tx_position = h_builder.tx_position(:,i_mobile);
        rx_position = h_builder.rx_positions(:,i_mobile);

        % There are 2 code paths - Planar waves (3GPP baseline) and Spherical waves
        if use_3GPP_baseline % Planar waves

            % Get the angles of the subpaths and perform random coupling.
            % Remove values that have 0-Power using "iPath" and "iClst"
            [ aod,eod,aoa,eoa ] = get_subpath_angles( h_builder, i_mobile, use_laplacian_pas );
            aod = double( aod(:,iPath) ); eod = double( eod(:,iPath) );
            aoa = double( aoa(:,iPath) ); eoa = double( eoa(:,iPath) );

            % Get the RX position relative to the track start
            tmp = h_builder.rx_track(1,i_mobile).positions;
            dist_rx = sqrt( sum([ tmp(1,:) - tmp(1,1) ; tmp(2,:) - tmp(2,1) ; tmp(3,:) - tmp(3,1) ].^2 ) );

        else % Spherical waves
            % Calculate the scatterer positions for each sub-path
            fbs_pos = double( h_builder.fbs_pos(:,iPath,i_mobile) );
            lbs_pos = double( h_builder.lbs_pos(:,iPath,i_mobile) );
        end

        % In case of ground reflection, store the exact reflection points for each snapshot
        if use_ground_reflection
            gr_pos = zeros( 3, n_snapshots );
        end

        for i_snapshot = 1 : n_snapshots

            if ~use_3GPP_baseline  % Spherical waves - update TX and RX position
                rx_position = h_builder.rx_positions(:,i_mobile) + h_builder.rx_track(1,i_mobile).positions(:,i_snapshot);
                if h_builder.tx_track(1,i_mobile).no_snapshots > 1
                    tx_position = h_builder.tx_position(:,i_mobile) + h_builder.tx_track(1,i_mobile).positions(:,i_snapshot);
                end

                % Set FBS/LBS of the LOS path to the half-way point between TX and RX
                fbs_pos(:,1) = tx_position + 0.5 * (rx_position - tx_position);
                lbs_pos(:,1) = fbs_pos(:,1);
            end

            % Calculate the total path length for each sub-path
            if i_snapshot == 1 || ~use_3GPP_baseline
                dist_rx_tx = sqrt(sum( ( rx_position - tx_position ).^2));
                path_length = double( clst_expand( dist_rx_tx + h_builder.taus( i_mobile, iClst ) * 299792458, n_subpaths ) );
            end

            % Special case "ground_reflection"
            if use_ground_reflection
                % Update the FBS/LBS position to the new reflection point
                r = [rx_position(1:2);-rx_position(3)] - tx_position;   % Direction vector
                t = -tx_position(3) / r(3);                             % Intersect parameter

                fbs_pos(1,2) = tx_position(1) + t * r(1);               % Update FBS position
                fbs_pos(2,2) = tx_position(2) + t * r(2);
                fbs_pos(3,2) = 0;
                lbs_pos(:,2) = fbs_pos(:,2);                            % Update LBS position

                gr_pos(:,i_snapshot) = fbs_pos(:,2);                    % Store ground reflection position

                path_length(2) = sqrt(sum(r.^2));                       % Update path length
                theta_r = acos( r(3) / path_length(2) ) - pi/2;         % Angle between plane and GR

                epsilon_r = h_builder.gr_epsilon_r( i_mobile );         % Relative permittivity
                Z         = sqrt( epsilon_r - (cos(theta_r)).^2 );
                R_par     = (epsilon_r * sin(theta_r) - Z) ./ (epsilon_r * sin(theta_r) + Z);
                R_per     = ( sin(theta_r) - Z) ./ ( sin(theta_r) + Z);

                % Read the path power scaling that was used in "generate_initial_paths.m"
                P_LOS = h_builder.pow(i_mobile,1);
                P_GR  = h_builder.pow(i_mobile,2);
                if P_GR < 1e-10
                    Sl  = 1;
                    gSv = 0;
                    gSh = 0;
                else
                    % Compensate for the power scaling in "generate_initial_paths.m"
                    Rsq   = 2 * P_GR / (P_LOS + P_GR);
                    gSv = sqrt(2/Rsq) * R_par;       % GR path vertical pol.
                    gSh = sqrt(2/Rsq) * R_per;       % GR path horizontal pol.
                    if P_LOS < 1e-10
                        Sl = 0;
                    else
                        Sl = 1 / sqrt( 1-Rsq/2 );    % LOS path
                    end
                end

                % Update Polarization coupling matrix
                M(1,1) = Sl;
                M(7,1) = -Sl;
                M(1,2) = real(gSv);
                M(2,2) = imag(gSv);
                M(7,2) = -real(gSh);
                M(8,2) = -imag(gSh);
            end

            % Get the MIMO channel coefficients from "quadriga_lib", outputs have size [n_rx_ports, n_tx_ports, n_paths]
            if ~use_3GPP_baseline % Spherical waves - claculate coefficients for each snapshot
                [coeff_re, coeff_im, delay_re] = quadriga_lib.get_channels_spherical( tx_ant, rx_ant, ...
                    fbs_pos, lbs_pos, path_gain, path_length, M, ...
                    tx_position, tx_orientation(:,i_snapshot), rx_position, rx_orientation(:,i_snapshot), center_frequency, use_absolute_delays );

            elseif i_snapshot == 1 % Planar waves - claculate coefficients only once
                [coeff_re, coeff_im, delay_re, rx_Doppler] = quadriga_lib.get_channels_planar( tx_ant, rx_ant, ...
                    aod, eod, aoa, eoa, path_gain, path_length, M, ...
                    tx_position, tx_orientation(:,1), rx_position, rx_orientation(:,1), center_frequency, use_absolute_delays );
            end

            % Reshape data to [ n_rx_ports * n_tx_ports, n_paths ] and apply random initial phases
            if i_snapshot == 1 || ~use_3GPP_baseline
                delay_re = reshape( delay_re, n_links, n_paths );
                cp = complex( reshape( coeff_re, n_links, n_paths ), reshape( coeff_im, n_links, n_paths ) ) .* ...
                    repmat( exp(-1j*pin), n_links, 1 );
            end

            if use_3GPP_baseline % Planar waves - Generate rotating Dopplers
                Doppler = exp( 1j * wave_no * rx_Doppler * dist_rx(i_snapshot) );
                ccp = cp .* repmat( Doppler, n_links, 1 );
            else % Spherical waves - Doppler already included in coefficients
                ccp = cp;
            end

            % Sum over the sub-paths in a cluster
            ls = 1;
            for l = 1 : n_clusters
                le = ls + n_subpaths(l) - 1;
                if le ~= ls
                    if use_laplacian_pas
                        tmp = ccp(:,ls:le) .* repmat( laplacian_weights{n_subpaths(l)}, n_links, 1 );
                        ppat(:,l,i_snapshot)  = sum( abs(tmp).^2,2 );
                        coeff(:,l,i_snapshot) = sum( tmp,2 );
                    else
                        ppat(:,l,i_snapshot)  = sum( abs(ccp(:,ls:le)).^2,2 );
                        coeff(:,l,i_snapshot) = sum( ccp(:,ls:le),2 );
                    end
                    delay(:,l,i_snapshot) = sum( delay_re(:,ls:le),2 ) / n_subpaths(l);
                else
                    ppat(:,l,i_snapshot)  = abs(ccp(:,ls)).^2;
                    coeff(:,l,i_snapshot) = ccp(:,ls);
                    delay(:,l,i_snapshot) = delay_re(:,ls);
                end
                ls = le + 1;
            end
        end

        % The path powers
        p_cl = repmat( h_builder.pow(i_mobile,iClst), n_links, 1 );

        % The powers of the antenna patterns at the given angles (power-sum)
        p_pat = sum( ppat,3 ) ./ size(ppat,3);

        % The powers in the current channel coefficients (complex sum)
        p_coeff = sum( abs(coeff).^2, 3 ) ./ size(coeff,3);

        % Correct the powers
        p_correct = sqrt( p_cl .* p_pat ./ p_coeff ./ repmat( n_subpaths, n_links, 1 ) );
        p_correct( p_pat < 1e-30 ) = 0; % Fix NaN caused by 0/0
        coeff = coeff .* repmat( p_correct, [1,1,n_snapshots] );

        % Now we apply the K-Factor and the shadowing profile
        if use_3GPP_baseline || isempty( h_builder.sos )

            % Get the PL for the initial position only
            [ ~, ~, path_loss , scale_sf ] = h_builder.get_pl( h_builder.rx_track(1,i_mobile),...
                [],h_builder.tx_track(1,i_mobile) );
            rx_power = -path_loss + 10*log10( h_builder.sf(1,i_mobile) ) .* scale_sf;
            rx_power = sqrt( 10.^( 0.1 * rx_power ) );

            % The initial KF is already applied in path powers. Here,
            % we only need to apply the SF and the path loss.
            coeff = coeff * rx_power;
        else

            % Calculate the path gain along the track segment
            path_gain = -h_builder.get_pl( h_builder.rx_track(1,i_mobile), [], h_builder.tx_track(1,i_mobile) );

            % We have a dynamic SF and KF profile that varies over the positions on the track.
            % Get shadowing profile along the track from the SOS generators.
            [sf,kf] = h_builder.get_sf_profile( h_builder.rx_track(1,i_mobile), h_builder.tx_track(1,i_mobile) );

            % When changing the cluster powers (e.g., by "add_paths"), the SF changes as well. We
            % obtain the difference of the SF by readig the initial SF values from the builder and
            % the dynamic SF profile.
            sf_init_builder = h_builder.sf(1,i_mobile);
            sf_init_sos     = sf( initial_pos );

            % We now correct the dynamic SF values (linear scale).
            sf = sf .* sf_init_builder/sf_init_sos;

            % Calculate the Rx power (sum-power of all clusters)
            rx_power = path_gain + 10*log10( sf );
            rx_power = 10.^( 0.1 * rx_power );
            rx_power = permute( rx_power, [1,3,2] );

            % Get the KF scaling
            if kf( initial_pos ) < 1e-10
                kf_scale = ones( size( kf ) );
            else
                kf_scale = kf ./ kf( initial_pos );
            end
            kf_scale = permute( kf_scale, [1,3,2] );

            % Get the normalized power for the LOS and NLOS componenets ( p_los + p_nlos = 1 )
            if use_ground_reflection
                p_los  = h_builder.pow( i_mobile,1 ) + h_builder.pow( i_mobile,2 );
                p_nlos = sum( h_builder.pow(i_mobile,3:end) );
            else
                p_los  = h_builder.pow( i_mobile,1 );
                p_nlos = sum( h_builder.pow(i_mobile,2:end) );
            end

            % Adjust the path powers to apply the varying KF along the track segment
            if p_los > 1e-4 && p_nlos > 1e-4 && any(kf_scale(:) ~= 1)

                % Adjust the power of the LOS component to match the target KF
                coeff(:,1,:) = coeff(:,1,:) .* repmat( sqrt(kf_scale),[n_links,1,1] );
                if use_ground_reflection
                    coeff(:,2,:) = coeff(:,2,:) .* repmat( sqrt(kf_scale),[n_links,1,1] );
                end

                % The power adjustment of the LOS component changes the total RX power
                % This needs to be compensated in the total RX power
                rx_power = rx_power ./ ( p_los .* kf_scale + p_nlos );
            end

            % Adjust the overall power of the channel coefficients
            rx_power = sqrt( rx_power );
            coeff = coeff .* repmat( rx_power, [n_links,n_clusters,1] );
        end

        % Save channels
        if use_3GPP_baseline % Planar waves - Average delays over all MIMO links
            delay = sum( delay,1 ) / n_links;
            delay = reshape( delay, n_clusters, n_snapshots );
        else % Spherical waves
            delay = reshape( delay , n_rx_ports , n_tx_ports , n_clusters , n_snapshots );
        end
        coeff = reshape( coeff , n_rx_ports , n_tx_ports , n_clusters , n_snapshots );
        h_channel(1,i_mobile) = qd_channel( coeff , delay , initial_pos );

        % This is important because the merger uses the name string to connect the channels.
        channel_name = h_builder.name;
        if isempty( channel_name ) || isempty( regexp( channel_name , '_', 'once' ) )
            channel_name = 'Scen_*';
        end
        % The "*" is added when there are multiple tx positions in the builder
        tmp = regexp( channel_name , '\*' );
        if ~isempty( tmp )
            channel_name = [ channel_name(1:tmp-1), h_builder.tx_track(1,i_mobile).name ];
        end
        channel_name = [ channel_name ,'_', h_builder.rx_track(1,i_mobile).name ]; %#ok
        h_channel(1,i_mobile).name = channel_name;
        h_channel(1,i_mobile).rx_position = h_builder.rx_track(1,i_mobile).positions_abs;
        h_channel(1,i_mobile).tx_position = h_builder.tx_track(1,i_mobile).positions_abs;
        h_channel(1,i_mobile).rx_orientation = h_builder.rx_track(1,i_mobile).orientation;
        h_channel(1,i_mobile).tx_orientation = h_builder.tx_track(1,i_mobile).orientation;
        h_channel(1,i_mobile).center_frequency = h_builder.simpar(1,1).center_frequency(1,1);

        % Add metadata to channel object
        clear par_struct
        if use_ground_reflection
            par_struct.has_ground_reflection = 1;
            par_struct.gr_pos = gr_pos;
        end
        par_struct.ds_parset = h_builder.ds( i_mobile ); % [s]
        par_struct.kf_parset = 10*log10( h_builder.kf( i_mobile ) ); % [db]
        if use_3GPP_baseline
            par_struct.pg_parset = 10*log10( rx_power.^2 ); % [db]
        else
            par_struct.pg_parset = 10*log10( mean(rx_power(:)).^2 ); % [db]
            par_struct.pg = 10*log10(abs( reshape( mean(mean(rx_power,1),2) , 1,[] ) ).^2);
        end
        par_struct.sf_parset = 10*log10( h_builder.sf( i_mobile ));
        par_struct.asD_parset = h_builder.asD( i_mobile ); % [deg]
        par_struct.asA_parset = h_builder.asA( i_mobile ); % [deg]
        par_struct.esD_parset = h_builder.esD( i_mobile ); % [deg]
        par_struct.esA_parset = h_builder.esA( i_mobile ); % [deg]
        if ~isempty(h_builder.xpr)
            par_struct.XPR_parset = 10*log10( h_builder.xpr( i_mobile ) ); % [db]
        end

        % Save the individual per-path values
        par_struct.AoD_cb = h_builder.AoD( i_mobile,iClst ) * 180/pi; % [deg]
        par_struct.AoA_cb = h_builder.AoA( i_mobile,iClst ) * 180/pi; % [deg]
        par_struct.EoD_cb = h_builder.EoD( i_mobile,iClst ) * 180/pi; % [deg]
        par_struct.EoA_cb = h_builder.EoA( i_mobile,iClst ) * 180/pi; % [deg]
        par_struct.pow_cb = h_builder.pow( i_mobile,iClst );          % [W]
        par_struct.gain_cb = h_builder.gain( i_mobile,iClst );          % [W]

        % Calculate the spreads at the output of the builder
        par_struct.ds_cb  = qf.calc_delay_spread( h_builder.taus( i_mobile,iClst ), h_builder.pow( i_mobile, iClst ) );
        par_struct.asD_cb = qf.calc_angular_spreads( h_builder.AoD( i_mobile,iClst ), h_builder.pow( i_mobile, iClst ) ) * 180/pi;
        par_struct.asA_cb = qf.calc_angular_spreads( h_builder.AoA( i_mobile,iClst ), h_builder.pow( i_mobile, iClst ) ) * 180/pi;
        par_struct.esD_cb = qf.calc_angular_spreads( h_builder.EoD( i_mobile,iClst ), h_builder.pow( i_mobile, iClst ) ) * 180/pi;
        par_struct.esA_cb = qf.calc_angular_spreads( h_builder.EoA( i_mobile,iClst ), h_builder.pow( i_mobile, iClst ) ) * 180/pi;

        if ~use_3GPP_baseline
            par_struct.NumSubPaths = h_builder.NumSubPaths(1,iClst);
            par_struct.fbs_pos = fbs_pos;
            par_struct.lbs_pos = lbs_pos;
        end

        % Save update rate
        mp = h_builder.rx_track(1,i_mobile).movement_profile;
        if n_snapshots > 1 && all(size(mp)==[2,2]) && mp(1,1)==0 && mp(2,1)==0 && ...
                abs(mp(2,2)-get_length(h_builder.rx_track(1,i_mobile))) < 1e-6
            par_struct.update_rate = mp(1,2)/n_snapshots;
        end
        h_channel(1,i_mobile).par = par_struct;
    end
end

% Fix for octave
if ~only_coeff && numel( h_channel ) == 1
    h_channel = h_channel(1,1);
end

if verbose && nargin == 1
    fprintf('] %5.0f seconds\n',round( etime(clock, tStart) ));
end

end
