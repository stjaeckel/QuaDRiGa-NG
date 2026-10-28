function init_sos( h_builder, usage )
%INIT_SOS Initializes all SOS random number generators
%
% Calling object:
%   Object array
%
% Description:
%   This method initializes all random number generators of the "qd_builder" object array. The
%   results are stored in the properties "lsf_generator", "path_sos", "xpr_sos", "pin_sos", and
%   "clst_dl_sos". Once initialized, the "qd_builder" object will generate identical parameters and
%   channels. Note that this does not work for 3GPP baseline simulations.
%
%   The phases of the large-scale generators are drawn by quadriga_lib.lsf_init_phases with a seed
%   from MATLAB's "rand"; the sample mean of each normal field over the Rx positions (and Tx
%   positions in dual mobility) is removed. The sum-of-sinusoids preset is given by
%   "qd_simulation_parameters.autocorrelation_function" ('Comb256' for 'Disable', which sets the
%   decorrelation distances of the ground permittivity and the absolute TOA to 0). With
%   "qd_simulation_parameters.lsf_model = 'legacy'", the large-scale parameters are correlated in 3D
%   space, otherwise in the horizontal plane.
%
% Input:
%   usage
%   Controls the behavior of the method: If set to 0, all existing SOS generators will be discarded
%   and the method exits. By default (1), new SOS generators will be created, existing ones will be
%   replaced. It set to 2, existing SOS generators will be reused and missing ones will be created.
%
% 
% QuaDRiGa Copyright (C) 2011-2020
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

% Parse input variables.
if ~exist( 'usage','var' ) || isempty( usage )
    usage = 1;
end

if numel(h_builder) > 1
    
    % Recursive call for all objects in the array
    sic = size( h_builder );
    for i_cb = 1 : numel(h_builder)
        [ i1,i2,i3,i4 ] = qf.qind2sub( sic, i_cb );
        init_sos( h_builder( i1,i2,i3,i4 ), usage );
    end
    
else
    
    % Fix for octave 4.0 (conversion from object-array to single object)
    h_builder = h_builder(1,1);
    
    % Remove existing parameters
    if usage == 0 || usage == 1
        h_builder.lsf_generator = [];
        h_builder.path_sos = [];
        h_builder.xpr_sos = [];
        h_builder.pin_sos = [];
        h_builder.clst_dl_sos = [];
        h_builder.subpath_coupling = [];
    end
    
    % Dual Mobility indicator
    if h_builder.dual_mobility == -1 
        h_builder.check_dual_mobility( false );
    end
    dual_mobility = h_builder.dual_mobility;
    
    % Check if we need to do anything else
    if usage == 0
        return
    end
    
    % Set the number of clusters
    L = h_builder.scenpar.NumClusters;
    
    % Number of paths per NLOS cluster
    M = h_builder.scenpar.NumSubPaths;
    
    % Spatial consistency decorrelation distance in [m]
    SC_lambda = h_builder.scenpar.SC_lambda;
    
    % ACF
    acf = h_builder.simpar(1,1).autocorrelation_function;
    if strcmp( acf, 'Disable' )     % Disable spatial consistency for SSF
        acf = 'Comb256';
        SC_lambda = 0;
    end
    
    use_ground_reflection = logical( h_builder.scenpar.GR_enabled );
    if use_ground_reflection && L == 1
        error('You need at least 2 clusters to enable the ground reflection option.');
    end
    
    % Number of NLOS clusters
    Ln = L - 1 - double( use_ground_reflection );
    
    % Number of frequencies
    F = numel( h_builder.simpar(1,1).center_frequency );
    
    % Random generators of the large-scale parameters
    if isempty( h_builder.lsf_generator )
        init_lsf_generator( h_builder, dual_mobility );
    end
    
    % Initialize SOS generators for the generation of multipath components
    if ~h_builder.simpar(1,1).use_3GPP_baseline && Ln ~= 0 && SC_lambda > 0 && isempty( h_builder.path_sos )
        path_sos = qd_sos([]);
        for i_cluster = 1 : Ln
            % Delays (are identical if Tx and Rx positions are swapped)
            path_sos(i_cluster,1) = qd_sos( acf, 'Uniform', SC_lambda ); % Delays
            path_sos(i_cluster,1).sos_phase(:,2) = path_sos(i_cluster,1).sos_phase(:,1);
            
            % Azimuth angles (departure become arrival angles if positions are swapped)
            path_sos(i_cluster,2) = qd_sos( acf, 'Uniform', SC_lambda ); % ASD
            path_sos(i_cluster,3) = qd_sos( acf, 'Uniform', SC_lambda ); % ASA
            path_sos(i_cluster,2).sos_phase(:,1) = path_sos(i_cluster,3).sos_phase(:,2);
            path_sos(i_cluster,3).sos_phase(:,1) = path_sos(i_cluster,2).sos_phase(:,2);
            
            % Elevation angles (departure become arrival angles if positions are swapped)
            path_sos(i_cluster,4) = qd_sos( acf, 'Uniform', SC_lambda ); % ASD
            path_sos(i_cluster,5) = qd_sos( acf, 'Uniform', SC_lambda ); % ASA
            path_sos(i_cluster,4).sos_phase(:,1) = path_sos(i_cluster,5).sos_phase(:,2);
            path_sos(i_cluster,5).sos_phase(:,1) = path_sos(i_cluster,4).sos_phase(:,2);
        end
        h_builder.path_sos = path_sos;
    end
    
    % Initialize SOS generators for the XPR
    if ~h_builder.simpar(1,1).use_3GPP_baseline && Ln ~= 0 && SC_lambda > 0 && isempty( h_builder.xpr_sos )
        xpr_sos = qd_sos([]);
        for i_freq = 1 : F
            for i_cluster = 1 : Ln
                for i_sub = 1 : M+1
                    xpr_sos(i_cluster,i_sub,i_freq) = qd_sos( acf, 'Normal', SC_lambda );
                    if i_freq == 1 % Fix for Octave
                        xpr_sos(i_cluster,i_sub).sos_phase(:,2) = xpr_sos(i_cluster,i_sub).sos_phase(:,1);
                    else
                        xpr_sos(i_cluster,i_sub,i_freq).sos_phase(:,2) = xpr_sos(i_cluster,i_sub,i_freq).sos_phase(:,1);
                    end
                end
            end
        end
        h_builder.xpr_sos = xpr_sos;
    end
    
    % Initialize SOS generators for initial phases
    if ~h_builder.simpar(1,1).use_3GPP_baseline && Ln ~= 0 && SC_lambda > 0 && isempty( h_builder.pin_sos )
        pin_sos = qd_sos([]);
        for i_freq = 1 : F
            for i_cluster = 1 : Ln
                for i_sub = 1 : M
                    pin_sos(i_cluster,i_sub,i_freq) = qd_sos( acf, 'Uniform', SC_lambda );
                    if i_freq == 1 % Fix for Octave
                        pin_sos(i_cluster,i_sub).sos_phase(:,2) = pin_sos(i_cluster,i_sub).sos_phase(:,1);
                    else
                        pin_sos(i_cluster,i_sub,i_freq).sos_phase(:,2) = pin_sos(i_cluster,i_sub,i_freq).sos_phase(:,1);
                    end
                end
            end
        end
        h_builder.pin_sos = pin_sos;
    end
    
    % Initialize SOS generators for mmMAGIC intra-cluster delay offsets
    if ~h_builder.simpar(1,1).use_3GPP_baseline && Ln ~= 0 && SC_lambda > 0 && isempty( h_builder.clst_dl_sos )
        clst_dl_sos = qd_sos([]);
        for i_cluster = 1 : Ln
            for i_sub = 1 : M
                clst_dl_sos(i_cluster,i_sub) = qd_sos( acf, 'Uniform', SC_lambda );
                clst_dl_sos(i_cluster,i_sub).sos_phase(:,2) = clst_dl_sos(i_cluster,i_sub).sos_phase(:,1);
            end
        end
        h_builder.clst_dl_sos = clst_dl_sos;
    end
    
    % Initialize subpath coupling
    NumSubPaths = Ln * M + L - Ln;
    if isempty( h_builder.subpath_coupling )
        if ~isempty( h_builder.NumSubPaths )
            NumSubPaths = sum( h_builder.NumSubPaths );
        end
        if h_builder.simpar(1,1).use_3GPP_baseline
            h_builder.subpath_coupling = rand( 4 , NumSubPaths , h_builder.no_rx_positions );
        else
            h_builder.subpath_coupling = rand( 4 , NumSubPaths , F );
        end
    elseif size( h_builder.subpath_coupling,2 ) > NumSubPaths
        h_builder.subpath_coupling = h_builder.subpath_coupling( :,1:NumSubPaths,: );
    end
end
end

function init_lsf_generator( h_builder, dual_mobility )
% Draws the phases of the large-scale generators

acf = h_builder.simpar(1,1).autocorrelation_function;
disable_sc = strcmp( acf, 'Disable' );
if disable_sc
    acf = 'Comb256';
end

% Names of previous versions map to 128, 256, 512 and 1024 sinusoids
tok = regexp( acf, '^(Exp|Gauss|Comb)(150|300|500|1000)$', 'tokens', 'once' );
if ~isempty( tok )
    new_no = {'128','256','512','1024'};
    acf = [ tok{1}, new_no{ strcmp( tok{2}, {'150','300','500','1000'} ) } ];
end

scp = h_builder.scenpar;
names = { 'DS_lambda','KF_lambda','SF_lambda','AS_D_lambda','AS_A_lambda','ES_D_lambda',...
    'ES_A_lambda','XPR_lambda','absTOA_lambda','SC_lambda' };
lambda = struct();
for n = 1 : numel( names )
    lambda.( names{n} ) = scp.( names{n} );
end
if disable_sc
    lambda.absTOA_lambda = 0;
    lambda.SC_lambda = 0;
end

horizontal = ~strcmp( h_builder.simpar(1,1).lsf_model, 'legacy' );
seed = floor( rand * 2^52 );

% Positions for the bias removal
rx_pos = [];
tx_pos = [];
if h_builder.no_rx_positions > 2
    rx_pos = h_builder.rx_positions.';
    tx_pos = h_builder.tx_position.';
end

h_builder.lsf_generator = struct( 'phase', ...
    quadriga_lib.lsf_init_phases( lambda, seed, tx_pos, rx_pos, logical( dual_mobility ), true, acf, horizontal ), ...
    'preset', acf, 'horizontal', horizontal, 'lambda', lambda );
end
