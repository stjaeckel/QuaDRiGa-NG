function gen_lsf_parameters( h_builder, usage, ~ )
%GEN_LSF_PARAMETERS Generates the large-scale-fading parameters for all terminals and frequencies
%
% Calling object:
%   Object array
%
% Description:
%   This method generates correlated large scale parameters for each user position. Those
%   parameters are needed by the channel builder to calculate initial SSF parameters for each track
%   or segment which are then evolved into time varying channels. The LSF parameters require that
%   the random generators are initialized first. Hence, you need to call 'qd_builder.init_sos'
%   before calling 'qd_builder.gen_lsf_parameters'. The output of this method is written to the
%   object properties.
%
%   The normalized random values are drawn by 'quadriga_lib.lsf_draw' and transformed into the
%   LSPs by 'quadriga_lib.lsp_generate' (3GPP TR 38.901 clause 7.5, step 4). The setting
%   'qd_simulation_parameters.lsf_model' selects the model:
%   * '38.901': Cholesky factorization of the cross-correlation matrix, azimuth spreads limited to
%     104 deg and elevation spreads to 52 deg, shadow-fading standard deviation of the path-loss
%     model
%   * 'legacy': symmetric square root of the cross-correlation matrix, no limits, shadow-fading
%     standard deviation of the scenario parameters (of the path-loss model in 3GPP baseline mode)
%   The ground permittivity and the absolute TOA offset are spatially consistent with the
%   decorrelation distances 'SC_lambda' and 'absTOA_lambda'.
%
% LSF parameters:
%   ds              The RMS delay spread in [s]
%   kf              The Ricean K-Factor [linear scale]
%   sf              The shadow fading [linear scale]
%   asD             The azimuth spread of departure in [deg]
%   asA             The azimuth spread of arrival in [deg]
%   esD             The elevation spread of departure in [deg]
%   esA             The elevation spread of arrival in [deg]
%   xpr             The cross polarization ratio [linear scale]
%   gr_epsilon_r    The relative permittivity for the ground reflection (optional)
%   absTOA_offset   The absolute time-of-arrival offset in [s] (optional)
%
% Input:
%   usage
%   Controls the behavior of the method: If set to 0, all existing LSF parameters will be discarded
%   and the method exits. By default (1), new LSF parameters will be created, existing ones will be
%   replaced. It set to 2, existing LSF parameters will be reused and missing ones will be created.
%
%   The third input is accepted for compatibility and has no effect.
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

if ~exist( 'usage','var' ) || isempty( usage )
    usage = 1;
end

if numel(h_builder) > 1
    
    % Recursive call for all objects in the array
    sic = size( h_builder );
    for i_cb = 1 : numel(h_builder)
        [ i1,i2,i3,i4 ] = qf.qind2sub( sic, i_cb );
        gen_lsf_parameters( h_builder( i1,i2,i3,i4 ), usage );
    end
    
else
    
    % Fix for octave 4.0 (conversion from object-array to single object)
    h_builder = h_builder(1,1);
    
    % Remove existing parameters
    if usage == 0 || usage == 1
        h_builder.ds = [];
        h_builder.kf = [];
        h_builder.sf = [];
        h_builder.asD = [];
        h_builder.asA = [];
        h_builder.esD = [];
        h_builder.esA = [];
        h_builder.xpr = [];
        h_builder.gr_epsilon_r = [];
        h_builder.absTOA_offset = [];
    end
    
    if usage == 0 || h_builder.no_rx_positions == 0
        return
    end
    
    if isempty( h_builder.lsf_generator )
        error('QuaDRiGa:qd_builder:gen_lsf_parameters','Random generators are not initialized.');
    end
    if ~h_builder.lsp_xcorr_chk
        error('QuaDRiGa:qd_builder:gen_lsf_parameters',...
            ['LSP cross-correlation matrix of "',h_builder.name,'" is not positive definite.']);
    end
    use_ground_reflection = logical( h_builder.scenpar.GR_enabled );
    if use_ground_reflection && h_builder.scenpar.NumClusters == 1
        error('QuaDRiGa:qd_builder:gen_parameters:gr_2_clusters',...
            'You need at least 2 clusters to enable the ground reflection option.');
    end
    if h_builder.dual_mobility == -1
        h_builder.check_dual_mobility( false );
    end
    
    [ lsp, abs_toa, gr_eps ] = lsf_eval( h_builder, h_builder.tx_position, h_builder.rx_positions, ...
        h_builder.dual_mobility );
    
    % Write the parameters that are not set yet; dimensions [ n_freq, n_rx ]
    names = { 'ds','kf','sf','asD','asA','esD','esA','xpr' };
    for n = 1 : 8
        if isempty( h_builder.( names{n} ) )
            h_builder.( names{n} ) = lsp(:,:,n).';
        end
    end
    
    if use_ground_reflection && isempty( h_builder.gr_epsilon_r )
        h_builder.gr_epsilon_r = gr_eps.';
    elseif use_ground_reflection && numel( h_builder.gr_epsilon_r ) == 1
        % A single value applies to all frequencies and positions
        h_builder.gr_epsilon_r = ones( size( gr_eps.' ) ) * h_builder.gr_epsilon_r;
    end
    
    if h_builder.scenpar.absTOA_mu > -30 && isempty( h_builder.absTOA_offset )
        h_builder.absTOA_offset = abs_toa.';
    end
end

end
