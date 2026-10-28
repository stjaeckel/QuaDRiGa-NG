function [ loss, scale_sf, loss_init, scale_sf_init ] = get_pl( h_builder, evaltrack, alt_plpar, txpos )
%GET_PL Path loss and shadow-fading scaling of the path-loss model
%
% Calling object:
%   Single object
%
% Description:
%   Evaluates the path-loss model of the scenario ('qd_builder.plpar', read from the
%   configuration file) with 'quadriga_lib.path_loss'. If 'evaltrack' carries an O2I penetration
%   loss ('evaltrack.par.o2i_loss'), it is added to the path loss of all positions of the track.
%   For an indoor terminal, the path-loss model uses the geometric distance between Tx and Rx
%   (3GPP TR 38.901 eq. 7.4-1). With qd_simulation_parameters.lsf_model = 'legacy', it uses the
%   outdoor distance, which is reduced by the indoor 2D distance 'evaltrack.par.o2i_d2din'
%   (minimum outdoor distance 10 m).
%
% Input:
%   evaltrack
%   A 'qd_track' object for which the PL should be calculated. If 'evaltrack' is not given, then
%   the path loss is calculated for each Rx position without O2I loss. A [3 x N] matrix of
%   positions is also accepted (no O2I loss).
%
%   alt_plpar
%   An optional alternative plpar which is used instead of 'qd_builder.plpar'.
%
%   txpos
%   The corresponding tx positions. This variable can be provided by:
%   (1) A 'qd_track' object having the same number of snapshots as the 'evaltrack';
%   (2) A 3-element vector containing the fixed x,y, and z coordinate of the tx;
%   (3) A [3 x N] matrix containing one tx position for each position on the 'evaltrack';
%   (4) The variable is unprovided or empty. In this case, the tx-positions are obtained from the
%       'qd_builder' object.
%
% Output:
%   loss
%   The path loss in [dB]; size [ n_freq, n_pos ]
%
%   scale_sf
%   Ratio of the shadow-fading standard deviation of the path-loss model and the scenario value
%   'SF_sigma + SF_delta * log10( f_GHz )'; it scales the shadow fading in [dB]. For a track, the
%   ratio is averaged over the track positions.
%
%   loss_init
%   The path loss in [dB] for the initial position of the "evaltrack" object. Returns an empty array
%   if no "evaltrack" is given.
%
%   scale_sf_init
%   "scale_sf" for the initial position of the "evaltrack" object. Returns an empty array
%   if no "evaltrack" is given.
%
%
% QuaDRiGa Copyright (C) 2011-2019
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

if numel( h_builder ) > 1
    error('QuaDRiGa:qd_builder:ObjectArray','??? "get_pl" is only defined for scalar objects.')
else
    h_builder = h_builder(1,1); % workaround for octave
end

fc = h_builder.simpar(1,1).center_frequency;
nF = numel( fc );
o2i_loss = zeros( nF,1 );
d2d_in = [];

if ~exist( 'evaltrack','var' ) || isempty( evaltrack )
    rxpos = h_builder.rx_positions;
    use_track = false;

elseif isa( evaltrack, 'qd_track' )
    % Track positions followed by the initial position (first track of an array)
    evaltrack = evaltrack(1,1);
    initial_position = evaltrack.initial_position;
    rxpos = [ evaltrack.positions + initial_position * ones(1,evaltrack.no_snapshots), initial_position ];
    use_track = true;

    trk_par = evaltrack.par;
    if ~isempty( trk_par ) && isfield( trk_par,'o2i_loss' ) && ~isempty( trk_par.o2i_loss )
        if evaltrack.no_segments > 2
            error('QuaDRiGa:qd_builder:get_pl','Too many segments in the evaltrack');
        end
        o2i_loss = trk_par.o2i_loss( :,evaltrack.no_segments,: );
        if size( o2i_loss,1 ) ~= 1
            error('QuaDRiGa:qd_builder:get_pl','O2I-loss is defined for more than on BS.');
        end
        if size( o2i_loss,3 ) ~= nF
            error('QuaDRiGa:qd_builder:get_pl','O2I-loss must be given for each frequency.');
        end
        o2i_loss = permute( o2i_loss,[3,1,2] );
        if strcmp( h_builder.simpar(1,1).lsf_model, 'legacy' ) && isfield( trk_par,'o2i_d2din' )
            d2d_in = trk_par.o2i_d2din( :,evaltrack.no_segments );
            if numel( d2d_in ) ~= 1
                error('QuaDRiGa:qd_builder:get_pl','O2I indoor distance is defined for more than on BS.');
            end
            d2d_in = d2d_in * ones( size(rxpos,2),1 );
        end
    end

else % A 3xN matrix of positions
    rxpos = evaltrack;
    use_track = false;
end
nP = size( rxpos,2 );
oP = ones( 1,nP );

if exist( 'alt_plpar','var' ) && ~isempty( alt_plpar )
    par = alt_plpar;
else
    par = h_builder.plpar;
end

% Tx positions, one column or one column per Rx position
if ~exist('txpos','var') || isempty( txpos )
    txpos = h_builder.tx_position;
    if any( size( txpos ) ~= [3,1] ) && any( size( txpos ) ~= [3,nP] )
        error('QuaDRiGa:qd_builder:get_pl','Tx-position is ambiguous');
    end

elseif isa( txpos, 'qd_track')
    txpos = txpos(1,1);
    if txpos.no_snapshots == 1 && nP > 2
        txpos = txpos.initial_position;
    elseif txpos.no_snapshots ~= nP-1
        error('QuaDRiGa:qd_builder:get_pl',...
            'Number of positions in tx track does not match number of positions in evaltrack.');
    else
        initial_position = txpos.initial_position;
        txpos = [ txpos.positions + initial_position * oP(1:end-1) , initial_position ];
    end

elseif any( size( txpos ) ~= [3,1] ) && any( size( txpos ) ~= [3,nP] )
     error('QuaDRiGa:qd_builder:get_pl','Tx-position is ambiguous.');
end

if nP == 0
    loss = zeros( nF,0 );
    scale_sf = loss;
else
    [ loss, sf_sigma ] = quadriga_lib.path_loss( par, h_builder.scenpar, fc, txpos.', rxpos.', d2d_in );
    loss = loss.' + o2i_loss * oP;
    sf_sigma = sf_sigma.';
    if use_track
        sf_sigma = mean( sf_sigma,2 ) * oP;
    end

    SF_sigma_scenpar = h_builder.scenpar.SF_sigma + h_builder.scenpar.SF_delta * log10( fc.' / 1e9 );
    if any( sf_sigma(:) ~= 0 ) && all( SF_sigma_scenpar ~= 0 )
        scale_sf = sf_sigma ./ ( SF_sigma_scenpar * oP );
    else
        scale_sf = ones( size( sf_sigma ) );
    end
end

if use_track
    loss_init       = loss(:,end);
    scale_sf_init   = scale_sf(:,end);
    loss            = loss(:,1:end-1);
    scale_sf        = scale_sf(:,1:end-1);
else
    loss_init       = [];
    scale_sf_init   = [];
end

end
