function o2i_loss_dB = gen_o2i_loss( h_layout, method, low_loss_fraction, SC_lambda, max_indoor_dist, gen_per_tx )
%GEN_O2I_LOSS Generates the outdoor-to-indoor penetration loss
%
% Calling object:
%   Single object
%
% Description:
%   This method generates the outdoor-to-indoor (O2I) penetration loss in [dB] and the indoor 2D
%   distance of each terminal. The O2I-loss is specific for each terminal. In other words, if a MT
%   is served by several BSs, the same O2I-loss applies for each BSs. The values therefore need to
%   be generated before any other LSF or SSF parameters are generated. The automatic channel
%   generation in 'qd_layout.get_channels' will call all functions in the correct order. The
%   values are stored with the tracks in 'qd_layout.rx_track.par.o2i_loss' and
%   'qd_layout.rx_track.par.o2i_d2din'. These two variables are then read again by
%   'qd_builder.get_pl' and applied to the path-loss in the specific scenario. Each segment of a
%   track gets its own values.
%
%   The models of 3GPP TR 38.901 v19.5.0, clause 7.4.3 are evaluated by 'quadriga_lib.o2i_loss'.
%   The random variables are spatially consistent with the decorrelation distances of TR 38.901
%   clause 7.6.3.3: 50 m for the choice of the model (building type), 25 m for the indoor
%   distance and 10 m for the deviation sigma_P.
%
% Input:
%   method
%   The O2I-loss model:
%   * '3GPP_38.901': the low-loss and high-loss models of TR 38.901 Table 7.4.3-2, mixed by
%     'low_loss_fraction'
%   * One or more models of 'quadriga_lib.o2i_loss' ('low', 'high', 'low_A', 'legacy_6GHz', 'car',
%     'car_metal'); string or cell array of strings, mixed by the weights 'low_loss_fraction'
%   * 'mmMAGIC': the model of H2020-ICT-671650-mmMAGIC/D2.2, Sec. 4.3
%
%   low_loss_fraction
%   For '3GPP_38.901' and 'mmMAGIC': the likelihood of the low-loss model, from 0 (high-loss only)
%   to 1 (low-loss only). Default-value is 0.5. For a list of models: the probability of each
%   model; default: equal probabilities.
%
%   SC_lambda
%   Decorrelation distance in [m] of the random variables. Empty (default): the distances of TR
%   38.901 clause 7.6.3.3 ('mmMAGIC': 50 m). 0: uncorrelated random variables. Any other value
%   applies to all random variables.
%
%   max_indoor_dist
%   The maximum indoor distance between the building wall and the UE. Default is 25 m.
%
%   gen_per_tx
%   Logical variable to enable or disable (default) per-BS parameter generation. By default, the
%   same O2I-loss applies for each BS. However, this can be changed so that the O2I-loss is
%   generated independently for each BS by setting this variable to 'true'.
%
% Output:
%   o2i_loss_dB
%   A cell array containing the results for each MT. These values are identical to the ones stored
%   in 'qd_layout.rx_track.par.o2i_loss'. Each cell contains an array of values, where the dimensions
%   correspond to: [ BS, Segment, Frequency ]. For outdoor-scenarios, an empty array is returned.
%
%
% QuaDRiGa Copyright (C) 2011-2021
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

if numel( h_layout ) > 1
   error('QuaDRiGa:qd_layout:gen_o2i_loss','gen_o2i_loss not defined for object arrays.');
else
    h_layout = h_layout(1,1); % workaround for octave
end

if ~exist( 'method', 'var' ) || isempty( method )
    error('QuaDRiGa:qd_layout:gen_o2i_loss','O2I loss method is not defined.');
end
if ~exist( 'low_loss_fraction', 'var' )
    low_loss_fraction = [];
end
if ~exist( 'SC_lambda', 'var' )
    SC_lambda = [];
end
if ~exist( 'max_indoor_dist', 'var' ) || isempty( max_indoor_dist )
    max_indoor_dist = 25;
end
if ~exist( 'gen_per_tx', 'var' ) || isempty( gen_per_tx )
    gen_per_tx = false;
end

no_rx   = h_layout.no_rx;
no_tx   = h_layout.no_tx;
fc      = h_layout.simpar(1,1).center_frequency;
no_freq = numel( fc );

% Segment start positions of all MTs
[ pos, pos_tx, rx_ind ] = parse_positions( h_layout );
no_seg = size( pos,2 );

if ischar( method ) && strcmp( method, 'mmMAGIC' )
    if isempty( low_loss_fraction )
        low_loss_fraction = 0.5;
    end
    if isempty( SC_lambda )
        SC_lambda = 50;
    end
    [ PL, d_2d_in ] = o2i_mmmagic( fc/1e9, low_loss_fraction, SC_lambda, max_indoor_dist, gen_per_tx, pos, pos_tx );

else
    if ischar( method ) && strcmp( method, '3GPP_38.901' )
        if isempty( low_loss_fraction )
            low_loss_fraction = 0.5;
        end
        models = { 'low', 'high' };
        weights = [ low_loss_fraction, 1-low_loss_fraction ];
    else
        models = cellstr( method );
        if isempty( low_loss_fraction )
            weights = ones( 1,numel(models) ) / numel(models);
        else
            weights = low_loss_fraction;
        end
    end

    PL = zeros( no_tx, no_seg, no_freq );
    d_2d_in = zeros( no_tx, no_seg );
    if no_seg > 0
        if gen_per_tx
            n_gen = no_tx;
        else
            n_gen = 1;
        end
        for t = 1 : n_gen
            [ loss, d2d_in ] = quadriga_lib.o2i_loss( fc, o2i_random_values( pos, SC_lambda ), ...
                models, weights, max_indoor_dist );
            if gen_per_tx
                PL( t,:,: ) = reshape( loss, 1, no_seg, no_freq );
                d_2d_in( t,: ) = d2d_in.';
            else
                PL = repmat( reshape( loss, 1, no_seg, no_freq ), [no_tx,1,1] );
                d_2d_in = ones( no_tx,1 ) * d2d_in.';
            end
        end
    end
end

% Write the values to the tracks
o2i_loss_dB = cell( 1, no_rx );
for r = 1 : no_rx
    pl_tmp = PL( :, rx_ind(1,:) == r, : );
    o2i_loss_dB{ r } = pl_tmp;
    if any( abs( pl_tmp(:) ) ~= 0 )
        par_tmp = h_layout.rx_track(1,r).par;
        par_tmp.o2i_loss  = pl_tmp;
        par_tmp.o2i_d2din = d_2d_in( :, rx_ind(1,:) == r );
        h_layout.rx_track(1,r).par = par_tmp;
    end
end

end

function rnd = o2i_random_values( pos, SC_lambda )
% Random inputs of quadriga_lib.o2i_loss: model choice, 2 x indoor distance, deviation; [n_pos, 4]
n_pos = size( pos,2 );
if isempty( SC_lambda )
    lambda = [ 50, 25, 25, 10 ];        % TR 38.901 clause 7.6.3.3
elseif SC_lambda == 0
    rnd = [ rand( n_pos,3 ), randn( n_pos,1 ) ];
    return
else
    lambda = SC_lambda * ones( 1,4 );
end
rnd = [ qd_sos.rand( lambda(1:3), pos ).', qd_sos.randn( lambda(4), pos ).' ];
end

function [ PL, d_2d_in ] = o2i_mmmagic( f_GHz, low_loss_fraction, SC_lambda, max_indoor_dist, gen_per_tx, pos, pos_tx )
% O2I loss of H2020-ICT-671650-mmMAGIC/D2.2, Sec. 4.3

no_tx   = size( pos_tx,3 );
no_seg  = size( pos,2 );
no_freq = numel( f_GHz );
o_tx    = ones( no_tx,1 );
o_frq   = ones( 1,no_freq );

% Material penetration losses and external wall loss of TR 38.901 v14.1.0, Table 7.4.3-1
L_glass    =  2    + 0.2  * f_GHz;
L_IRglass  = 23    + 0.3  * f_GHz;
L_concrete =  5    + 4    * f_GHz;
PL_tw_low  = 5 - 10*log10( 0.3 * 10.^( -L_glass/10 )   + 0.7 * 10.^( -L_concrete/10 ) );
PL_tw_high = 5 - 10*log10( 0.7 * 10.^( -L_IRglass/10 ) + 0.3 * 10.^( -L_concrete/10 ) );

% Low-loss or high-loss building
if SC_lambda == 0
    randC = rand( 1, no_seg );
else
    randC = qd_sos.rand( SC_lambda , pos );
end
use_low_loss  = randC <= low_loss_fraction;
PL_tw = zeros( no_tx, no_seg, no_freq );
for f = 1 : no_freq
    PL_tw( :,use_low_loss,f )  = PL_tw_low(1,f);
    PL_tw( :,~use_low_loss,f ) = PL_tw_high(1,f);
end

% Indoor 2D distance from the CDF of the minimum of two uniform random variables
if SC_lambda == 0
    if gen_per_tx
        randC = rand( no_tx, no_seg );
    else
        randC = rand( 1, no_seg );
    end
elseif gen_per_tx
    randC = qd_sos.rand( SC_lambda*ones(1,no_tx) , pos );
else
    randC = qd_sos.rand( SC_lambda , pos );
end
val = 0:0.1:1;
cdf = [ 0, 0.193, 0.36, 0.514, 0.642, 0.749, 0.84, 0.911, 0.962, 0.991, 1 ];
d_2d_in = quadriga_lib.interp( cdf, [], val, randC );
if gen_per_tx
    d_2d_in = reshape( d_2d_in, size(randC) );
else
    d_2d_in = o_tx * d_2d_in;
end
d_2d_in = d_2d_in * max_indoor_dist;

% 2D distance and elevation angle
d_2d  = zeros( no_tx, no_seg );
theta = zeros( no_tx, no_seg );
for t = 1 : no_tx
    pos_tmp = pos - pos_tx(:,:,t);
    d_2d(t,:)  = sqrt(sum( abs( pos_tmp(1:2,:) ).^2,1 ));
    theta(t,:) = atand( pos_tmp(3,:) ./ d_2d(t,:) );
end
d_2d_in( d_2d_in > d_2d ) = d_2d( d_2d_in > d_2d );

% Indoor loss
if SC_lambda == 0
    randC = rand( 1, no_seg );
else
    randC = qd_sos.rand( SC_lambda , pos );
end
alpha = 0.5 + randC;
if SC_lambda == 0
    randC = rand( 1, no_seg );
else
    randC = qd_sos.rand( SC_lambda , pos );
end
d_room = 2 + 2*randC;
PL_in = alpha(o_tx,:) .* max( d_2d_in - d_room(o_tx,:) , 0 );
PL_in = PL_in( :,:,o_frq );

% Elevation-dependent loss
PL_el = 20*abs(theta/90);
PL_el = PL_el(:,:,o_frq);

% Deviation
sig_P = 4 + 0.1 * f_GHz;
if SC_lambda == 0
    if gen_per_tx
        randC = rand( no_tx, no_seg );
    else
        randC = o_tx * rand( 1, no_seg );
    end
elseif gen_per_tx
    randC = qd_sos.rand( SC_lambda*ones(1,no_tx) , pos );
else
    randC = o_tx * qd_sos.rand( SC_lambda , pos );
end
PL_sf = zeros( no_tx, no_seg, no_freq );
for f = 1 : no_freq
    PL_sf( :,:,f ) = randC * sig_P(1,f);
end

PL = PL_tw + PL_el + PL_in + PL_sf;
end
