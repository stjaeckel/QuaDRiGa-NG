function [ lsp, abs_toa, gr_eps ] = lsf_eval( h_builder, tx_pos, rx_pos, dual_mobility, i_ref )
%LSF_EVAL Large-scale parameters at a set of links
%
%   Evaluates the large-scale generators of the builder ('lsf_generator') at the links with
%   'quadriga_lib.lsf_draw' and generates the large-scale parameters with
%   'quadriga_lib.lsp_generate'. With 'lsf_model = 38.901' and in 3GPP baseline mode, the shadow-
%   fading standard deviation is given by the path-loss model. If 'i_ref' is given, the
%   standard deviation and the environment height of link 'i_ref' apply to all links (e.g. the
%   segment start of a track).
%
% Input:
%   tx_pos          Tx positions in [m]; [3, 1] or [3, n_link]
%   rx_pos          Rx positions in [m]; [3, n_link]
%   dual_mobility   true: the Tx is a second mobile device
%   i_ref           Optional index of the reference link
%
% Output:
%   lsp             DS [s], KF, SF (linear), ASD, ASA, ESD, ESA [deg], XPR (linear); [n_link, n_freq, 8]
%   abs_toa         Absolute TOA offset in [s]; [n_link, 1]
%   gr_eps          Relative permittivity of the ground; [n_link, n_freq]
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

fc = h_builder.simpar(1,1).center_frequency;
n_link = size( rx_pos,2 );
if n_link == 0
    lsp = zeros( 0, numel(fc), 8 );
    abs_toa = zeros( 0, 1 );
    gr_eps = zeros( 0, numel(fc) );
    return
end

gen = h_builder.lsf_generator;
rnd = quadriga_lib.lsf_draw( gen.lambda, gen.phase, tx_pos.', rx_pos.', [], logical( dual_mobility ), ...
    gen.preset, gen.horizontal );

use_tr = ~strcmp( h_builder.simpar(1,1).lsf_model, 'legacy' );
sf_sigma = [];
if use_tr || h_builder.simpar(1,1).use_3GPP_baseline
    rnd_he = rnd(:,11);
    if exist( 'i_ref','var' )
        rnd_he(:) = rnd_he( i_ref );
    end
    [ ~, sf_sigma ] = quadriga_lib.path_loss( h_builder.plpar, h_builder.scenpar, fc, tx_pos.', rx_pos.', [], rnd_he );
    if exist( 'i_ref','var' )
        sf_sigma = ones( n_link,1 ) * sf_sigma( i_ref,: );
    end
end

[ lsp, abs_toa, gr_re, gr_im ] = quadriga_lib.lsp_generate( h_builder.scenpar, fc, tx_pos.', rx_pos.', ...
    rnd, sf_sigma, use_tr, use_tr );
gr_eps = complex( gr_re, gr_im );

end
