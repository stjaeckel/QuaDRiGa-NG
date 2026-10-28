function [ lsp_mu, lsp_sigma ] = lsp_distribution( h_builder, tx_pos, rx_pos )
%LSP_DISTRIBUTION Mean and standard deviation of the large-scale parameters at a set of links
%
%   Returns the distribution of the large-scale parameters in the log domain (lg for DS, ASD, ASA,
%   ESD, ESA; dB for KF, SF, XPR) from 'quadriga_lib.lsp_generate', including the distance,
%   height, elevation and frequency terms of the scenario.
%
% Input:
%   tx_pos          Tx positions in [m]; [3, 1] or [3, n_link]
%   rx_pos          Rx positions in [m]; [3, n_link]
%
% Output:
%   lsp_mu          Mean; [n_link, n_freq, 8]
%   lsp_sigma       Standard deviation; [n_link, n_freq, 8]
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
    lsp_mu = zeros( 0, numel(fc), 8 );
    lsp_sigma = lsp_mu;
    return
end
use_tr = ~strcmp( h_builder.simpar(1,1).lsf_model, 'legacy' );
[ ~, ~, ~, ~, ~, lsp_mu, lsp_sigma ] = quadriga_lib.lsp_generate( h_builder.scenpar, fc, ...
    tx_pos.', rx_pos.', zeros( n_link, 11 ), [], use_tr, use_tr );

end
