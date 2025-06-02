function [ directivity_dBi, gain_dBi ] = calc_gain( h_qd_arrayant, i_element )
%CALC_GAIN Calculates the gain in dBi of the array antenna
%
% Calling object:
%   Single object
%
% Input:
%   i_element
%   A list of element indices.
%
% Output:
%   directivity_dBi
%   Normalized gain of the antenna in dBi.
%
%   gain_dBi
%   Maximum gain of the antenna in dBi (gain = directivity - losses)
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

if numel( h_qd_arrayant ) > 1
    error('QuaDRiGa:qd_arrayant:calc_gain','calc_gain not definded for object arrays.');
else
    h_qd_arrayant = h_qd_arrayant(1,1); % workaround for octave
end

if nargin == 1
    i_element = 1:h_qd_arrayant.no_elements;
end

if ~(any(size(i_element) == 1) && isnumeric(i_element) ...
        && isreal(i_element) && all(mod(i_element, 1) == 0) && all(i_element > 0))
    error('??? "i_element" must be integer and > 0')
elseif any(i_element > h_qd_arrayant.no_elements)
    error('??? "i_element" exceeds "no_elements"')
end

if isa( h_qd_arrayant.PFa, 'single' )
    single( h_qd_arrayant );
else % double
    double( h_qd_arrayant );
end

directivity_dBi = zeros(numel(i_element),1);
gain_dBi = zeros(numel(i_element),1);

az = h_qd_arrayant.azimuth_grid;
el = h_qd_arrayant.elevation_grid;

for n = 1 : numel(i_element)
    
    % Read the qd_arrayant elements
    Fa = h_qd_arrayant.PFa(:, :, i_element(n));
    Fb = h_qd_arrayant.PFb(:, :, i_element(n));
    
    directivity_dBi(n) = quadriga_lib.arrayant_calc_directivity( real(Fa), imag(Fa), real(Fb), imag(Fb), az, el );
    
    if nargout > 1
        P = abs(Fa).^2 + abs(Fb).^2;
        P_max = max( P(:) );
        gain_dBi(n) = 10*log10(P_max);
    end
end

end

