function [ beamwidth_az, beamwidth_el, az_point_ang, el_point_ang ] =...
    calc_beamwidth( h_qd_arrayant, i_element, thres_dB )
%CALC_BEAMWIDTH Calculates the beam width for each antenna element in [deg]
%
% Calling object:
%   Single object
%
% Description:
%   This method calculates the beamwidth in azimuth and elevation direction as well as the pointing
%   angles of each element of the array antenna. Interpolation is used to achieve a higher
%   precision as provided by the sampling angle grid.
%
% Input:
%   i_element
%   A list of element indices. Default: 1 ... no_elements
%
%   thres_dB
%   The threshold in dB (Default: 3 dB, equivalent to FWHM)
%
% Output:
%   beamwidth_az
%   The azimuth beamwidth for each element in [deg]
%
%   beamwidth_el
%   The elevation beamwidth for each element in [deg]
%
%   az_point_ang
%   The azimuth pointing angle for the main beam in [deg]
%
%   el_point_ang
%   The elevation pointing angle for the main beam in [deg]
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

if ~exist('i_element','var') || isempty(i_element)
    i_element = 1:h_qd_arrayant.no_elements;
elseif ~(any(size(i_element) == 1) && isnumeric(i_element) ...
        && isreal(i_element) && all(mod(i_element, 1) == 0) && all(i_element > 0))
    error('??? "i_element" must be integer and > 0')
elseif any(i_element > h_qd_arrayant.no_elements)
    error('??? "i_element" exceeds "no_elements"')
end

if ~exist('thres_dB','var') || isempty(thres_dB)
    thres_dB = 3;
end

[ beamwidth_az, beamwidth_el, az_point_ang, el_point_ang ] = quadriga_lib.arrayant_calc_beamwidth( ...
    real(h_qd_arrayant.Fa), imag(h_qd_arrayant.Fa), real(h_qd_arrayant.Fb), imag(h_qd_arrayant.Fb), ...
    h_qd_arrayant.azimuth_grid, h_qd_arrayant.elevation_grid, i_element, thres_dB );

end
