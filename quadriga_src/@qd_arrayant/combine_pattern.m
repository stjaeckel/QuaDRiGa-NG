function combine_pattern( h_qd_arrayant, center_frequency )
%COMBINE_PATTERN Calculates a virtual pattern of the given array
%
% Calling object:
%   Single object
%
% Description:
%   When the inputs of an array antenna are coupled (i.e. fed with the same signal), then it is
%   possible to combine the elements of the array. This function calculates the virtual pattern by
%   using the QuaDRiGa simulator. Individual coupling weights can be set in the "coupling" property
%   of the qd_arrayant object. Phase offsets of the individual antenna elements due to their
%   positions in the array ("element_position" property of the calling qd_arrayant object) are
%   calculated for the phase center of the array.
%
% Input:
%   center_frequency
%   The center frequency in [Hz]. If this input variable is not given, it is assumed that the
%   element spacings in the "element_position" property of the calling arrayant object are given in
%   multiples of the carrier wavelength.
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
    error('QuaDRiGa:qd_arrayant:combine_pattern','combine_pattern not definded for object arrays.');
else
    h_qd_arrayant = h_qd_arrayant(1,1); % workaround for octave
end

% We assume that spacings can be given in units "lambda".
% The base-frequency is therefor ~ 300 MHz , so the wavelength is 1 m
if ~exist('center_frequency','var')
    center_frequency = h_qd_arrayant.center_frequency;
end

% Call quadriga-lib library function
pat = quadriga_lib.arrayant_combine_pattern([], center_frequency, [], [], ...
        real(h_qd_arrayant.PFa), imag(h_qd_arrayant.PFa), ...
        real(h_qd_arrayant.PFb), imag(h_qd_arrayant.PFb), ...
        h_qd_arrayant.azimuth_grid, h_qd_arrayant.elevation_grid, ...
        h_qd_arrayant.element_position, real(h_qd_arrayant.coupling), imag(h_qd_arrayant.coupling), ...
        h_qd_arrayant.center_frequency, h_qd_arrayant.name);

% Write the output pattern
no_elements = size( pat.e_theta_re, 3 );
h_qd_arrayant.PFa = complex( pat.e_theta_re, pat.e_theta_im );
h_qd_arrayant.PFb = complex( pat.e_phi_re, pat.e_phi_im );
h_qd_arrayant.Pelement_position = zeros(3, no_elements);
h_qd_arrayant.Pcoupling = eye( no_elements );
h_qd_arrayant.Pphase_diff = [];
h_qd_arrayant.Pno_elements = no_elements;

end
