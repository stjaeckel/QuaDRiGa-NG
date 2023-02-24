function [ V, H, dist, azimuth, elevation ] = interpolate( h_qd_arrayant, azimuth, elevation, i_element,...
    orientation, threshold, use_gpu )
%INTERPOLATE Interpolates the array antenna field patterns
%
% Calling object:
%   Single object
%
% Description:
%   This function interpolates the polarimetric antenna field patterns for a given set of azimuth
%   and elevation angles. Interpolation of the antenna field patterns is very computing intensive.
%   It must be performed several thousands of times during a simulation run. Therefore, the
%   function implements two interpolation methods: 2D linear interpolation and spheric
%   interpolation. Linear interpolation is used for real-valued field patterns and for complex-
%   valued field pattern that have a small phase-variation between neighboring elements. However,
%   linear interpolation performs very poorly when there is a large phase difference between two
%   neighboring samples of the pattern. In this case, spheric interpolation is used.
%
% Input:
%   azimuth
%   A vector of azimuth angles in [rad]. The default dimensions are: [ 1, no_ang ]. It is possible
%   to provide a different angle for each element of the array antenna. In this case, the
%   dimensions are [ no_elements, no_ang ], where 'no_elements' corresponds to the number of
%   entries in 'i_element' or the number of elements in the array antenna if 'i_element' is not
%   given.
%
%   elevation
%   A vector of elevation angles in [rad]. The dimensions are: [ 1, no_ang ] or in case of per-
%   element angles [ no_elements, no_ang ].
%
%   i_element
%   The element indices for which the interpolation is done. If no element index is given, the
%   interpolation is done for all elements in the array. Dimensions: [ 1, no_elements ]
%
%   orientation
%   This (optional) 3-element vector describes the orientation of the array antenna. The The first
%   value describes the ”bank angle”, the second value describes the ”tilt angle”, (positive values
%   point upwards), the third value describes the bearing or ”heading angle”, in mathematic sense.
%   East corresponds to 0, and the angles increase counter-clockwise, so north is 90 degrees, south
%   is -90 degree, and west is 180 degree. All values are given in [rad]. By default, the
%   orientation is [0;0;0], i.e. the broadside of the antenna points at the horizon towards the
%   East.
%
%   threshold
%   The maximum phase difference in [deg] between two neighboring entries in the antenna field
%   pattern for which linear interpolation is allowed. Linear interpolation is much faster, but
%   also inaccurate for large phase offsets. By default, a 0 degree threshold is used, i.e.,
%   spheric polarization is used for all complex-valued patterns.
%
%   use_gpu
%   Enables or disables (default) GPU acceleration. This requires a compatible GPU and the
%   "Parallel Computing Toolbox" for MATLAB. In Octave, GPU acceleration is available through the
%   "ocl"-package (https://octave.sourceforge.io/ocl).
%
% Output:
%   V
%   The interpolated vertical field pattern (E-θ-component). Dimensions: [ no_elements, no_ang ]
%
%   H
%   The interpolated horizontal field pattern (E-ϕ-component).  Dimensions: [ no_elements, no_ang ]
%
%   dist
%   The effective distances between the antenna elements when seen from the direction of the
%   incident path. The distance is calculated by an projection of the array positions on the normal
%   plane of the incident path. This is needed for the planar wave approximation.
%   Dimensions: [ no_elements, no_ang ]
%
%   azimuth
%   The azimuth angles in [rad] for the local antenna coordinate system, i.e., after applying the
%   'orientation'. If no orientation vector is given, these angles are identical to the input
%   azimuth angles.
%
%   elevation
%   The elevation angles in [rad] for the local antenna coordinate system, i.e., after applying the
%   'orientation'. If no orientation vector is given, these angles are identical to the input
%   elevation angles.
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

if numel( h_qd_arrayant ) > 1
    error('QuaDRiGa:qd_arrayant:interpolate','interpolate not definded for object arrays.');
else
    h_qd_arrayant = h_qd_arrayant(1,1); % workaround for octave
end

if ~exist( 'i_element','var' ) || isempty( i_element )
    i_element = 1 : h_qd_arrayant.Pno_elements;
end

i_element = uint32( i_element );

if ~exist( 'orientation','var' ) || isempty( orientation )
    orientation = [];
elseif all( abs(orientation) < 1e-6 )
    orientation = [];
end

if ~exist( 'threshold','var' ) || isempty( threshold )
    threshold = 0; % Always use circular interpolarion
end

if ~exist( 'use_gpu','var' ) || isempty( use_gpu )
    use_gpu = 0;
elseif logical( use_gpu ) && ~qd_simulation_parameters.has_gpu
    use_gpu = 0;
end

% Adjust single/double if needed
use_double = isa( h_qd_arrayant.PFa, 'double' );
if ~isa( h_qd_arrayant.PFb, 'double' ) && use_double
    h_qd_arrayant.PFb = double( h_qd_arrayant.PFb );
end
if ~isa( h_qd_arrayant.azimuth_grid, 'double' ) && use_double
    h_qd_arrayant.azimuth_grid = double( h_qd_arrayant.azimuth_grid );
end
if ~isa( h_qd_arrayant.elevation_grid, 'double' ) && use_double
    h_qd_arrayant.elevation_grid = double( h_qd_arrayant.elevation_grid );
end
if ~isa( h_qd_arrayant.element_position, 'double' ) && use_double
    element_position = double( h_qd_arrayant.element_position(:,i_element) );
else
    element_position = h_qd_arrayant.element_position(:,i_element);
end

% Match azimuth angle dimensions
no_element = numel( i_element );
no_ang = size( azimuth,2 );
per_element_interpol = size(azimuth,1) == no_element;
if ~per_element_interpol
    azimuth = azimuth(:)';
    elevation = elevation(:)';
end

% Adjust input types
if ~isa( azimuth, 'double' ) && use_double
    azimuth = double( azimuth );
end
if ~isa( elevation, 'double' ) && use_double
    elevation = double( elevation );
end
if ~isa( orientation, 'double' ) && use_double
    orientation = double( orientation );
end

if nargout < 4
    [Vr, Vi, Hr, Hi, dist] = quadriga_lib.arrayant_interpolate( real(h_qd_arrayant.PFa), imag(h_qd_arrayant.PFa), ...
        real(h_qd_arrayant.PFb), imag(h_qd_arrayant.PFb), ...
        h_qd_arrayant.azimuth_grid, h_qd_arrayant.elevation_grid, azimuth, elevation, ...
        i_element, orientation, element_position );
else
    [Vr, Vi, Hr, Hi, dist, azimuth, elevation ] = quadriga_lib.arrayant_interpolate( real(h_qd_arrayant.PFa), imag(h_qd_arrayant.PFa), ...
        real(h_qd_arrayant.PFb), imag(h_qd_arrayant.PFb), ...
        h_qd_arrayant.azimuth_grid, h_qd_arrayant.elevation_grid, azimuth, elevation, ...
        i_element, orientation, element_position );
end
V = complex(Vr,Vi);
H = complex(Hr,Hi);

if nargout > 3 % Return transformed azimuth and elevation angles
    azimuth = reshape( azimuth, [], no_ang );
    elevation = reshape( elevation, [], no_ang );
end

end
