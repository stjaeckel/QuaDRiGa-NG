function [vertices, faces, faceColors] = visualize_pattern_geometry( h_qd_arrayant, ...
    i_element, position, orientation, scale, range_db, no_subdiv, cmap, edge_color, create_plot )
%VISUALIZE_PATTERN_GEOMETRY Plots the power pattern of an antenna
%
% Calling object:
%   Single object
%
% Description:
%   This function creates a plot of the antenna beam pattern. It can be used to visualize the
%   location and orientation of the BS antennas for a given simulation. The plot is 
%
% Input:
%   i_element
%   The element index for which the plot is created. If no element index are given, a plot is
%   created for the first element in the array. 
%
%   position
%   This position is given in global Cartesian coordinates (x,y, and z-component) in units of [m].
%   If no position is given, the plot is created at the origin of the cordinate system.
%
%   orientation
%   This 3-element vector describes the orientation of the antenna. The reference system for
%   aircraft principal axes is used. The first value describes the "roll angle", i.e. the rotation
%   around an axis drawn through the body of the vehicle from tail to nose in the normal direction
%   of movement. Positive rotation is clockwise (seen from the pilot/drivers perspective). The
%   second value describes the "pitch angle", i.e. the vertical (tilt) angle relative to the
%   horizontal plane; positive rotation is up. The third value describes the bearing or "yaw angle",
%   in mathematic sense. East corresponds to 0°, and the angles increase counter-clockwise, so north
%   is 90° degrees, south is -90°, and west is 180°. All values are given in [rad]. Note that 
%   by default, QuaDRiGa antennas face east (roll = 0, pitch = 0, yaw = 0).
%
%   scale
%   Scaling of the plot in units of [m]. Default: 1 m
%
%   range_db
%   The power range in dB for which the pattern is plotted, relative to tis maximum directivity.
%   Default: 20 dB
%
%   no_subdiv
%   Number of subdivisions in the generating Icosphere. Integer value. Default = 6
%
%   cmap
%   Colormat used for the face colors. Can be either a string or an array with 3 columns conraining
%   the RGB values in range 0 to 1. Default: 'jet'
%
%   edge_color
%   Edge color, specified as a string of GB value. The default edge color is black with a value of
%   [0 0 0]. 
%
%   create_plot
%   If set to 1, the output is plotted to the last opened figure. Default = 1
%
%
% Output:
%   vertices
%   Polygon vertices for the pattern plot
%
%   faces
%   Face definitions for the polygon plot
%
%   faceColors
%   Face colors as RGB values ranging from 0-1
%   
%   
%
% QuaDRiGa Copyright (C) 2011-2024
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
    error('QuaDRiGa:qd_arrayant:visualize_pattern_geometry',...
        '??? visualize_pattern_geometry is not definded for object arrays.');
else
    h_qd_arrayant = h_qd_arrayant(1,1); % workaround for octave
end

if ~exist('i_element','var') || isempty(i_element)
    i_element = 1;
elseif ~( size(i_element,1) == 1 && isnumeric(i_element) && all(size(i_element) == [1 1]) ...
        &&  all( mod(i_element,1)==0 ) && min(i_element) > 0 && max(i_element) <= h_qd_arrayant.no_elements )
    error('QuaDRiGa:qd_arrayant:visualize_pattern_geometry',...
        '??? "i_element" must be scalar, integer > 0 and can not exceed the number of antenna elements')
end

if ~exist('position','var') || isempty(position)
    position = [0;0;0];
elseif ~( size(position,1) == 3 && isnumeric(position) && all(size(position) == [3 1]) )
    error('QuaDRiGa:qd_arrayant:visualize_pattern_geometry',...
        '??? "orientation" must be a 3 element vector')
end

if ~exist('orientation','var') || isempty(orientation)
    orientation = [0;0;0];
elseif ~( size(orientation,1) == 3 && isnumeric(orientation) && all(size(orientation) == [3 1]) )
    error('QuaDRiGa:qd_arrayant:visualize_pattern_geometry',...
        '??? "orientation" must be a 3 element vector')
end

if ~exist('scale','var') || isempty(scale)
    scale = 1;
elseif ~( size(scale,1) == 1 && isnumeric(scale) && all(size(scale) == [1 1]) &&  scale > 0  )
    error('QuaDRiGa:qd_arrayant:visualize_pattern_geometry',...
        '??? "scale" must be scalar, integer > 0 and can not exceed the number of antenna elements')
end

if ~exist('range_db','var') || isempty(range_db)
    range_db = 20;
elseif ~( size(range_db,1) == 1 && isnumeric(range_db) && all(size(range_db) == [1 1]) &&  range_db > 0  )
    error('QuaDRiGa:qd_arrayant:visualize_pattern_geometry',...
        '??? "scale" must be scalar, integer > 0 and can not exceed the number of antenna elements')
end
min_value = -range_db;

if ~exist('no_subdiv','var') || isempty(no_subdiv)
    no_subdiv = 6;
elseif ~( size(no_subdiv,1) == 1 && isnumeric(no_subdiv) && all(size(no_subdiv) == [1 1]) &&  no_subdiv > 0  )
    error('QuaDRiGa:qd_arrayant:visualize_pattern_geometry',...
        '??? "no_subdiv" must be scalar, integer > 0')
end

if ~exist('cmap','var') || isempty(cmap)
    cmap_values = jet(256);
elseif ischar(cmap)
    cmap_values = feval(cmap, 256);
elseif isnumeric(cmap) && size(cmap,2) == 3 && min(cmap(:)) >= 0 && max(cmap(:)) <= 1
    cmap_values = cmap;
else
    error('QuaDRiGa:qd_arrayant:visualize_pattern_geometry',...
        '??? "cmap" must be a character array or a [256,3] matrix having values between 0 and 1.')
end

if ~exist('edge_color','var') || isempty(edge_color)
    edge_color = 'k';
end

if ~exist('create_plot','var') || isempty(create_plot)
    create_plot = 1;
end

% Generate Icosphere
[ center, ~, vert ] = quadriga_lib.icosphere( no_subdiv, 1 );
no_faces = size(center,1);

vertices = [ center + vert(:,1:3); center + vert(:,4:6); center + vert(:,7:9) ];
faces = reshape(1:3*no_faces,[],3);

% Get power values
a = h_qd_arrayant.sub_array( i_element );

% Get angles for pattern interpoaltion
[azimuth, elevation] = quadriga_lib.cart2geo(vertices');

% Interpoalte pattersn
[V_re, V_im, H_re, H_im] = ...
    quadriga_lib.arrayant_interpolate( real(a.Fa), imag(a.Fa), real(a.Fb), imag(a.Fb), ...
    a.azimuth_grid, a.elevation_grid, azimuth', elevation', 1, orientation, a.element_position );

% Calc power scaling
P = V_re.^2 + V_im.^2 + H_re.^2 + H_im.^2;
P = P./max(P);
P = 10*log10(P);
P(P < min_value) = min_value;
P = (P - min_value) ./ (0 - min_value) .* scale;

% Scale verices
vertices = vertices .* P';
vertices = vertices + position' + a.element_position';

% Get colormap
[azimuth, elevation] = quadriga_lib.cart2geo(center');
[V_re, V_im, H_re, H_im] = ...
    quadriga_lib.arrayant_interpolate( real(a.Fa), imag(a.Fa), real(a.Fb), imag(a.Fb), ...
    a.azimuth_grid, a.elevation_grid, azimuth', elevation', 1, orientation, a.element_position );
P = V_re.^2 + V_im.^2 + H_re.^2 + H_im.^2;
P = P./max(P);
P = 10*log10(P);
P(P < min_value) = min_value;

no_cmap_values = size(cmap_values,1) - 1;

P = (P + range_db) * (no_cmap_values / range_db) + 1;
P = round(P);

faceColors = cmap_values( P,: );

% Create the figure and use the patch function to render the mesh with different colors for each face
if create_plot
    patch('Vertices', vertices, 'Faces', faces, 'FaceVertexCData', faceColors, 'FaceColor', 'flat', 'EdgeColor', edge_color);
end

end
