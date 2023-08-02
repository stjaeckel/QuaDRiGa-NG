function rotate_pattern( h_qd_arrayant, deg, rotaxis, i_element, usage )
%ROTATE_PATTERN Rotates antenna patterns
%
% Calling object:
%   Single object
%
% Description:
%   Pattern rotation provides the option to assemble array antennas out of single elements. By
%   setting the 'element_position' property of an array object, elements can be placed at different
%   coordinates. In order to freely design arbitrary array configurations, however, elements often
%   need to be rotated (e.g. to assemble a +/- 45° crosspolarized array out of single dipoles).
%   This functionality is provided here.
%
% Input:
%   deg
%   The rotation angle in [degrees] ranging from -180° to 180°
%
%   rotaxis
%   The rotation axis specified by the string 'x','y', 'z', or 'xyz'. In case of 'xyz', three
%   rotations are performed: one around the 'x' axis, followed by one around the 'y' axis and
%   around the 'z' axis. In this case, the input 'deg' must be a 3-element vector containing the 3
%   rotation angles.
%
%   i_element
%   The element indices for which the rotation is done. If no element index is given, the rotation
%   is done for all elements in the array.
%
%   usage
%   The optional parameter 'usage' can limit the rotation procedure either to the pattern or
%   polarization. Possible values are:
%      * 0: Rotate both (pattern+polarization) - default
%      * 1: Rotate only pattern
%      * 2: Rotate only polarization
%      * 3: Same as (0), but without grid interpolation
%
% QuaDRiGa Copyright (C) 2011-2023
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
   error('QuaDRiGa:qd_arrayant:rotate_pattern','calc_gain not definded for object arrays.');
else
    h_qd_arrayant = h_qd_arrayant(1,1); % workaround for octave
end

% Parse input arguments
if exist('rotaxis','var')
    if ~( ischar(rotaxis) && ...
            (strcmp(rotaxis,'x') || strcmp(rotaxis,'y') || strcmp(rotaxis,'z') || strcmp(rotaxis,'xyz') ) )
        error('QuaDRiGa:qd_arrayant:rotate_pattern','??? "rotaxis" can only be x,y or z.')
    end
else
    rotaxis = 'y';
end

if exist('deg','var')
    if ~( isnumeric(deg) && isreal(deg) )
        error('QuaDRiGa:qd_arrayant:rotate_pattern','??? "deg" must be real-valued')
    end
else
    deg = 0;
end

rotate_all = false;
if exist('i_element','var') && ~isempty( i_element )
    if ~( size(i_element,1) == 1 && isnumeric(i_element) ...
            &&  all( mod(i_element,1)==0 ) && min(i_element) > 0 && max(i_element)<=h_qd_arrayant.no_elements )
        error('QuaDRiGa:qd_arrayant:rotate_pattern',...
            '??? "i_element" must be integer > 0 and can not exceed the number of elements')
    end
    i_element = unique( i_element );
    if numel(i_element) == h_qd_arrayant.no_elements && all( ((1:h_qd_arrayant.no_elements) - (i_element)) == 0)
        rotate_all = true;
    end
else
    i_element = 1:h_qd_arrayant.no_elements;
    rotate_all = true;
end

if exist('usage','var')
    if ~( all(size(usage) == [1 1]) && isnumeric(usage) ...
            && any(usage == [0,1,2,3]) )
        error('QuaDRiGa:qd_arrayant:rotate_pattern','??? "usage" must be 0,1,2 or 3')
    end
else
    usage = 0;
end

% Determine if we use single or double precision
if isa( h_qd_arrayant.Fa, 'single' )
    single( h_qd_arrayant );
else
    double( h_qd_arrayant );
end

zrot = 0;
yrot = 0;
xrot = 0;
switch rotaxis
    case 'x'
        xrot = deg;
    case 'y'
        yrot = deg;
    case 'z'
        zrot = deg;
    case 'xyz'
        xrot = deg(1);
        yrot = deg(2);
        zrot = deg(3);
end

% Call quadriga-lib library function
V = h_qd_arrayant.PFa(:,:,i_element);
H = h_qd_arrayant.PFb(:,:,i_element);

[ Vr, Vi, Hr, Hi, azimuth_grid, elevation_grid, element_pos] ...
    = quadriga_lib.arrayant_rotate_pattern( real(V), imag(V), real(H), imag(H), ...
    h_qd_arrayant.azimuth_grid, h_qd_arrayant.elevation_grid, h_qd_arrayant.Pelement_position(:,i_element), ...
    xrot, yrot, zrot, usage);

clear V;
clear H;

% It is not possible to rotate a subset of elements if the angular grid needs to be interpolated.
if ~rotate_all && (h_qd_arrayant.no_az ~= numel(azimuth_grid) || h_qd_arrayant.no_el ~= numel(elevation_grid) )
    error('QuaDRiGa:qd_arrayant:rotate_pattern',...
        'Pattern requires interpolation of the angle grid. You cannot select individual elements in this case.')
end

if h_qd_arrayant.no_az ~= numel(azimuth_grid) || ...
        h_qd_arrayant.no_el ~= numel(elevation_grid) || ...
        any(h_qd_arrayant.azimuth_grid - azimuth_grid' ~= 0) || ...
        any(h_qd_arrayant.elevation_grid - elevation_grid' ~= 0) 
    h_qd_arrayant.set_grid(azimuth_grid, elevation_grid, 0);
end

h_qd_arrayant.PFa(:,:,i_element) = complex( Vr, Vi );
h_qd_arrayant.PFb(:,:,i_element) = complex( Hr, Hi );
h_qd_arrayant.Pelement_position(:,i_element) = element_pos;
h_qd_arrayant.Pphase_diff = [];

end
