function gain = diff_trans( h_mesh, orig, dest, center_frequency, lod, obj_id, verbose )
%DIFF_TRANS Calculates the gain including diffraction and transmission losses
%
% Calling object:
%   Single object
%
% Description:
%   This function implements a diffraction and transmission loss model for the direct path from
%   transmitter to receiver. Free-Space path-loss is not included in the calculated gain. Upon
%   availability, this function uses GPU-acceleration to improve the computing performance.
%
% Input:
%   orig
%   Ray origins in 3D Cartesian coordinates; Dimensions: (3xN) or (3x1)
%
%   dest
%   Ray destinations in 3D Cartesian coordinates; Dimensions: (3xN) or (3x1)
%
%   center_frequency
%   Center frequency in [Hz]
%
%   lod
%   Level of detail, scalar value 0-6

%   obj_id
%   A vector containing the object indices that should be included. Default: All
%
%   verbose
%   Enables (1, default) or disables (0) progress report.
%
% Output:
%   gain
%   Gain caused by diffraction and transmission effects; linear scale; Dimensions: (1xN)
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

if numel( h_mesh ) > 1
    error('QuaDRiGa:qd_mesh:diff_trans','intersect not definded for object arrays.');
else
    h_mesh = h_mesh(1,1); % workaround for octave
end

if ~exist( 'orig','var' ) || size(orig,1) ~= 3
    error('QuaDRiGa:qd_mesh:diff_trans','"orig" is not given or has wrong format');
end

if ~exist( 'dest','var' ) || size(dest,1) ~= 3
    error('QuaDRiGa:qd_mesh:diff_trans','"dest" is not given or has wrong format');
end

use_single = false;
if isa(orig,'single')
    use_single = true;
    dest = single(dest);
end

if size(orig,2) == 1
    orig = repmat(orig, 1, size(dest,2));
end

if ~exist( 'center_frequency','var' ) || size(center_frequency,1) ~= 1
    error('QuaDRiGa:qd_mesh:diff_trans','"center_frequency" is not given or has wrong format');
end

if ~exist( 'lod','var' ) || isempty(lod)
   lod = 2;
end

if ~exist('obj_id','var') || isempty( obj_id )
    obj_id = true(size(h_mesh.obj_index));
    obj_id_in = [];
else
    obj_id_in = obj_id;
    obj_index = uint64( h_mesh.obj_index );
    ii = false( size( obj_index ) );
    for n = 1 : numel( obj_id )
        ii = ii | obj_index == uint64( obj_id(n) );
    end
    obj_id = ii;
end

if ~exist( 'verbose','var' ) || isempty( verbose )
    verbose = 1;
end

% Read the vertices from the mesh
if use_single
    mesh = single( h_mesh.mesh(obj_id,:) );
    mtl_prop = single( h_mesh.mtl_prop )';
else
    mesh = double( h_mesh.mesh(obj_id,:) );
    mtl_prop = double( h_mesh.mtl_prop )';
end
mtl_prop = mtl_prop( h_mesh.mtl_index(obj_id), : );

if isempty( h_mesh.Psub_mesh_index )
    gain = quadriga_lib.calc_diffraction_gain( orig', dest', mesh,  mtl_prop, center_frequency, lod, verbose );
else
    gain = quadriga_lib.calc_diffraction_gain( orig', dest', mesh,  mtl_prop, center_frequency, lod, verbose, h_mesh.Psub_mesh_index - 1 );
end

end
