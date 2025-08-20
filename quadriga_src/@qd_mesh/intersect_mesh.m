function [ islos, no_trans, fbs, sbs, iFBS, iSBS ] = intersect_mesh( h_mesh, orig, dest, obj_id )
%INTERSECT_MESH Calculates ray-mesh intersections
%
% Calling object:
%   Single object
%
% Description:
%   A ray is defined by its origin 3D-position and its destination 3d-position. Along the line from
%   origin to destination, there may or may not be an object that blocks that line. This functions
%   calculates for a set of N rays if the line-of-sight (LOS) is blocked. If it is blocked, the
%   coordinates at which the interaction happened are calculates as well. 
%
% Input:
%   orig
%   Ray origins in 3D Cartesian coordinates; Dimensions: (3xN) or (3x1)
%
%   dest
%   Ray destinations in 3D Cartesian coordinates; Dimensions: (3xN) or (3x1)
%
%   obj_id
%   A vector containing the object indices that should be considered. Default: All
%
%
% Output:
%   islos
%   Logical vector indicating if an intersection happened; Dimensions: (1xN)
%
%   no_trans
%   Number of intersections for each ray; uint64 vector; Dimensions: (1xN)
%
%   fbs
%   First interaction point; single precision; Dimensions: (3xN)
%
%   sbs
%   Second interaction point; single precision; Dimensions: (3xN)
%
%   iFBS
%   Index of the first mesh element that was hit by the ray; uint64; Dimensions: (1xN)
%
%   iFBS
%   Index of the second mesh element that was hit by the ray; uint64; Dimensions: (1xN)
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
    error('QuaDRiGa:qd_mesh:intersect_mesh','intersect not definded for object arrays.');
else
    h_mesh = h_mesh(1,1); % workaround for octave
end

if ~exist( 'orig','var' ) || size(orig,1) ~= 3
    error('QuaDRiGa:qd_mesh:intersect','"orig" is not given or has wrong format');
end

if ~exist( 'dest','var' ) || size(dest,1) ~= 3
    error('QuaDRiGa:qd_mesh:intersect','"dest" is not given or has wrong format');
end

use_single = false;
if isa(orig,'single')
    use_single = true;
    dest = single(dest);
end

if size(orig,2) == 1
    orig = repmat(orig, 1, size(dest,2));
end

use_object_id = false;
if ~exist('obj_id','var') || isempty( obj_id )
    obj_id = true(size(h_mesh.obj_index));
else
    use_object_id = true;
    obj_index = uint64( h_mesh.obj_index );
    ii = false( size( obj_index ) );
    for n = 1 : numel( obj_id )
        ii = ii | obj_index == uint64( obj_id(n) );
    end
    obj_id = ii;
end

% Read the vertices from the mesh
if use_object_id
    mesh = h_mesh.mesh(obj_id,:);
else
    mesh = h_mesh.mesh;
end

if use_single
    mesh = single( mesh );
else
    mesh = double( mesh );
end

if use_object_id || isempty( h_mesh.Psub_mesh_index )
    [ fbs, sbs, no_trans, iFBS, iSBS ] = quadriga_lib.ray_triangle_intersect( orig', dest', mesh );
else
    [ fbs, sbs, no_trans, iFBS, iSBS ] = quadriga_lib.ray_triangle_intersect( orig', dest', mesh, h_mesh.Psub_mesh_index - 1 );
end

islos = no_trans == uint64(0);

end
