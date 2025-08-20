function generate_index( h_mesh, target_size, vec_size )
%GENERATE_INDEX Generates the sub-mesh index
%
% Calling object:
%   Single object
%
% Description:
%   This function generates the sub-mesh index. By splitting a triangle mesh into sub-meshes, a 
%   significant increase in computing performance can be obtained. If a ray does not interact with 
%   face within a sub-mesh, the entire sub-mesh can be skipped. However, this requires that the 
%   order of the faces in the mesh is changed.
%
% Input:
%   target_size
%   Target value for the sub-mesh size, default = max( 1024, sqrt(no_face) )
%
%   vec_size
%   Vector size for SIMD processing (e.g. 8 for AVX2, 32 for CUDA), default = autodetect
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
    error('QuaDRiGa:qd_mesh:diff_trans','generate_index is not definded for object arrays.');
else
    h_mesh = h_mesh(1,1); % workaround for octave
end

if ~exist( 'target_size','var' ) || isempty(target_size)
    target_size = max( ceil( sqrt( h_mesh.no_face )*1.5 ), 1024 );
end

if ~exist( 'vec_size','var' ) || isempty(vec_size)
    vec_size = 1;
    s = qd_simulation_parameters;
    if strcmp( s.quadriga_lib_version(end-3:end), 'AVX2' )
        vec_size = 8;
    end
end

if h_mesh.no_face <= target_size
    h_mesh.Psub_mesh_index = [];

else
    [ msh, si, mi ] = quadriga_lib.triangle_mesh_segmentation( h_mesh.mesh, target_size, vec_size );
    % Note: si is 0-based and mi is 1-based

    % Add Padding
    xi = [si ;uint64( numel(mi) )];
    for n = 1 : numel(si)

        % Face indes of the last element in the current sub-mesh
        last_sub_face = mi( xi(n+1) );
        if last_sub_face == 0
            h_mesh.Pvert = [ h_mesh.Pvert, msh( xi(n+1) , 1:3 )' ];
            h_mesh.Pface = [ h_mesh.Pface, uint64( h_mesh.no_vert([1;1;1]) ) ];
            h_mesh.Pobj_index = [ h_mesh.Pobj_index, 1 ];
            h_mesh.Pmtl_index = [ h_mesh.Pmtl_index, 1 ];
            
            tmp = uint64( h_mesh.no_face );
            for m = xi(n+1) - vec_size + 2 : xi(n+1)
                if mi(m) == 0
                    mi(m) = tmp;
                end
            end
        end
    end

    % Reorder mesh
    h_mesh.Pface = h_mesh.Pface( :,mi );
    h_mesh.Pobj_index = h_mesh.Pobj_index( mi );
    h_mesh.Pmtl_index = h_mesh.Pmtl_index( mi );
    h_mesh.Psub_mesh_index = si' + 1;
end
end

