function read_obj( h_mesh, fname )
%READ_OBJ Reads mesh data from Wavefront .obj file format
%
% Calling object:
%   Single object
%
% Description:
%   The OBJ file format is a simple data-format that represents 3D geometry - namely, the position
%   of each vertex and the faces that make each polygon defined as a list of vertices. This method
%   parses the OBJ file and reads the relevant mesh data into the calling 'qd_mesh' object. The
%   linked material library file name (mtllib) is read from the OBJ file and parsed separately. The
%   following conversions are made:
%
%   * Vertices (v) are stored as 'qd_mesh.vert'
%   * Face ids (f) are stored as 'qd_mesh.face'; texture coordinates and vertex normal are ignored
%   * Object names (o) are stored as 'qd_mesh.obj_name'
%   * The corresponding face ids belonging to this object as are stored as 'qd_mesh.obj_index'
%   * Materials (usemtl) are allocated to all faces belonging to an object
%   * The diffuse color (Kd) is read from the MTL file an stored as 'qd_mesh.mtl_color'
%   * The refraction index (Ni) is read from the MTL file and its squared value is used for the
%     relative permittivity. Conductivity is set to 0 and relative permeability is set to 1.
%   * Material thickness is set to 0.1
%
% Input:
%   fname
%   Path to the OBJ File
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

if numel( h_mesh ) > 1
    error('QuaDRiGa:qd_mesh:read_obj','intersect not definded for object arrays.');
else
    h_mesh = h_mesh(1,1); % workaround for octave
end

if ~exist( 'fname','var' ) || isempty( fname )
    error('QuaDRiGa:qd_mesh:read_obj','Filename is not given.');
end

[ ~, vert_list, face_ind, obj_ind, obj_name, mtl_ind, mtl_name, bsdf, csv_ind, ~, csv_prop ] = quadriga_lib.obj_file_read( fname );

% Load material properties, the diffuse color (Kd) is stored in the first 3 columns of the BSDF
no_existing_mtl = h_mesh.no_mtl;
no_new_mtl = numel(mtl_name);
mtl_color = ones(3,no_new_mtl) * 0.8;
if size( bsdf,1 ) == no_new_mtl
    mtl_color = bsdf(:,1:3)';
end

% Write materials to h_mesh
if no_new_mtl ~= 0
    h_mesh.mtl_name = cat(2, h_mesh.mtl_name, mtl_name' );
    h_mesh.mtl_color(:,no_existing_mtl+1:end) = mtl_color;
    h_mesh.mtl_thickness(:,no_existing_mtl+1:end) = ones(1,no_new_mtl) * 0.1;

    for n = 1 : no_new_mtl
        i_csv = csv_ind( find(mtl_ind == uint64(n),1) );
        if i_csv ~= 0
            h_mesh.mtl_prop(:,no_existing_mtl+n) = [ csv_prop.a(i_csv); csv_prop.b(i_csv); csv_prop.c(i_csv); csv_prop.d(i_csv); csv_prop.att(i_csv) ];
        end
    end
end
mtl_ind = mtl_ind + uint64(no_existing_mtl);

% Write mesh
no_exisiting_vert = uint64( h_mesh.no_vert );
no_existing_face  = uint64( h_mesh.no_face );
no_existing_obj   = uint64( h_mesh.no_obj );

h_mesh.vert = [ h_mesh.vert, vert_list' ];
h_mesh.face = [ h_mesh.face,  face_ind' + no_exisiting_vert ];

% Write objects to mesh
if ~isempty( obj_name )
    h_mesh.obj_name = cat( 2, h_mesh.obj_name, obj_name' );
    h_mesh.obj_index( no_existing_face+1:end ) = obj_ind' + no_existing_obj;
end

% Write material index
if ~isempty( mtl_ind )
    h_mesh.mtl_index( no_existing_face+1:end ) = mtl_ind';
end

% Reset sib-mesh index
h_mesh.Psub_mesh_index = [];

end
