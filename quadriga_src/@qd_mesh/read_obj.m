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

if numel( h_mesh ) > 1
    error('QuaDRiGa:qd_mesh:read_obj','intersect not definded for object arrays.');
else
    h_mesh = h_mesh(1,1); % workaround for octave
end

use_octave = isempty( strfind( version,'R20' ) ); %#ok

if ~exist( 'fname','var' ) || isempty( fname )
    error('QuaDRiGa:qd_mesh:read_obj','Filename is not given.');
end

[ ~, mtl_prop, vert_list, face_ind, obj_ind, mtl_ind, obj_name, mtl_name ] = quadriga_lib.obj_file_read( fname );

% Read the name of the material library
i_obj = 0;
obj_mtl = {};
mtl_file_name = '';
fid = fopen(fname, 'r');        % Open File
while ~feof(fid)
    l = fgets(fid);
    if numel(l) > 7 && strcmp( l(1:6),'mtllib' ) % Look for mtllib
        mtl_file_name = sscanf( l(8:end),'%s' );
        break
    end
end
fclose( fid );

% Load material properties
no_existing_mtl = h_mesh.no_mtl;
no_new_mtl = numel(mtl_name);
mtl_color = ones(3,no_new_mtl) * 0.8;
mtl_Ni = ones(1,no_new_mtl) * 1.45;

if ~isempty( mtl_file_name )
    fid = fopen(fullfile( fileparts(fname),mtl_file_name), 'r');            % Open MTL File

    i_mtl = 0;
    while ~feof(fid)
        l = fgets(fid);
        if l(1) == 'n' && numel(l) > 7 && strcmp( l(1:6),'newmtl' ) % Look for mtllib
            name = sscanf( l(8:end),'%s' );
            ind = strcmp( mtl_name, name );
            if sum(ind) == 1
                i_mtl = find(ind);
            else
                i_mtl = 0;
            end
        elseif l(1) == 'K' && l(2) == 'd' && i_mtl ~= 0
            mtl_color(:,i_mtl) = sscanf( l,'Kd %f %f %f' );
        elseif l(1) == 'N' && l(2) == 'i' && i_mtl ~= 0
            mtl_Ni(i_mtl) = sscanf( l,'Ni %f' );
        end
    end
    fclose( fid );
end

% Write materials to h_mesh
if no_new_mtl ~= 0
    h_mesh.mtl_name = cat(2, h_mesh.mtl_name, mtl_name' );
    h_mesh.mtl_color(:,no_existing_mtl+1:end) = mtl_color;
    h_mesh.mtl_thickness(:,no_existing_mtl+1:end) = ones(1,no_new_mtl) * 0.1;

    for n = 1 : no_new_mtl
        h_mesh.mtl_prop(:,no_existing_mtl+n) = mtl_prop( find(mtl_ind == uint32(n),1) ,: )';
    end
end
mtl_ind = mtl_ind + uint32(no_existing_mtl);

% Write mesh
no_exisiting_vert = uint32( h_mesh.no_vert );
no_existing_face  = uint32( h_mesh.no_face );
no_existing_obj   = uint32( h_mesh.no_obj );

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
