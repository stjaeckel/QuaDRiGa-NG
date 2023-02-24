function h_array = xml_read( fn, ignore_layout )
%XML_READ Reads antenna patterns from a QDANT XML file
%
% Calling object:
%   None (static method)
%
% Description:
%   The QuaDRiGa array antenna exchange format (QDANT) is a file format used to store antenna
%   pattern data in XML and load them into QuaDRiGa. The file format specification is described in
%   the documentation. This method loads the correctly formatted XML file into a qd_arrayant object
%   array.
%
% Input:
%   fn
%   Filename of the QDANT XML file.
%
%   ignore_layout
%   Boolean value. By default (0), the layout of multiple qd_arrayant objects is stored in the
%   QDANT file. This layout is restored by xml_read. Setting ignore_layout to 1 loads all
%   qd_arrayant objects in a [1 x N] object array.
%
% Output:
%   h_array
%   Array of qd_arrayant objects.
%
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

if ~exist( 'fn','var' ) || isempty( fn )
    error('QuaDRiGa:qd_arrayant:xml_read:filename_not_given','You did not specify a filename.');
end

if ~exist( 'ignore_layout','var' ) || isempty( ignore_layout )
    ignore_layout = false;
end

[e_theta_re, e_theta_im, e_phi_re, e_phi_im, azimuth_grid, elevation_grid, element_pos, ...
    coupling_re, coupling_im, center_frequency, name, a_ind] = quadriga_lib.arrayant_qdant_read(fn,1);

a_cnt = max(a_ind(:));
if ignore_layout
    a_ind = 1:a_cnt;                                            % Option to ignore the layout
end

% Read array antennas from file
a = qd_arrayant([]);
for n = 1 : a_cnt
    if n > 1
        [e_theta_re, e_theta_im, e_phi_re, e_phi_im, azimuth_grid, elevation_grid, element_pos, ...
            coupling_re, coupling_im, center_frequency, name] = quadriga_lib.arrayant_qdant_read(fn,n);
    end
    
    % Build qd_arrayant object
    a(1,n) = qd_arrayant([]);
    a(1,n).name = name;
    a(1,n).center_frequency = center_frequency;
    a(1,n).elevation_grid = elevation_grid;
    a(1,n).azimuth_grid = azimuth_grid;
    a(1,n).no_elements = size(element_pos,2);
    a(1,n).element_position = element_pos;
    a(1,n).Fa = complex( e_theta_re, e_theta_im );
    a(1,n).Fb = complex( e_phi_re, e_phi_im );
    a(1,n).coupling = complex( coupling_re, coupling_im );
end

% Format output data
h_array = qd_arrayant([]);
for n = 1 : size( a_ind,1 )
    for m = 1 : size( a_ind,2 )
        h_array(n,m) = a(a_ind(n,m));
    end 
end

end
