function xml_write( h_array, fn, id )
%XML_WRITE Writes antenna patterns into a QDANT XML file
%
% Calling object:
%   Object array
%
% Description:
%   The QuaDRiGa array antenna exchange format (QDANT) is a file format used to store antenna
%   pattern data in XML. The file format specification is described in the documentation. This
%   method saves a qd_arrayant object array to a XML file.
%
% Input:
%   fn
%   Filename of the QDANT XML file.
%
%   id
%   Integer number defining the array antenna ID (optional). If multiple array antennas are stored
%   in the same file, each antenna must be identified by an unique ID.
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

if ~exist( 'id','var' ) || isempty( id )
    id = 1;
end

% Process array indices
a = qd_arrayant([]);                                            % Empy array of qd_arrayant objects
a_cnt = 0;                                                      % Counter
a_ind = zeros( size( h_array) );                                % Index list for the array antennas
for ifrq = 1 : size( h_array,1 )
    for iue = 1 : size( h_array,2 )
        t_array = h_array( ifrq,iue );                          % Copy handle
        iseq = qf.eqo( t_array, a );                            % Check if antenna already exists
        if any( iseq )
            a_ind( ifrq,iue ) = find( iseq );                   % Save index
        else
            a_cnt = a_cnt + 1;                                  % Increase counter
            a(1,a_cnt) = t_array;                               % Copy handle
            a_ind( ifrq,iue ) = a_cnt;                          % Save index
        end
    end
end

for n = 1 : a_cnt
    
    if a_cnt > 1
        id = n;
    end
    
    azimuth_grid    = double( a(1,n).azimuth_grid );
    elevation_grid  = double( a(1,n).elevation_grid );
    e_theta_re      = double( real( a(1,n).Fa ));
    e_theta_im      = double( imag( a(1,n).Fa ));
    e_phi_re        = double( real( a(1,n).Fb ));
    e_phi_im        = double( imag( a(1,n).Fb ));
    element_pos     = double( a(1,n).element_position );
    coupling_re     = double( real( a(1,n).coupling ));
    coupling_im     = double( imag( a(1,n).coupling ));
    center_freq     = double( a(1,n).center_frequency );
    name            = a(1,n).name;
    
    if n < a_cnt
        id_file = quadriga_lib.arrayant_qdant_write( e_theta_re, e_theta_im, e_phi_re, e_phi_im, ...
            azimuth_grid, elevation_grid, element_pos, coupling_re, coupling_im, center_freq, name,...
            fn, id);
    else
        id_file = quadriga_lib.arrayant_qdant_write( e_theta_re, e_theta_im, e_phi_re, e_phi_im, ...
            azimuth_grid, elevation_grid, element_pos, coupling_re, coupling_im, center_freq, name,...
            fn, id, a_ind);
    end
    if id_file ~= id
        error('QuaDRiGa:qd_arrayant:xml_write','Something is wrong.');
    end
end

end

