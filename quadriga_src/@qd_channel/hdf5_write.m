function hdf5_write( h_channel, fn, layout, D1, D2, D3, D4 )
%HDF5_WRITE Save data to HDF5 file
%
% QuaDRiGa can store the generated channels in a compact HDF5-based file format which allows partial
% read and write access to the stored channels. The method 'hdf5_write' saves the current array of
% channel objects to a HDF5 file. It also allows to append new channels to an existing file which
% can significantly reduce the memory requirements to run QuaDRiGa. There are several options for
% using this method:
%
% HDF5_WRITE( fn )
%   This saves the complete array of channel objects to the file. You can reshape the array before
%   saving such that e.g. the rows represent the receivers and the columns represent the
%   transmitters. This allows you to easily load only parts of the file later on.
%
%
% HDF5_WRITE( fn, layout, D1, D2, D3, D4 )
%   This writes the channels at the specified positions in the file. The parameter "layout" defines
%   the dimensions of the channel object in the file. The parameters  'D1 ... D4" specify the
%   indices where the data should be written. For example: Your simulation has 100 receivers and 5
%   transmitters. You can create the channels for each transmitter and each receiver one after the
%   other. If you want to save the channel Tx3-Rx20 to the file, you can use
%
%       hdf5_write( fn, [100,5], 20, 3 )
%
%   If you call this function for the first time, the file will be created. You can also write
%   multiple channels at the same time. For example: If you have a 1x2 array containing the channels
%   for Tx3-Rx20 and Tx4-Rx20, you can write them to the file with
%
%       hdf5_write( fn, [100,5], 20, [3,4] )
%
%
% Note: Once the layout is written to the file, the total number of channels that can be stored in
% this file is fixed. However, you can reshape the layout by using:
%   quadriga_lib.hdf5_reshape_layout(fn, layout)
%
%
% Input:
%   filename
%       The path to the HDF5 file.
%
%   layout [optional]
%       Defines the storage layout in the file. If no layout is specified, the defualt is set as
%       follows: 
%
%       [65536,1,1,1] or [numel(channel_obj),1,1,1] if the channel array has only one dimension
%       [1024,64,1,1] if the channel array has two dimensions
%       [256,16,16,1] if the channel array has three dimensions
%       [128,8,8,8] if the channel array has four dimensions
%
%   D1 - D4 [optional]
%       Location in the file. Number of locations must match the size of the channel array.
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

if ~exist( 'fn','var' ) 
    error('QuaDRiGa:Channel:wrongInputValue','??? Filename is missing')
end

if exist( 'layout','var' ) && isnumeric(layout) && ~isempty(layout)
    if numel(layout) > 4
        error('QuaDRiGa:Channel:wrongInputValue','??? Storage layout cannot have more that 4 dimensions')
    elseif numel(layout) < 4
        layout = [layout(:); ones(4-numel(layout),1)]';
    end
elseif exist( 'layout','var' ) && ~isnumeric(layout)
    error('QuaDRiGa:Channel:wrongInputValue','??? Storage layout must be numeric')
else
    layout = [];
end

if isempty(layout)
    % Try loading layout from file
    layout = quadriga_lib.hdf5_read_layout(fn);

    if layout(1) == 0 % File does not exist
        layout = [ size(h_channel,1), size(h_channel,2), size(h_channel,3), size(h_channel,4) ];
        if layout(4) > 1
            layout = max( [ layout;[128,8,8,8] ], [], 1 );
        elseif layout(3) > 1
            layout = max( [ layout;[256,16,16,1] ], [], 1 );
        elseif layout(2) > 1
            layout = max( [ layout;[1024,64,1,1] ], [], 1 );
        else
            layout = max( [ layout;[65536,1,1,1] ], [], 1 );
        end
    end
end
layout = uint32(layout);

if exist( 'D1','var' )
    D1 = uint32( D1(:)' );
else
    D1 = uint32( 1 : size(h_channel,1) );
end

if exist( 'D2','var' )
    D2 = uint32( D2(:)' );
else
    D2 = uint32( 1 : size(h_channel,2) );
end

if exist( 'D3','var' )
    D3 = uint32( D3(:)' );
else
    D3 = uint32( 1 : size(h_channel,3) );
end

if exist( 'D4','var' )
    D4 = uint32( D4(:)' );
else
    D4 = uint32( 1 : size(h_channel,4) );
end

% Check dimensions
if size(h_channel,1) ~= numel(D1) || size(h_channel,2) ~= numel(D2) ||...
        size(h_channel,3) ~= numel(D3) || size(h_channel,4) ~= numel(D4)
    error('QuaDRiGa:Channel:wrongInputValue','??? Channel dimensions do not match requested storage dimensions');
end

% Check bounds
if any(D1 == 0) || max(D1) > layout(1) || any(D2 == 0) || max(D2) > layout(2) ||...
        any(D3 == 0) || max(D3) > layout(3) || any(D4 == 0) || max(D4) > layout(4)
    error('QuaDRiGa:Channel:wrongInputValue','??? Storage location out of bound.');
end

% Load the stored dimensions in the file
layout_file = quadriga_lib.hdf5_read_layout(fn);

% Create file
if layout_file(1) == uint32(0)
    layout_file = quadriga_lib.hdf5_create_file(fn,layout);
end
if any( layout_file ~= layout )
    error('QuaDRiGa:Channel:wrongInputValue','??? Storage layout in file does not match given layout')
end

% Save the channel objects
sic = size( h_channel );
for n = 1 : prod( sic )
    [ i1,i2,i3,i4 ] = qf.qind2sub( sic, n );
    location = [D1(i1), D2(i2), D3(i3), D4(i4)];
    quadriga_lib.hdf5_write_channel( fn, location, h_channel(i1,i2,i3,i4).par, h_channel(i1,i2,i3,i4).rx_position, ...
       h_channel(i1,i2,i3,i4).tx_position, real(h_channel(i1,i2,i3,i4).coeff), imag(h_channel(i1,i2,i3,i4).coeff), h_channel(i1,i2,i3,i4).delay,...
       h_channel(i1,i2,i3,i4).center_frequency, h_channel(i1,i2,i3,i4).name, h_channel(i1,i2,i3,i4).initial_position );
end

end
