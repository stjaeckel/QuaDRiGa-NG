function [ h_channel, layout, has_data ] = hdf5_read( fn, I1, I2, I3, I4, IS )
%HDF5_READ Load data from stored HDF5 file
%
% QuaDRiGa can store the generated channels in a compact HDF5-based file format which allows partial
% read and write access to the stored channels. The static method 'hdf5_read' allows access to these
% stored channels. In addition, it is also possible to only load channels object partially in order
% to save memory. There are several options for using this method:
%
% h_channel = HDF5_READ( fn )
%   reads the complete file with all stored channels.
%
% [ ~, dims, has_data ] = HDF5_READ( fn, 0 )
%   no channels will be loaded and "h_channel" will returned as an empty array. However, the output
%   variable "dims" will contain the dimensions of the stored channel objects.
%
% h_channel = HDF5_READ( filename, rowInd, colInd, ... )
%   "rowInd" and "colInd" (and optional "d3Ind" and "d4Ind") are given. Provided that the channel
%   objects are stored as a matrix in the HDF5 file, you can select the rows and columns to load.
%   For example, if the rows represent the receiver index and the columns represent the transmitter
%   index, you can use
%
%       h_channel = qd_channel.hdf5_read( filename, 2:3, 4 )
%
%   to load the channels for Tx4-Rx2 and Tx4-Rx3. A 2x1 array of channel objects will be returned.
%
% h_channel = HDF5_READ( filename, [], colInd, ... )
%   "rowInd" is empty (i.e. []) and "colInd" (and optional "d3Ind" and "d4Ind") are given. In this
%   case, all rows for the given "colInd" are returned. For example, if the rows represent the
%   receiver index and the columns represent the transmitter index, you can use
%
%       h_channel = qd_channel.hdf5_read( filename, [], 4 )
%
%   to load the channels for all receivers belonging to Tx4.
%
% h_channel = qd_channel.hdf5_read( filename,rowInd,colInd,[],[],snap_ind )
%   This option allows to only load selected snapshots of the channel objects.
%
% Input:
%   filename
%       The path to the HDF5 file which contains the stored channel data.
%
%   I1
%       This is an optional parameter which can be used for loading only parts of the file. See
%       description for details.
%
%   I2
%       The index of the columns of the channel objects which are loaded from the file. This
%       variable is only used if "rowInd" is given. "colInd" can also be empty. In this case, all
%       columns are returned (same as case 5 above).
%
%   I3
%       The index of the 3rd dimension of the channel objects which are loaded from the file. This
%       variable is only used if "rowInd" and "colInd" are given.
%
%   I4
%       The index of the 4th dimension of the channel objects which are loaded from the file. This
%       variable is only used if "rowInd", "colInd" and "d3Ind" are given.
%
%   IS
%       Index of the snapshots to be loaded.
%
% Output:
%   h_channel
%       An array of 'qd_channel' objects
%
%   layout
%       The storage layout of the data in the file
%
%   has_data
%       Array indicating if data was stored in that location
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

if ~exist( 'fn','var' )
    error('QuaDRiGa:Channel:wrongInputValue','??? Filename is missing')
end

% Read layout from file
[ layout, has_data ] = quadriga_lib.hdf5_read_layout(fn);
if layout(1) == 0
    error('QuaDRiGa:Channel:wrongInputValue','??? File does not exist');
end

if sum(has_data(:)) == 0 % File has no stored channels
    I1 = uint32(0);
    I2 = uint32(0);
    I3 = uint32(0);
    I4 = uint32(0);
end

if ~exist( 'I1','var' )
    I1 = [];
else
    I1 = uint32( I1(:) );
end
if ~exist( 'I2','var' )
    I2 = [];
else
    I2 = uint32( I2(:) );
end
if ~exist( 'I3','var' )
    I3 = [];
else
    I3 = uint32( I3(:) );
end
if ~exist( 'I4','var' )
    I4 = [];
else
    I4 = uint32( I4(:) );
end

u1 = isempty(I1);
u2 = isempty(I2);
u3 = isempty(I3);
u4 = isempty(I4);

h_channel = qd_channel([]);
h_channel.name = 'empty';
if u1 || I1(1) ~= 0

    % There are 16 possible combinations
    if u1 && u2 && u3 && u4             % 0000
        I1 = find( sum(sum(sum( has_data,4),3),2) ~= 0 );
        I2 = find( sum(sum(sum( has_data,4),3),1) ~= 0 );
        I3 = find( sum(sum(sum( has_data,4),2),1) ~= 0 );
        I4 = find( sum(sum(sum( has_data,3),2),1) ~= 0 );

    elseif u1 && u2 && u3 && ~u4        % 0001
        I1 = find( sum(sum(sum( has_data(:,:,:,I4),4),3),2) ~= 0 );
        I2 = find( sum(sum(sum( has_data(:,:,:,I4),4),3),1) ~= 0 );
        I3 = find( sum(sum(sum( has_data(:,:,:,I4),4),2),1) ~= 0 );

    elseif u1 && u2 && ~u3 && u4        % 0010
        I1 = find( sum(sum(sum( has_data(:,:,I3,:),4),3),2) ~= 0 );
        I2 = find( sum(sum(sum( has_data(:,:,I3,:),4),3),1) ~= 0 );
        I4 = find( sum(sum(sum( has_data(:,:,I3,:),3),2),1) ~= 0 );

    elseif u1 && u2 && ~u3 && ~u4       % 0011
        I1 = find( sum(sum(sum( has_data(:,:,I3,I4),4),3),2) ~= 0 );
        I2 = find( sum(sum(sum( has_data(:,:,I3,I4),4),3),1) ~= 0 );

    elseif u1 && ~u2 && u3 && u4        % 0100
        I1 = find( sum(sum(sum( has_data(:,I2,:,:),4),3),2) ~= 0 );
        I3 = find( sum(sum(sum( has_data(:,I2,:,:),4),2),1) ~= 0 );
        I4 = find( sum(sum(sum( has_data(:,I2,:,:),3),2),1) ~= 0 );

    elseif u1 && ~u2 && u3 && ~u4       % 0101
        I1 = find( sum(sum(sum( has_data(:,I2,:,I4),4),3),2) ~= 0 );
        I3 = find( sum(sum(sum( has_data(:,I2,:,I4),4),2),1) ~= 0 );

    elseif u1 && ~u2 && ~u3 && u4       % 0110
        I1 = find( sum(sum(sum( has_data(:,I2,I3,:),4),3),2) ~= 0 );
        I4 = find( sum(sum(sum( has_data(:,I2,I3,:),3),2),1) ~= 0 );

    elseif u1 && ~u2 && ~u3 && ~u4      % 0111
        I1 = find( sum(sum(sum( has_data(:,I2,I3,I4),4),3),2) ~= 0 );

    elseif ~u1 && u2 && u3 && u4        % 1000
        I2 = find( sum(sum(sum( has_data(I1,:,:,:),4),3),1) ~= 0 );
        I3 = find( sum(sum(sum( has_data(I1,:,:,:),4),2),1) ~= 0 );
        I4 = find( sum(sum(sum( has_data(I1,:,:,:),3),2),1) ~= 0 );

    elseif ~u1 && u2 && u3 && ~u4       % 1001
        I2 = find( sum(sum(sum( has_data(I1,:,:,I4),4),3),1) ~= 0 );
        I3 = find( sum(sum(sum( has_data(I1,:,:,I4),4),2),1) ~= 0 );

    elseif ~u1 && u2 && ~u3 && u4       % 1010
        I2 = find( sum(sum(sum( has_data(I1,:,I3,:),4),3),1) ~= 0 );
        I4 = find( sum(sum(sum( has_data(I1,:,I3,:),3),2),1) ~= 0 );

    elseif ~u1 && u2 && ~u3 && ~u4      % 1011
        I2 = find( sum(sum(sum( has_data(I1,:,I3,I4),4),3),1) ~= 0 );

    elseif ~u1 && ~u2 && u3 && u4       % 1100
        I3 = find( sum(sum(sum( has_data(I1,I2,:,:),4),2),1) ~= 0 );
        I4 = find( sum(sum(sum( has_data(I1,I2,:,:),3),2),1) ~= 0 );

    elseif ~u1 && ~u2 && u3 && ~u4      % 1101
        I3 = find( sum(sum(sum( has_data(I1,I2,:,I4),4),2),1) ~= 0 );

    elseif ~u1 && ~u2 && ~u3 && u4      % 1110
        I4 = find( sum(sum(sum( has_data(I1,I2,I3,:),3),2),1) ~= 0 );
    end

    if isempty(I1); I1 = 1; end
    if isempty(I2); I2 = 1; end
    if isempty(I3); I3 = 1; end
    if isempty(I4); I4 = 1; end
    
    I1 = uint32( I1(:) );
    I2 = uint32( I2(:) );
    I3 = uint32( I3(:) );
    I4 = uint32( I4(:) );

    if exist( 'IS','var' )
        IS = uint32( IS(:) );
    else
        IS = uint32( [] );
    end

    % Load data from file
    ow = 0; 
    for iw = I4'
        ow = ow + 1;
        oz = 0;
        for iz = I3'
            oz = oz + 1;
            oy = 0;
            for iy = I2'
                oy = oy + 1; 
                ox = 0;
                for ix = I1'
                    ox = ox + 1;
                    h_channel(ox,oy,oz,ow) = qd_channel([]);
                    [ h_channel(ox,oy,oz,ow).par, rx_position, tx_position,...
                        coeff_re, coeff_im, delay, center_frequency, h_channel(ox,oy,oz,ow).name,...
                        initial_position, path_gain, path_length, path_polarization, path_angles,...
                        path_coord, rx_orientation, tx_orientation ] = quadriga_lib.hdf5_read_channel( fn, [ix,iy,iz,iw], IS );

                    if ~isempty(coeff_re) && ~isempty(coeff_im)
                        h_channel(ox,oy,oz,ow).coeff = complex(coeff_re, coeff_im);
                    end
                    if ~isempty(delay)
                        h_channel(ox,oy,oz,ow).delay = delay;
                    end
                    if ~isempty(rx_position)
                        h_channel(ox,oy,oz,ow).rx_position = rx_position;
                    end
                    if ~isempty(tx_position)
                        h_channel(ox,oy,oz,ow).tx_position = tx_position;
                    end
                    if ~isempty(center_frequency)
                        h_channel(ox,oy,oz,ow).center_frequency = center_frequency(1);
                    end
                    h_channel(ox,oy,oz,ow).initial_position = initial_position;
                    if ~isempty(path_gain)
                        h_channel(ox,oy,oz,ow).par.path_gain = path_gain;
                    end
                    if ~isempty(path_length)
                        h_channel(ox,oy,oz,ow).par.path_length = path_length;
                    end
                    if ~isempty(path_polarization)
                        h_channel(ox,oy,oz,ow).par.path_polarization = path_polarization;
                    end
                    if ~isempty(path_angles)
                        h_channel(ox,oy,oz,ow).par.path_angles = path_angles;
                    end
                    if ~isempty(path_coord)
                        h_channel(ox,oy,oz,ow).par.path_coord = path_coord;
                    end
                    if ~isempty(rx_orientation)
                        h_channel(ox,oy,oz,ow).par.rx_orientation = rx_orientation;
                    end
                    if ~isempty(tx_orientation)
                        h_channel(ox,oy,oz,ow).par.tx_orientation = tx_orientation;
                    end
                end
            end
        end
    end
end

end
