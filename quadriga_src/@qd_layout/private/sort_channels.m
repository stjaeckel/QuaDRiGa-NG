function h_channel_sorted = sort_channels( h_layout, h_channel )
%SORT_CHANNELS Sorts channels to that thier order matches the order of the TX/RXs in the layout
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

no_rx = h_layout.no_rx;
no_tx = h_layout.no_tx;
n_freq = numel( h_layout.simpar(1,1).center_frequency );

if numel(h_channel) == no_rx * no_tx * n_freq
    % Reorder the channels such that they match the order in the layout
    h_channel = qf.reshapeo( h_channel, [ no_rx, no_tx, n_freq ] );
    
    rx_name = cell(1,size( h_channel,1 ));                  % Get the correct RX order
    for ir = 1 : size( h_channel,1 )
        ii = regexp( h_channel(ir,1,1).name,'_' );
        rx_name{ir} = h_channel(ir,1,1).name(ii(1)+1:end);
    end
    rx_order = zeros( 1,h_layout.no_rx );
    for ir = 1 : h_layout.no_rx
        rx_order( ir ) = find( strcmp( h_layout.rx_track(1,ir).name, rx_name ) );
    end
    
    tx_name = cell(1,size( h_channel,2 ));                  % Get the correct TX order
    for it = 1 : size( h_channel,2 )
        ii = regexp( h_channel(1,it,1).name,'_' );
        if n_freq == 1
            tx_name{it} = h_channel(1,it).name(1:ii(1)-1);
        else
            tx_name{it} = h_channel(1,it,1).name(5:ii(1)-1);
        end
    end
    tx_order = zeros( 1,h_layout.no_tx );
    for it = 1 : h_layout.no_tx
        tx_order( it ) = find( strcmp( h_layout.tx_track(1,it).name, tx_name ) );
    end
                
    h_channel_sorted = h_channel( rx_order, tx_order, 1:n_freq );
else
    h_channel_sorted = h_channel;
end

end