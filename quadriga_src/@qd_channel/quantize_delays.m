function h_channel_quant = quantize_delays( h_channel, tap_spacing, max_no_taps,...
    i_rxant, i_txant, fix_taps, verbose ) %#ok<INUSD>
%QUANTIZE_DELAYS Fixes the path delays to a grid of delay bins
%
% Calling object:
%   Single object
%
% Description:
%   For some applications, e.g. channel emulation, it is not possible to achieve an infinite delay
%   accuracy. However, when the delays are rounded to a fixed grid of delay-bins (also refereed to
%   as "taps"), the time-evolving channel is no longer smooth. When a delay "jumps" from one delay-
%   bin to the next, e.g. when a MT is moving away from the BS, the phases in the frequency domain
%   representation of the channel will suddenly change as well. Multi-carrier communications
%   systems with closed-loop channel adaption (e.g. OFDM, WiFi, LTE, etc.) will exhibit poor
%   performance in this case. This method corrects this problem by approximating the "real" delay
%   value by two delays at a fixed spacing.  For example: when we assume that the required tap
%   spacing is 5 ns (which corresponds to a 200 MHz sample-rate) and the distance between BS and MT
%   is 10 m, the delay of the LOS path would be 33.4 ns. However, the fixed tap spacing only allows
%   values of 30 or 35 ns. This method approximates the LOS delay by two taps (one at 30 and one at
%   35 ns) and linear interpolation of the path power. Hence, in the frequency domain, the
%   transition from one tap to the next is smooth. But note: this only works when the bandwidth of
%   the communication system is significantly less than the sample-rate.
%
% Input:
%   tap_spacing
%   The spacing of the delay-bin in [seconds]. Default: 5 ns
%
%   max_no_taps
%   Limits the maximum number of taps. By default, this number is infinite. If the input is
%   provided, a mapping of paths to taps is done. If the maximum number of taps is too small to
%   export all paths, only the paths with the strongest power are exported. Interpolation is done
%   whenever possible, i.e., when there are sufficient taps.
%
%   i_rxant
%   A list of receive element indices. By default, all elements are exported.
%
%   i_txant
%   A list of transmit element indices. By default, all elements are exported.
%
%   fix_taps
%   An integer number from 0 to 3, indicating if same delays should be used for different antennas
%   or snapshots. The options are: 
%
%     0  Use different delays for each tx-rx pair and for each snapshot (default)
%     1  Use same delays for all antenna pairs and snapshots (least accurate)
%     2  Use same delays for all antenna pairs, but different delays for the snapshots
%     3  Use same delays for all snapshots, but different delays for each tx-rx pair
%
%   verbose
%   Not used. The argument is kept for compatibility with previous versions.
%
% Output:
%   h_channel_quant
%   A qd_channel object containing the approximated delays and channel coefficients
%
%
% QuaDRiGa Copyright (C) 2011-2020
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

% Test if we have a scalar channel object
if numel( h_channel ) > 1
    error('QuaDRiGa:qd_channel:quantize_delays','"quantize_delays" is only defined for scalar objects.')
else
    h_channel = h_channel(1,1); % workaround for octave
end

if ~exist( 'tap_spacing' , 'var' ) || isempty( tap_spacing )
    tap_spacing = 5e-9;
end

if ~exist( 'max_no_taps' , 'var' ) || isempty( max_no_taps ) || isinf( max_no_taps )
    max_no_taps = 0; % Quadriga-Lib uses 0 for an unlimited number of taps
end

if ~exist( 'i_txant' , 'var' ) || isempty( i_txant )
    i_txant = uint32( 1:h_channel.no_txant );
else
    i_txant = uint32( i_txant );
end

if ~exist( 'i_rxant' , 'var' ) || isempty( i_rxant )
    i_rxant = uint32( 1:h_channel.no_rxant );
else
    i_rxant = uint32( i_rxant );
end

if ~exist( 'fix_taps' , 'var' ) || isempty( fix_taps )
    fix_taps = 0;
end

% Read coefficients and delays for the selected antennas
if h_channel.individual_delays
    delay = h_channel.Pdelay( i_rxant, i_txant, :, : );
else % Delays are identical on all MIMO links, size [ 1, 1, n_path, n_snap ]
    delay = reshape( h_channel.Pdelay, 1, 1, h_channel.no_path, h_channel.no_snap );
end
coeff = h_channel.Pcoeff( i_rxant, i_txant, :, : );

% Map the delays to the tap grid, a power exponent of 0.5 interpolates the path power linearly
[ coeff_re, coeff_im, delay ] = quadriga_lib.quantize_delays( real(coeff), imag(coeff), delay, ...
    tap_spacing, max_no_taps, 0.5, fix_taps );

% Delays that are shared by all antennas are returned with size [ 1, 1, n_taps, n_snap ]
if size( delay,1 ) ~= size( coeff_re,1 ) || size( delay,2 ) ~= size( coeff_re,2 )
    delay = repmat( delay, [ size(coeff_re,1), size(coeff_re,2), 1, 1 ] );
end

% Create output channel object
h_channel_quant = qd_channel( complex( coeff_re, coeff_im ), delay );

% Copy remaining data from the input channel
h_channel_quant.name = h_channel.name;
h_channel_quant.center_frequency = h_channel.center_frequency;
h_channel_quant.par = h_channel.par;
h_channel_quant.tx_position = h_channel.tx_position;
h_channel_quant.rx_position = h_channel.rx_position;

end
