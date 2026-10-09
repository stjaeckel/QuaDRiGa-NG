function h_channel_quant = quantize_delays_new( h_channel, tap_spacing, max_no_taps,...
    i_rxant, i_txant, fix_taps, ~ )
%QUANTIZE_DELAYS_NEW Fixes the path delays to a grid of delay bins (quadriga_lib backend)
%
%   Drop-in replacement for qd_channel.quantize_delays that delegates the
%   heavy lifting to the MEX function quadriga_lib.quantize_delays.
%
%   See original quantize_delays.m for full documentation.

% --- Input validation (same as original) ------------------------------------

if numel( h_channel ) > 1
    error('QuaDRiGa:qd_channel:quantize_delays',...
        '"quantize_delays" is only defined for scalar objects.')
else
    h_channel = h_channel(1,1);
end

if ~exist( 'tap_spacing' , 'var' ) || isempty( tap_spacing )
    tap_spacing = 5e-9;
end

if ~exist( 'max_no_taps' , 'var' ) || isempty( max_no_taps )
    max_no_taps = Inf;
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

% --- Ensure individual delays -----------------------------------------------

had_individual = h_channel.individual_delays;
if ~had_individual
    h_channel = copy( h_channel );
    h_channel.individual_delays = true;
end

% --- Extract raw data from channel object -----------------------------------

coeff = h_channel.Pcoeff( i_rxant, i_txant, :, : );   % Complex coefficients
delay = h_channel.Pdelay( i_rxant, i_txant, :, : );   % Delays [s]

% Extract real/imaginary parts; quadriga_lib accepts any numeric type
coeff_re = real( coeff );
coeff_im = imag( coeff );

% --- Check if shared delays can be used ------------------------------------
%   If the original channel had shared delays, pass [1,1,n_path,n_snap]

if ~had_individual
    delay = delay(1,1,:,:);
end

% --- Map max_no_taps for the MEX API (0 = unlimited) -----------------------

if isinf( max_no_taps )
    mex_max_taps = 0;          % quadriga_lib convention: 0 = unlimited
else
    mex_max_taps = max_no_taps;
end

% --- Call quadriga_lib MEX --------------------------------------------------

[ coeff_re_q, coeff_im_q, delay_q ] = quadriga_lib.quantize_delays( ...
    coeff_re, coeff_im, delay, tap_spacing, mex_max_taps, 0.5, fix_taps );

% --- Build output channel object --------------------------------------------

h_channel_quant = qd_channel( complex( coeff_re_q, coeff_im_q ), delay_q );

% Copy metadata
h_channel_quant.name             = h_channel.name;
h_channel_quant.center_frequency = h_channel.center_frequency;
h_channel_quant.par              = h_channel.par;
h_channel_quant.tx_position      = h_channel.tx_position;
h_channel_quant.rx_position      = h_channel.rx_position;

end
