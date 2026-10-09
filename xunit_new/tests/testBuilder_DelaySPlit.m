function testBuilder_DelaySPlit
%%
b = qd_builder('3GPP_38.901_UMi_LOS_GR');
b.scenpar.NumClusters = 5;
b.simpar.show_progress_bars = false;

a = qd_arrayant('ula2');    % V-Polarization
a.element_position(:) = 0;
b.rx_array = a;
b.tx_array = a;

d = b.simpar.wavelength * 200;     % LOS path length = 2D distance
x = b.simpar.wavelength * 200.5;   % GR path length (0.5 lambda longer)
h = 0.5*sqrt(x^2-d^2);      % Height

b.tx_position = [0;0;h];
b.rx_positions = [d,0,h ; 1e4,0,h]';

gen_parameters(b);
c = get_channels(b);
