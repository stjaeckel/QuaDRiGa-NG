function testLayout_3GPP_baseline_multifreq

s = qd_simulation_parameters;               % Set general simulation parameters
s.center_frequency = [ 6e9, 30e9, 70e9 ];   % Set center frequencies for the simulations
s.use_3GPP_baseline = 1;                    
s.show_progress_bars = 0;

l = qd_layout(s);
l.no_rx = 2;
l.randomize_rx_positions(100,1.5,1.5,0);
l.set_scenario('3GPP_38.901_UMi');

c = l.get_channels;

assertEqual( size(c), [2 1 3] )
assertEqual( c(1,1,1).name, 'F01-Tx0001_Rx0001' )
assertEqual( c(1,1,2).name, 'F02-Tx0001_Rx0001' )
assertEqual( c(2,1,3).name, 'F03-Tx0001_Rx0002' )



