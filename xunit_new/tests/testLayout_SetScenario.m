function testLayout_SetScenario
%%
l = qd_layout;
l.simpar.show_progress_bars = 0;
l.no_tx = 2;
l.no_rx = 2;
l.randomize_rx_positions( 100, 1.5, 1.5, 20 );
l.rx_track(1,2).segment_index = [1 300];

% Set single scenario that is supported by the builder
l.set_scenario('Freespace');
assertEqual( size( l.rx_track(1,1).scenario) , [2 1] )
assertEqual( size( l.rx_track(1,2).scenario) , [2 2] )

% Try to set a scenario that is not supported
f = @() l.set_scenario('I_do_not_exist');
assertExceptionThrown( f , 'QuaDRiGa:qd_layout:set_scenario:scenario_not_supported');

% Set complex scenario with LOS probability
l.set_scenario('3GPP_3D_UMi');
assertEqual( size( l.rx_track(1,1).scenario) , [2 1] )
assertEqual( size( l.rx_track(1,2).scenario) , [2 2] )

indoor_rx = l.set_scenario('3GPP_3D_UMa',[],[],1);
assertEqual( indoor_rx , true(1,2) );
assertEqual( size( l.rx_track(1,1).scenario) , [2 1] )
assertEqual( size( l.rx_track(1,2).scenario) , [2 2] )

l.set_scenario('mmMAGIC_UMi');
assertEqual( size( l.rx_track(1,1).scenario) , [2 1] )
assertEqual( size( l.rx_track(1,2).scenario) , [2 2] )

l.set_scenario('mmMAGIC_Indoor');
assertEqual( size( l.rx_track(1,1).scenario) , [2 1] )
assertEqual( size( l.rx_track(1,2).scenario) , [2 2] )

