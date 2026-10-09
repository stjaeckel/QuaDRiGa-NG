function testBuilder_LOS_channels
%%
% Setup basic test scenario
b = qd_builder('3GPP_38.901_UMi_LOS');
b.rx_array = qd_arrayant('ula4');
b.tx_array = qd_arrayant('ula8');
b.name = 'Tx1';
b.simpar.show_progress_bars = false;
b.tx_position = [0;0;25];
b.rx_positions = [10,0,1.5;0,10,2.5]';
b.rx_track = qd_track.generate('linear',1,0);
b.rx_track.initial_position = b.rx_positions(:,1);
b.rx_track.name = 'Rx1';
b.rx_track(1,2) = qd_track.generate('linear',1,pi/2);
b.rx_track(1,2).initial_position = b.rx_positions(:,2);
b.rx_track(1,2).name = 'Rx2';

c = b.get_los_channels;
assertEqual( numel(c) , 1 );
assertEqual( c.no_rxant , b.rx_array.no_elements );
assertEqual( c.no_txant , b.tx_array.no_elements  );
assertEqual( c.no_snap , b.no_rx_positions );
assertEqual( c.tx_position , b.tx_position );
assertEqual( c.rx_position , b.rx_positions );
assertEqual( c.individual_delays , false );

cr = b.get_los_channels('coeff');
assertEqual( size(cr) , [ b.rx_array(1,1).no_elements , b.tx_array(1,1).no_elements , b.no_rx_positions ] );
assertEqual( permute( c.coeff, [1,2,4,3] ) , cr );

% All phases of MT1 mut be idenical (arrays aligned)
assertTrue( all( abs( reshape( cr(:,:,1),1,[] ) - cr(1,1,1) ) < 1e-12 ) );

% Rx antennas of the second MT mut have same phase, but tx must be different
assertTrue( all( reshape( abs( ones(4,1) * cr(1,:,2) - cr(:,:,2)) , 1 ,[] ) < 1e-12 ) );

c = b.get_los_channels([],[1,2,3]);
assertTrue( isa( c.coeff , 'double' ) );
assertEqual( c.no_txant , 3  );

cr = b.get_los_channels('raw');
assertEqual( size(cr) , [ 2 , b.tx_array(1,1).no_elements , b.no_rx_positions ] );
assertTrue( all( abs(cr(1,:)) - ones(1,b.no_rx_positions*b.tx_array(1,1).no_elements) < 1e-12 ));
assertTrue( all(  abs(cr(2,:)) < 1e-12  ));

b.tx_array(1,1).coupling(2,2) = 2;
cr = b.get_los_channels('raw');
assertTrue(all(abs(  reshape( abs(cr(1,2,:)) , 1, [] ) - [2 2]   )<1e-12));

b.tx_array(1,1).coupling(1,3) = 2;
cr = b.get_los_channels('raw');
assertTrue(all(abs(    cr(1,:,1) - [1 2 3 1 1 1 1 1]     )<1e-12));

b.tx_array(1,1).coupling = eye(b.tx_array(1,1).no_elements);
b.rx_array(1,1).coupling(2,2) = 2;
cr = b.get_los_channels('coeff');

% Second row must be 2*first row
assertTrue( all( reshape( abs( abs(cr(2,:,:)) - 2*abs(cr(1,:,:)) ) , 1,[] ) < 1e-12 ) );

b.tx_array(1,1).coupling(2,2) = 2;
cr = b.get_los_channels('coeff',1:2);

assertTrue( all( abs( reshape( abs(cr(:,:,1))./ abs(cr(1,1,1)) - [1,2;2,4;1,2;1,2] ,1, [] ) ) <1e-12 ) );

% Test successive builder array generation
b(1,2) = qd_builder('LOSonly');
b(1,2).tx_position = [0;0;25];
b(1,2).rx_positions = [25,0,0]';

c = get_los_channels( b );
assertEqual( size(c), size(b) );

% Test warning message
b(1,1).rx_array(1,2) = qd_arrayant('omni');

f = @() get_los_channels( b(1,1) );
assertExceptionThrown( f , 'QuaDRiGa:qf_builder:get_los_channels:Rx_array_ambiguous');


