function testLayout_Randomize
%%
l = qd_layout;
assertEqual(  l.rx_name  ,  {'Rx0001'});

l.no_rx = 3;
assertEqual(  numel(l.rx_track)  , 3 );
assertEqual(  l.rx_name  ,  {'Rx0001','Rx0002','Rx0003'});

l.randomize_rx_positions;
assertEqual(  l.rx_name  ,  {'Rx0001','Rx0002','Rx0003'});
assertTrue( all( abs( l.rx_position(:) ) ~= 0 ) );

% Set-function for the rx-positoins
l.rx_position = ones(3); % write
rxpos = l.rx_position; % read
assertEqual( rxpos  , ones(3) );

l.randomize_rx_positions(100);
l.randomize_rx_positions(100,2);
l.randomize_rx_positions(100,2,2);
l.randomize_rx_positions(100,2,2,0);
l.randomize_rx_positions(100,2,2,5);
l.randomize_rx_positions(100,2,2,1,[],100,[pi/4,0,pi/2]');

assertTrue( all( abs( l.rx_position(1,:).^2 +l.rx_position(2,:).^2 - 100^2 ) < 1e-11 ) )

assertTrue( all( l.rx_track(1,1).orientation(1,:) - pi/4 < 1e-13 ) )
assertTrue( all( l.rx_track(1,1).orientation(3,:) - pi/2 < 1e-13 ) )

assertTrue( all( l.rx_track(1,3).orientation(1,:) - pi/4 < 1e-13 ) )
assertTrue( all( l.rx_track(1,3).orientation(3,:) - pi/2 < 1e-13 ) )

l.randomize_rx_positions(100,2,2,0,[],30,-pi/2);
assertTrue( l.rx_track(1,1).orientation(3,:) - pi/2 < 1e-13 )
assertTrue( l.rx_track(1,3).orientation(3,:) - pi/2 < 1e-13 )

