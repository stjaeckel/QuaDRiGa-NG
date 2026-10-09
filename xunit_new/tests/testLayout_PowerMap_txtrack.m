function testLayout_PowerMap_txtrack
%%

a = qd_arrayant('custom',10,10,0.1);
gain = a.calc_gain;

l = qd_layout;
l.tx_array = a.copy;
l.tx_array.rotate_pattern(90,'z');
l.tx_array.rotate_pattern(-45,'x');
l.tx_position = [0;0;25];
l.tx_track(1,1).orientation = [0;0;0];

[map1 , posx , posy ] = power_map( l , 'LOSonly' , 'quick' , 5, -100,100,-100,100 );

% Check if gain is included correctly in the map
ix = posx < 0.5 & posx > -0.5;
iy = posy < 25.5 & posx > 24.5;
assertTrue( abs( 10*log10(map1{1}(iy,ix)) - gain ) < 1e-4 );

l.tx_array = a;
l.tx_track(1,1).orientation = [0;-45*pi/180;90*pi/180];
map2 = power_map( l , 'LOSonly' , 'quick' , 5, -100,100,-100,100 );

assertTrue( abs( 10*log10(map2{1}(iy,ix)) - gain ) < 1e-4 );


assertTrue(  sum(abs( map1{1}(:) - map2{1}(:) ) < 1e-3) / numel(map2{1}) > 0.99 )


