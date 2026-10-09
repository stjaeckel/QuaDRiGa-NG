function testChannel_copy

c = qd_channel( randn( 2,12,10,1 ), ones(10,1)  );
c.name = 'bla';
c.version = '1.1.1-1';
c.par.test = [1 2 3];
c.tx_position = [0,1,2]';
c.rx_position = [0,2,4]';
c(1,1,2) = qd_channel( randn( 2,12,10,2 ) , ones(10,2)*2 );
c(1,1,2).initial_position = 2;

d = copy( c );
assertEqual( size(c), size(d) );
assertEqual( d(1,1,1).no_snap , 1 );
assertEqual( size( d(1,1,1).delay) , [10 1] );
assertEqual( d(1,1,1).delay(1) , 1 );
assertEqual( d(1,1,2).delay(1) , 2 );
assertEqual( c(1,1,1).par.test , [1 2 3] );
assertEqual( c(1,1,1).name , 'bla' );
assertEqual( c(1,1,1).version , '1.1.1-1' );
assertEqual( c(1,1,2).initial_position , 2 );
assertEqual( c(1,1,1).tx_position , [0,1,2]' );
assertEqual( c(1,1,1).rx_position , [0,2,4]' );

