function testChannel_split_snap
%%

G = rand( 2,3,4,20 );
D = rand(4,20) * 5 * 1e-7;
c = qd_channel(G, D);
c.name = 'Tx1_Rx2';

c.rx_position = [ 1:20 ; 1:20 ; 1:20 ];
c.tx_position = (1:3)';

tmp = [1;1]*(1:10);
c.par.cluster_ind = tmp(:)';
c.par.bla = 1:10;

d = c.split_snap( 1:10, 11:2:19 );

assertEqual( d(1,1).coeff, G(:,:,:,1:10) )
assertEqual( d(1,2).coeff, G(:,:,:,11:2:19) )

assertEqual( d(1,1).delay, D(:,1:10) )
assertEqual( d(1,2).delay, D(:,11:2:19) )

assertEqual( d(1,1).tx_position, (1:3)' )
assertEqual( d(1,2).tx_position, (1:3)' )

assertEqual( d(1,1).rx_position, [1;1;1]*(1:10) )
assertEqual( d(1,2).rx_position, [1;1;1]*(11:2:19) )


assertEqual( d(1,1).par.cluster_ind, [1 1 2 2 3 3 4 4 5 5] )
assertEqual( d(1,2).par.cluster_ind, [6 7 8 9 10] )

assertEqual( d(1,1).par.bla, 1:5 )
assertEqual( d(1,2).par.bla, 6:10 )

c.individual_delays = 1;

d = c.split_snap( 1:10, 11:2:19 );

assertEqual( permute( d(1,1).delay(1,1,:,:),[3,4,1,2] ), D(:,1:10) )
assertEqual( permute( d(1,2).delay(1,1,:,:),[3,4,1,2] ), D(:,11:2:19) )

