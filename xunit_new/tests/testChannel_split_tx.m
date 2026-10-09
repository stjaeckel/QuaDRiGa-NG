function testChannel_split_tx
%%
% Generate sorted array of channel objects
G = rand( 2,12,10,3 );
D = rand(10,3) * 5 * 1e-7;
c = qd_channel(G, D);
c.name = 'Tx1_Rx1';

c(2,1) = qd_channel(G+1, D);
c(2,1).name = 'Tx1_Rx2';

c(3,1) = qd_channel(G+2, D);
c(3,1).name = 'Tx1_Rx3';

c(1,2) = qd_channel(G+3, D);
c(1,2).name = 'Tx2_Rx1';

c(2,2) = qd_channel(G+4, D);
c(2,2).name = 'Tx2_Rx2';

c(3,2) = qd_channel(G+5, D);
c(3,2).name = 'Tx2_Rx3';

cs = split_tx( c, {1:5,6:10} , {1:4,5:9,10} );

% Check oder of the outputs
assertEqual( size(cs),[3,5] )
assertEqual( cs(1,1).name , 'Tx1s1_Rx1' );
assertEqual( cs(1,2).name , 'Tx1s2_Rx1' );
assertEqual( cs(1,3).name , 'Tx2s1_Rx1' );
assertEqual( cs(2,4).name , 'Tx2s2_Rx2' );
assertEqual( cs(3,5).name , 'Tx2s3_Rx3' );

% Check correctnes of the split
assertEqual( cs(1,1).coeff , G(:,1:5,:,:) );
assertEqual( cs(1,2).coeff , G(:,6:10,:,:) );
assertEqual( cs(1,3).coeff , G(:,1:4,:,:)+3 );
assertEqual( cs(2,4).coeff , G(:,5:9,:,:)+4 );
assertEqual( cs(3,5).coeff , G(:,10,:,:)+5 );

c(3,2).name = 'Tx2_Rx4';

cs = split_tx( c, {1:5,6:10} , {1:10} );
assertEqual( size(cs),[1,9] )

c = qd_channel(G, D);
c.name = 'Bla_Tx1_Rx2';
c.individual_delays = 1;
c(1,2) =  qd_channel(G, D);
c(1,2).name = 'Bla_Tx2_Rx2';

cs = split_tx( c, {1:5} );
assertEqual( cs(1,1).name , 'Tx1s1_Rx2' );
assertEqual( cs(1,2).name , 'Tx2s1_Rx2' );
assertEqual( cs(1,1).individual_delays , true );
assertEqual( cs(1,2).individual_delays , false );

