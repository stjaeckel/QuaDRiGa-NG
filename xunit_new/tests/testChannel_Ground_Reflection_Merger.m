function testChannel_Ground_Reflection_Merger
%%

C = ones(1,1,4,30);
C(:,:,1,:) = 4;
C(:,:,2,:) = 3;
C(:,:,3,:) = 2;
C(:,:,4,:) = 1;

D = ones(1,1,4,100);
D(:,:,1,:) = 1;
D(:,:,2,:) = 2;
D(:,:,3,:) = 3;
D(:,:,4,:) = 4;

P = ones(3,30);

c(1,1) = qd_channel( C(:,:,:,1:20),D(:,:,:,1:20) );
c(1,2) = qd_channel( C,D(:,:,:,11:40) );
c(1,3) = qd_channel( C,D(:,:,:,31:60) );
c(1,4) = qd_channel( C,D(:,:,:,51:80) );
c(1,5) = qd_channel( C,D(:,:,:,71:100) );

c(1,1).rx_position = P(:,1:20);
c(1,2).rx_position = 2*P;
c(1,3).rx_position = 3*P;
c(1,4).rx_position = 4*P;
c(1,5).rx_position = 2*P;

c(1,2).initial_position = 11;
c(1,3).initial_position = 11;
c(1,4).initial_position = 11;
c(1,5).initial_position = 11;

c(1,1).name = 'sc_tx1_rx1_seg1';
c(1,2).name = 'sc_tx1_rx1_seg2';
c(1,3).name = 'sc_tx1_rx1_seg3';
c(1,4).name = 'sc_tx1_rx1_seg4';
c(1,5).name = 'sc_tx1_rx1_seg5';

c(1,1).par(1).has_ground_reflection = 0;
c(1,2).par(1).has_ground_reflection = 1;
c(1,3).par(1).has_ground_reflection = 1;
c(1,4).par(1).has_ground_reflection = 0;
c(1,5).par(1).has_ground_reflection = 0;

c = c( randperm(5) );

d = merge(c,1,0);
assertEqual( d.no_snap ,100 )

d.individual_delays = 0;

cf = permute( d.coeff , [4,3,2,1] );
dl = d.delay.';

p = cf.^2 ./ (sum(cf.^2,2) * ones( 1,d.no_path ));
ds = sqrt( sum(p.*dl.^2,2) - sum((p.*dl),2).^2 );   

% Coefficients should be merged such that the DS remains the the same
assertTrue( all( abs( ds-ds(1) ) < 1e-14 ) );

% LOS should be the same for all channels
assertTrue( all( abs( cf(:,1) - 4 ) < 1e-14 ) );

% GR
assertTrue( all( abs( cf(1:10,2) ) < 1e-14 ) );
assertTrue( all( abs( cf(61:100,2) ) < 1e-14 ) );
assertTrue( all( abs( cf(21:50,2) - 3 ) < 1e-14 ) );
assertTrue( all( abs( dl(11:60,2) - 2 ) < 1e-14 ) );

