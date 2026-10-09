function testChannel_Merger
%%

warning('off','QuaDRiGa:qd_channel:merge')

C = ones(1,1,2,100);
D = ones(1,1,2,100);
P = ones(3,100);
D(1,1,1,:) = 0:99;
D(1,1,2,:) = 50:149;

c = qd_channel;
c(1,1) = qd_channel( C,D );
c(1,1).rx_position = P;
c(1,2) = qd_channel( 2*C,D+50 );
c(1,2).rx_position = 2*P;
c(1,3) = qd_channel( cat(3,3*C,4*C),cat(3,D+100,D+400) );
c(1,3).rx_position = 3*P;
c(1,4) = qd_channel( C,D );

c(1,2).initial_position = 51;
c(1,3).initial_position = 51;

c(1,1).name = 'sc_tx1_rx1_seg1';
c(1,2).name = 'sc_tx1_rx1_seg2';
c(1,3).name = 'sc_tx1_rx1_seg3';
c(1,4).name = 'a_b';

d = merge(c,0.2,0);

assertEqual( size(d) , [1,2] );                 % Vector of channel objects
assertEqual( d(1,1).name , 'a_b' );             % Alphabetic order
assertEqual( d(1,2).name , 'tx1_rx1' );         % Name string processing

assertEqual( size(d(1,2).coeff) , [1,1,5,200] );       % Merged coefficients
cf = reshape( d(1,2).coeff(1,1,1,:) , 1, []);
dl = reshape( d(1,2).delay(1,1,1,:) , 1, []);

assertTrue( all(abs( cf(1:89) - ones(1,89) ) < 1e-12) );
assertTrue( all(abs( cf(90:100) - (1+sin(pi/2*(1:11)/12).^2) ) < 1e-12) ); % LOS ramp
assertTrue( all(abs( cf(101:139) - ones(1,39)*2 ) < 1e-12) );
assertTrue( all(abs( cf(140:150) - (2+sin(pi/2*(1:11)/12).^2) ) < 1e-12) ); % LOS ramp
assertTrue( all(abs( cf(151:200) - ones(1,50)*3 ) < 1e-12) );
assertTrue( all(abs( dl - (0:199) ) < 1e-12) );

cf = reshape( d(1,2).coeff(1,1,2,:) , 1, []);
dl = reshape( d(1,2).delay(1,1,2,:) , 1, []);
assertTrue( all(abs( cf(1:89) - ones(1,89) ) < 1e-12) );
assertTrue( all(abs( dl(1:100) - (50:149) ) < 1e-12) );

cf = reshape( d(1,2).coeff(1,1,3,:) , 1, []);
dl = reshape( d(1,2).delay(1,1,3,:) , 1, []);
assertTrue( all(abs( cf(101:139) - ones(1,39)*2 ) < 1e-12) );
assertTrue( all(abs( dl(90:150) - (139:199) ) < 1e-12) );

