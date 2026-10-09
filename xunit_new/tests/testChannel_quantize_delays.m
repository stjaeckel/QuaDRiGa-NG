function testChannel_quantize_delays
%%

% Generate test channel object
C = randn( 2,3,4,5 ) + 1j*randn(2,3,4,5);
C(1,2,:,1) = [1,2,3,4];

D = (1:(2*3*4*5))*2.5e-9;
D = reshape( D,2,3,4,5);

c = qd_channel(C,D);
c.name = 'test';
c.par.test = 1.5;

c.tx_position = [1;2;3];
c.rx_position = rand(3,5);
c.center_frequency = 1e9;

d = c.quantize_delays([],[],[],[],0,0);

PC = abs( c.coeff ).^2;         % Power of the coefficients
PC = sum( PC,3 );               % Sum over all taps

PD = abs( d.coeff ).^2;         % Power of the coefficients
PD = sum( PD,3 );               % Sum over all taps

assertTrue( all(abs( PC(:) - PD(:) ) < 1e-5) );

assertEqual( c.tx_position, d.tx_position );
assertEqual( c.rx_position, d.rx_position );
assertEqual( c.par.test, d.par.test );
assertEqual( c.name, d.name );
assertEqual( c.center_frequency, d.center_frequency );

assertEqual( d.no_rxant, 2 );
assertEqual( d.no_txant, 3 );
assertEqual( d.no_path, 8 );
assertEqual( d.no_snap, 5 );

x = mod( d.delay(:)./5e-9,1 );
assertTrue( all( (x > -1e-5 & x < 1e-5) | (x > 1-1e-5 & x < 1+1e-5) ) ) % Single precision

assertTrue( abs( c.coeff(2,1,1,1) - d.coeff(2,1,1,1) ) < 1e-5 ) % Same path
assertTrue( abs( d.coeff(1,1,1,1) - d.coeff(1,1,2,1) ) < 1e-5 ) % One path approximated by 2 equal taps

% Identical coefficients if delays lay on sampling grid
x = c.coeff(2,1,:,:) - d.coeff(2,1,1:4,:);
assertTrue( all( all( abs( x(:) ) < 1e-5 ) )); % Single precision

d = c.quantize_delays([],4,[],[],0,0);
assertEqual( d.no_path, 4 );

% Coefficients should be identical 
x = c.coeff - d.coeff;
assertTrue( all( all( abs( x(:) ) < 1e-5 ) )); % Single precision

% Test already quantized channels
e = d.quantize_delays([],3,[],2,0,0);
assertEqual( e.no_path, 3 );
assertEqual( e.no_rxant, 2 );
assertEqual( e.no_txant, 1 );
x = e.coeff(1,1,:,1);
assertEqual( x(:), single([2;3;4]) );

% If there are less taps than paths, the trongest are returned
d = c.quantize_delays([],3,1,2,0,0);
x = d.coeff(1,1,:,1);
assertTrue( all( abs( x(:) - [2;3;4] ) < 1e-5 ) );

% If there are more taps than paths, the strongest shoule be split
d = c.quantize_delays([],5,1,2,0,0);
x = d.coeff(1,1,:,1);
assertTrue( all( abs( x(:) - [1;2;3;sqrt(0.5)*[4;4]] ) < 1e-5 ) );

d = c.quantize_delays([],6,1,2,0,0);
x = d.coeff(1,1,:,1);
assertTrue( all( abs( x(:) - [1;2;sqrt(0.5)*[3;3;4;4]] ) < 1e-5 ) );

% Test fixed delays
d = c.quantize_delays(2.5e-9,[],[],[],1,0);
assertEqual( numel( D ), d.no_path );          
assertEqual( sum(d.coeff(1,1,:,1)~=0), 4 );

d = c.quantize_delays(2.5e-9,[],[],[],2,0);     % Fixed delays for all antennas
assertEqual( d.no_path, 2*3*4 );   

d = c.quantize_delays(2.5e-9,[],[],[],3,0);     % Fixed delays for all snapshots
assertEqual( d.no_path, 5*4 );  





