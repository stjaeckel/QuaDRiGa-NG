function testChannel_Constructor
% Construct channel and check defaults

c = qd_channel;
assertEqual(  c.name  ,  'New_channel' );
assertEqual(  c.version  ,  qd_simulation_parameters.version );
assertEqual(  c.no_rxant , 0 );
assertEqual(  c.no_txant , 0 );
assertEqual(  c.no_path , 0 );
assertEqual(  c.no_snap , 0 );

% Set name
c.name = 'bla';
assertEqual(  c.name  ,  'bla' );

% Parse version into fields
v = c.get_version;
assertEqual( numel(v), 4 );

% Set coeff
X = rand(10,5,30,50);
c = qd_channel(X, [], 2);
assertEqual( c.coeff, X );
assertEqual( c.delay, zeros(30,50) );
assertEqual( c.individual_delays, false );
assertEqual(  c.no_rxant , 10 );
assertEqual(  c.no_txant , 5 );
assertEqual(  c.no_path , 30 );
assertEqual(  c.no_snap , 50 );
assertEqual(  c.initial_position , 2 );

% Set coeff and delay
Y = rand(30,50);
c = qd_channel(X,Y);
assertEqual( c.coeff, X );
assertEqual( c.delay, Y );
assertEqual( c.individual_delays, false );
assertEqual(  c.initial_position , 1 );

c.individual_delays = true;
assertEqual(  size( c.delay )  ,  [10,5,30,50] );
assertEqual( c.individual_delays, true );

c.individual_delays = false;
assertEqual(  size( c.delay )  ,  [30,50] );
assertEqual( c.individual_delays, false );

c = qd_channel(X,X);
assertEqual( c.individual_delays, true );

try
    c = qd_channel(X,rand(30,60)); % This should not work
    assertTrue( false );
end

c = qd_channel(Y,Y);
assertEqual( c.individual_delays, true );
c.individual_delays = 0;
assertEqual(  size( c.delay )  ,  [1,1] );
c.individual_delays = 1;
assertEqual(  size( c.delay )  ,  [30,50] );
