function testBuilder_check_dual_mobility
%%

warning('off','QuaDRiGa:qd_builder:check_dual_mobility:no_rx_antenna')
warning('off','QuaDRiGa:qd_builder:check_dual_mobility:no_tx_antenna')

b = qd_builder;
b.check_dual_mobility;
assertEqual( b.dual_mobility, -1 );

% Single statc tx and rx
b.rx_positions = [0;0;1];
b.tx_position = [0;0;25];
b.check_dual_mobility;
assertEqual( b.dual_mobility, false );
assertEqual( b.rx_track.initial_position, b.rx_positions );
assertEqual( b.tx_track.initial_position, b.tx_position );
assertEqual( b.tx_track.orientation, [0;0;0] );
assertTrue( isa( b.rx_array, 'qd_arrayant') );
assertTrue( isa( b.tx_array, 'qd_arrayant') );

% Rx-track only
b = qd_builder;
b.tx_position = [0;0;25];
b.rx_track = qd_track([]);
b.check_dual_mobility;
assertEqual( b.dual_mobility, false );
assertEqual( b.rx_track.initial_position, b.rx_positions );

% Two mobiles
b = qd_builder;
b.rx_positions = [0,0;0,0;1,2];
b.tx_position = [0;0;25 ];
b.check_dual_mobility;
assertEqual( b.dual_mobility, false );
assertEqual( numel( b.rx_track ), 2 );
assertEqual( numel( b.rx_array ), 2 );
assertTrue( qf.eqo( b.rx_track(1,1), b.rx_track(1,2) ) );
assertTrue( qf.eqo( b.rx_array(1,1), b.rx_array(1,2) ) );

% Two mobiles
b = qd_builder;
b.simpar.center_frequency(2) = 10e9;    % Dual-Freq.
b.rx_positions = [0,0;0,0;1,2];     % Two Rx
b.tx_position = [0;0;25 ];          % One Tx
b.rx_track = qd_track([]);          % One rx-track
b.rx_array = qd_arrayant([]);
b.rx_array(2,1) = qd_arrayant;
b.check_dual_mobility;
assertEqual( b.dual_mobility, false );
assertEqual( b.rx_track(1,1).orientation, [0;0;0] );
assertTrue( qf.eqo( b.rx_track(1,1), b.rx_track(1,2) ) );
assertTrue( qf.eqo( b.rx_array(1,1), b.rx_array(1,2) ) );
assertFalse( qf.eqo( b.rx_array(1,2), b.rx_array(2,2) ) );
assertTrue( qf.eqo( b.rx_array(1,1), b.rx_array(1,2) ) );
assertTrue( qf.eqo( b.tx_array(1,1), b.tx_array(2,2) ) );

% Two mobiles, two tx
b = qd_builder;
b.simpar.center_frequency(2) = 10e9;    % Dual-Freq.
b.rx_positions = [0,0;0,0;1,2];     % Two Rx
b.tx_position = [0,0;0,0;25,26 ];          % Two Tx
b.rx_array = qd_arrayant([]);
b.rx_array(1,2) = qd_arrayant;
b.tx_array = qd_arrayant([]);
b.tx_array(2,1) = qd_arrayant;
b.check_dual_mobility;
assertEqual( b.dual_mobility, true );
assertFalse( qf.eqo( b.rx_array(2,1), b.rx_array(2,2) ) );
assertTrue( qf.eqo( b.rx_array(1,1), b.rx_array(2,1) ) );
assertFalse( qf.eqo( b.tx_array(1,2), b.tx_array(2,2) ) );
assertTrue( qf.eqo( b.tx_array(1,1), b.tx_array(1,2) ) );

% Two mobiles, two (identical) tx - Option 1
b = qd_builder;
b.simpar.center_frequency(2) = 10e9;    % Dual-Freq.
b.rx_positions = [0,0;0,0;1,2];     % Two Rx
b.tx_position = [0,0;0,0;25,25 ];          % Two Tx
b.tx_array = qd_arrayant([]);
b.tx_array(1,2) = qd_arrayant;
b.check_dual_mobility;
assertEqual( b.dual_mobility, false );
assertFalse( qf.eqo( b.tx_array(2,1), b.tx_array(2,2) ) );
assertTrue( qf.eqo( b.tx_array(1,1), b.tx_array(2,1) ) );

% Two mobiles, two (identical) tx - Option 2
b = qd_builder;
b.rx_positions = [0,0;0,0;1,2];     % Two Rx
b.tx_track = qd_track( 'linear',0,pi/2 );
b.tx_track.initial_position = [0;0;25];
b.tx_track(1,2) = qd_track( 'linear',0,0 );
b.tx_track(1,2).initial_position = [0;0;25];
b.check_dual_mobility;
assertEqual( b.dual_mobility, false );

% One mobile, one (mobile) tx
b = qd_builder;
b.rx_track = qd_track('linear',2,pi/2);
b.tx_track = qd_track('linear',1,pi/2);
b.rx_track.initial_position(3)=1;
b.check_dual_mobility;
assertEqual( b.dual_mobility, true );

% Two mobiles, one (mobile) tx
b = qd_builder;
b.rx_track = qd_track('linear',2,pi/2);
b.rx_track(1,2) = qd_track('linear',2,pi/2);
b.tx_track = qd_track('linear',1,pi/2);
b.tx_track.initial_position(3)=1;

b.check_dual_mobility;
assertEqual( b.dual_mobility, true );
assertTrue( qf.eqo( b.tx_track(1,1), b.tx_track(1,2) ) );

% Test for errors
b = qd_builder;
b.simpar = [];
b.rx_positions = [0;0;1];
f = @() b.check_dual_mobility;
assertExceptionThrown( f , 'QuaDRiGa:qd_builder:check_dual_mobility:no_simpar');

b = qd_builder;
b.rx_track = qd_track([]);          % One rx-track
b.rx_track(2,1) = qd_track([]);          % One rx-track
f = @() b.check_dual_mobility;
assertExceptionThrown( f , 'QuaDRiGa:qd_builder:check_dual_mobility:rx_track_rows');

b = qd_builder;
b.rx_track = qd_arrayant([]);          % One rx-track
f = @() b.check_dual_mobility;
assertExceptionThrown( f , 'QuaDRiGa:qd_builder:check_dual_mobility:wrong_Rx_track_class');

b = qd_builder;
b.rx_positions = [0,0,0;0,0,0;1,2,3];     % Two Rx
b.rx_track = qd_track([]);          % One rx-track
b.rx_track(1,2) = qd_track([]);          % One rx-track
f = @() b.check_dual_mobility;
assertExceptionThrown( f , 'QuaDRiGa:qd_builder:check_dual_mobility:rx_track_size_mismatch');

b = qd_builder;
b.rx_positions = [0;0;1];
b.rx_array = qd_track([]);
f = @() b.check_dual_mobility;
assertExceptionThrown( f , 'QuaDRiGa:qd_builder:check_dual_mobility:wrong_Rx_array_class');

b = qd_builder;
b.rx_positions = [0,0,0;0,0,0;1,2,3];     % Two Rx
b.rx_array = qd_arrayant([]);          % One rx-track
b.rx_array(1,2) = qd_arrayant([]);          % One rx-track
f = @() b.check_dual_mobility;
assertExceptionThrown( f , 'QuaDRiGa:qd_builder:check_dual_mobility:rx_array_size_mismatch');

b = qd_builder;
b.rx_positions = [0;0;1];
b.tx_position = [0;0;25];
b.tx_array = qd_track([]);
f = @() b.check_dual_mobility;
assertExceptionThrown( f , 'QuaDRiGa:qd_builder:check_dual_mobility:wrong_Tx_array_class');

b = qd_builder;
b.rx_positions = [0;0;1];
b.tx_position = [0;0;25];
b.tx_array = qd_arrayant([]);
b.tx_array(2,1) = qd_arrayant([]);
f = @() b.check_dual_mobility;
assertExceptionThrown( f , 'QuaDRiGa:qd_builder:check_dual_mobility:tx_array_size_mismatch');

b = qd_builder;
b.rx_positions = [0;0;1];
f = @() b.check_dual_mobility;
assertExceptionThrown( f , 'QuaDRiGa:qd_builder:check_dual_mobility:tx_position_undefined');

b = qd_builder;
b.rx_positions = [0,0,0;0,0,0;1,2,3];     
b.tx_track = qd_track([]);         
b.tx_track(1,2) = qd_track([]);         
f = @() b.check_dual_mobility;
assertExceptionThrown( f , 'QuaDRiGa:qd_builder:check_dual_mobility:tx_track_size_mismatch');

b = qd_builder;
b.rx_positions = [0,0,0;0,0,0;1,2,3];     
b.tx_position = [0,0;0,0;25,25 ];          % Two Tx  
f = @() b.check_dual_mobility;
assertExceptionThrown( f , 'QuaDRiGa:qd_builder:check_dual_mobility:tx_position_size_mismatch');

b = qd_builder;
b.rx_positions = [0;0;1];     
b.tx_track = qd_track([]);         
b.tx_track(2,1) = qd_track([]);   
f = @() b.check_dual_mobility;
assertExceptionThrown( f , 'QuaDRiGa:qd_builder:check_dual_mobility:tx_track_rows');

b = qd_builder;
b.rx_positions = [0;0;1];     
b.tx_track = qd_arrayant([]);         
f = @() b.check_dual_mobility;
assertExceptionThrown( f , 'QuaDRiGa:qd_builder:check_dual_mobility:wrong_Tx_track_class');

b = qd_builder;
b.rx_positions = [0,0,0;0,0,0;1,2,3];     
b.tx_position = [0,0,0;0,0,0;1,2,3];  
b.tx_track = qd_track([]);         
b.tx_track(1,2) = qd_track([]);   
f = @() b.check_dual_mobility;
assertExceptionThrown( f , 'QuaDRiGa:qd_builder:check_dual_mobility:tx_track_size_mismatch2');

b = qd_builder;
b.tx_track = qd_track('linear',1,0);    
b.rx_track = qd_track('linear',1,0);  
b.rx_track.interpolate_positions(5);
f = @() b.check_dual_mobility;
assertExceptionThrown( f , 'QuaDRiGa:qd_builder:check_dual_mobility:tx_rx_track_lenght_mismatch');

b = qd_builder;
b.tx_track = qd_track('linear',1);    
b.rx_track = qd_track('linear',1);  
f = @() b.check_dual_mobility;
assertExceptionThrown( f , 'QuaDRiGa:qd_builder:check_dual_mobility:colocated_tx_rx');

b = qd_builder;
b.tx_track = qd_track('linear',0);    
b.rx_track = qd_track('linear',1);  
f = @() b.check_dual_mobility;
assertExceptionThrown( f , 'QuaDRiGa:qd_builder:check_dual_mobility:colocated_tx_rx');