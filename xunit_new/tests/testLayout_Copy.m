function testLayout_Copy
%% Test the copy function

a = qd_layout;
a.name = 'Bla1';
a.simpar.center_frequency = 1e6;
a.no_tx = 2;
a.no_rx = 3;
a.rx_array(1,3) = qd_arrayant('xpol');
a.tx_track(1,2) = qd_track('circular',1,1);

a(1,2) = qd_layout;         % Different layout
a(1,2).name = 'Bla2';
a(1,2).simpar(1,1).center_frequency = 2e6;

a(1,3) = a(1,1);            % Copy of handle

b = copy(a);

assertEqual( b(1,1).name  , 'Bla1'  );
assertEqual( b(1,1).simpar(1,1).center_frequency  , 1e6  );

assertEqual( b(1,2).name  , 'Bla2'  );
assertEqual( b(1,2).simpar(1,1).center_frequency  , 2e6  );

assertEqual( b(1,3).name  , 'Bla1'  );
assertEqual( b(1,3).simpar(1,1).center_frequency  , 1e6  );

assertEqual( qf.eqo( b(1,1).rx_array(1,1), b(1,1).rx_array ), [true true false] )
assertEqual( qf.eqo( b(1,3).rx_array(1,1), b(1,1).rx_array ), [true true false] )

assertEqual( qf.eqo( b(1,1).rx_array(1,1), b(1,3).rx_array ), [true true false] )

