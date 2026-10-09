function testLayout_KML_write_read

s = qd_simulation_parameters;
s.center_frequency = [1e9, 2e9];
s.use_absolute_delays = 1;
s.autocorrelation_function = 'Exp150';
s.use_random_initial_phase = 0;
s.show_progress_bars = 0;
s.sample_density = 0.03;

l = qd_layout(s);
l.name = 'MiMuMa';
l.update_rate = 0.200001;
l.tx_position = [10;10;7];
l.rx_track = qd_track('linear',50.00001,pi/2);
l.rx_track.orientation = rand(3,1)*ones(1,2);
l.rx_track.set_speed(1);
l.tx_track.orientation = rand(3,1);
l.rx_position = [-10;10;1.5];
l.tx_name = {'MyfancyTx'};
l.rx_name = {'MyfancyRx1'};

l.no_rx = 2;
l.rx_track(1,2) = qd_track('circular',2,0);
l.rx_track(1,2).orientation(3,:) = 0;
l.rx_track(1,2).initial_position = [20,5,0.5]';
l.rx_track(1,2).name = 'CircleCircle';
l.rx_track(1,2).set_speed(0.04);

l.set_scenario('3GPP_38.901_UMa_NLOS');
l.tx_array = qd_arrayant('xpol');
l.tx_array.set_grid( (-180:10:180)*pi/180, (-90:10:90)*pi/180 );
l.rx_array = qd_arrayant('dipole');
set_grid( l.rx_array, (-180:5:180)*pi/180, (-90:5:90)*pi/180 );

l.layout2kml('test.kml',[12,52],1,0,[5,20,10,0]);   % Write KML
k = qd_layout.kml2layout('test.kml');               % Read KML

% Compare
assertEqual( l.name, k.name );
assertEqual( l.simpar.center_frequency, k.simpar.center_frequency );
assertEqual( l.simpar.use_absolute_delays, k.simpar.use_absolute_delays );
assertEqual( l.simpar.use_3GPP_baseline, k.simpar.use_3GPP_baseline );
assertEqual( l.simpar.autocorrelation_function, k.simpar.autocorrelation_function );
assertEqual( l.simpar.show_progress_bars, k.simpar.show_progress_bars );
assertEqual( l.simpar.use_random_initial_phase, k.simpar.use_random_initial_phase );
assertEqual( l.update_rate, k.update_rate );
assertTrue( all( abs(l.tx_position - k.tx_position) < 1e-6 ) );
assertTrue( all( abs(l.tx_track(1,1).orientation - k.tx_track(1,1).orientation) < 1e-6 ) );
assertEqual( l.tx_name, k.tx_name );
assertEqual( l.rx_name, k.rx_name );
assertTrue( all( abs(l.simpar(1,1).sample_density - k.simpar(1,1).sample_density) < 1e-10 ) );
assertEqual( k.ReferenceCoord, [12 52] );

% Split segments
assertTrue( all( abs( get_length(l.rx_track) - get_length( k.rx_track )) < 1e-5 ) );
assertEqual( k.rx_track(1,1).no_snapshots, 10 );     % Start, End, 4 Segments, 4 Correct overlap
assertTrue( all( abs( k.rx_track(1,1).positions(:,1) ) < 1e-5 ) );
assertTrue( all( abs( k.rx_track(1,1).positions(1,:) ) < 1e-5 ) );
assertTrue( all( abs( k.rx_track(1,1).positions(3,:) ) < 1e-5 ) );
assertTrue( all( abs( k.rx_track(1,1).positions(2, k.rx_track(1,1).segment_index(2:end)-1 ) - [10,20,30,40] ) < 1e-4 ) );
assertTrue( all(all( abs( k.rx_track(1,1).orientation - l.rx_track(1,1).orientation(:,ones(1,10)) ) < 1e-7 )) );
assertTrue( all( strcmp( k.rx_track(1,1).scenario , '3GPP_38.901_UMa_NLOS' ) ) )

% Closed circular track with specific orientation
assertEqual( k.rx_track(1,2).no_snapshots, 129 );    % Circle
assertTrue( all(abs(l.rx_track(1,2).positions(:) - k.rx_track(1,2).positions(:)) < 1e-5) ); % Same positions
assertTrue( all(abs( k.rx_track(1,2).orientation(:) ) < 1e-5 )); % Same positions
assertTrue( k.rx_track(1,2).closed );

% Antennas
assertTrue( all( abs( l.tx_array.Fa(:) - k.tx_array.Fa(:) ) < 1e-3 ) );
assertTrue( all( abs( l.tx_array.Fb(:) - k.tx_array.Fb(:) ) < 1e-3 ) );
assertTrue( all( abs( l.rx_array(1,1).Fa(:) - k.rx_array(1,2).Fa(:) ) < 1e-3 ) );

%% Dual Mobility
l.no_rx = 1;
l.simpar(1,1).use_3GPP_baseline = 0;
l.tx_track = qd_track('linear',25,0);
l.tx_track.set_speed(0.5);
l.tx_track.interpolate('distance',1,[],[],1);
l.tx_track.orientation = rand(3,l.tx_track.no_snapshots);
l.tx_position = [0,0,1]';
l.tx_track.name = 'thx';

delete('test.kml');

l.layout2kml('test.kml',[12,52],0,1,[10,40,25,0]);   % Write KML --> use description field
k = qd_layout.kml2layout('test.kml');               % Read KML

assertEqual( exist('test.kml.qdant','file'),2 )     % QDANT file created

% Compare
assertEqual( l.name, k.name );
assertEqual( l.simpar.center_frequency, k.simpar.center_frequency );
assertEqual( l.simpar.use_absolute_delays, k.simpar.use_absolute_delays );
assertEqual( l.simpar.use_3GPP_baseline, k.simpar.use_3GPP_baseline );
assertEqual( l.simpar.autocorrelation_function, k.simpar.autocorrelation_function );
assertEqual( l.simpar.show_progress_bars, k.simpar.show_progress_bars );
assertEqual( l.simpar.use_random_initial_phase, k.simpar.use_random_initial_phase );
assertEqual( l.update_rate, k.update_rate );
assertTrue( all( abs(l.tx_position - k.tx_position) < 1e-6 ) );
assertTrue( all( abs(l.tx_track(1,1).orientation(:) - k.tx_track(1,1).orientation(:)) < 1e-6 ) );
assertEqual( l.tx_name, k.tx_name );
assertEqual( l.rx_name, k.rx_name );
assertTrue( all( abs(l.simpar(1,1).sample_density - k.simpar(1,1).sample_density) < 1e-10 ) );
assertEqual( k.ReferenceCoord, [12 52] );

% Split segments
assertTrue( all( abs(get_length(l.rx_track) - get_length(k.rx_track)) < 1e-5 ) );
assertEqual( k.rx_track.no_snapshots, 4 );     % Start, End, 1 Segment, 1 Correct overlap
assertTrue( all( abs( k.rx_track.positions(:,1) ) < 1e-5 ) );
assertTrue( all( abs( k.rx_track.positions(1,:) ) < 1e-5 ) );
assertTrue( all( abs( k.rx_track.positions(3,:) ) < 1e-5 ) );
assertTrue( all( abs( k.rx_track.positions(2, k.rx_track.segment_index(2:end)-1 ) - 25 ) < 1e-4 ) );
assertTrue( all(all( abs( k.rx_track(1,1).orientation - l.rx_track(1,1).orientation(:,ones(1,4)) ) < 1e-7 )) );
assertTrue( all( strcmp( k.rx_track.scenario , '3GPP_38.901_UMa_NLOS' ) ) )

% Antennas
assertTrue( all( abs( l.tx_array(1,1).Fa(:) - k.tx_array(1,1).Fa(:) ) < 1e-3 ) );
assertTrue( all( abs( l.tx_array(1,1).Fb(:) - k.tx_array(1,1).Fb(:) ) < 1e-3 ) );
assertTrue( all( abs( l.rx_array(1,1).Fa(:) - k.rx_array(1,1).Fa(:) ) < 1e-3 ) );

%% Channel generation

c = k.get_channels;

assertEqual( size(c), [1 1 2] );
c = c(1,1,1);
assertEqual( c.name, 'F01-thx_MyfancyRx1' );
assertEqual( c.center_frequency, 1e9 );
assertEqual( c.no_snap, 251 );
assertEqual( c.no_txant, 2 );
assertEqual( c.no_rxant, 1 );
assertTrue( all( abs( c.rx_position(2,:) - (10:0.2:60) ) < 1e-4 ) )
assertTrue( all( abs( c.tx_position(1,:) - (0:0.1:25) ) < 1e-4 ) )
assertFalse( any( isnan( c.coeff(:) ) ) )

delete('test.kml');
delete('test.kml.qdant');


