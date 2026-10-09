function testSatellite_LayoutKML

s = qd_satellite;
s.semimajor_axis = 12000;
s.eccentricity = 0.3;
s.inclination = 50;
s.lon_asc_node = 300;

t = 1500:100:19500;

l = qd_layout;
l.tx_track = s.init_tracks( [],t );

% Dummy Rx Track
rx_track = qd_track([]);                                    % New track
rx_track.name = 'Rx0001';                               
rx_track.positions = rx_track.positions(:,ones(1,numel(t))); 
rx_track.movement_profile = [ 0,t(end)-t(1) ; 1,numel(t) ];
l.rx_track = rx_track;

l.set_scenario( '5G-ALLSTAR_DenseUrban_LOS' );
l.rx_track(1,1).segment_index = [1 10];
l.simpar.show_progress_bars = 0;

% Check if movement profile is correct
len = get_length( l.tx_track(1,1) );
assertTrue( abs( l.tx_track(1,1).movement_profile(2,end) - len  ) < 1e-12 )

% Plots
% s.visualize_orbit;
% s.visualize_earth( t );
% s.visualize_lonlat( t );
% l.visualize;

% Test the ccordinate transformation functions
pos = l.tx_track.positions_abs;
[ lon, lat, hnn ] = l.call_private_fcn( 'trans_ue2global', l.tx_track.positions_abs );
[ ~, ~, ~, latO, ~, lonO ] = s.orbit_predictor( t );
posN = l.call_private_fcn( 'trans_global2ue', lon, lat, hnn );

assertTrue( all(abs( lon'-lonO ) < 1e-12) )
assertTrue( all(abs( lat'-latO ) < 1e-12) )
assertTrue( all(abs( pos(:)-posN(:) ) < 1e-8) )

% Save to KML
l.layout2kml('test.kml');

% Read from KML
k = qd_layout.kml2layout( 'test.kml' );

% Compare
A = l.tx_track(1,1).positions_abs;
B = k.tx_track(1,1).positions_abs;
assertTrue( all(abs( A(:)-B(:) ) < 1e-6) )

A = l.tx_track(1,1).orientation;
B = k.tx_track(1,1).orientation;
assertTrue( all(abs( A(:)-B(:) ) < 1e-6) )

assertEqual( l.tx_track(1,1).segment_index, k.tx_track(1,1).segment_index )

delete('test.kml');


