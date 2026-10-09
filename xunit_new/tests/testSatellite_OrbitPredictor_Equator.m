function testSatellite_OrbitPredictor_Equator

% Circular LEO orbit around Equator
s = qd_satellite;
s.semimajor_axis = 6500;
s.arg_periapsis = 10;
s.true_anomaly = 80;

% Predict orbit at start point
[ xyzI, xyzR, r, lat, lonI, lonR, pq ] = s.orbit_predictor( 0, 1 );

% Inertial and rotating frame are equal at the beginning of the simulation
assertTrue( all( abs( xyzI - xyzR ) < 1e-7 ) );

% Longitude positions must be 90 degree (arg_periapsis + true_anomaly)
assertTrue( all( abs( lonI - 90 ) < 1e-7 ) );
assertTrue( all( abs( lonR - 90 ) < 1e-7 ) );
assertTrue( all( abs( xyzI - [0;6500;0] ) < 1e-7 ) );

% Latitude must be 0
assertTrue( all( abs( lat ) < 1e-7 ) );

% Radius must be 6500 km
assertTrue( all( abs( r - 6500 ) < 1e-7 ) );

% Position in the orbital plane must be at 80 deg 
assertTrue( all( atand( pq(2)/pq(1) ) - 80 < 1e-7 ) );

% Predict orbit at half the orbital period
s.station_keeping = 1;  % Force station keeping
[ xyzI, xyzR, r, lat, lonI, lonR, pq ] = s.orbit_predictor( s.orbit_period/2, 1 );

assertTrue( all( abs( xyzI + [0;6500;0] ) < 1e-7 ) );
assertTrue( abs( lonI + 90 ) < 1e-7 );
assertTrue( abs( lonR - lonI + 360/(24*3600) * s.orbit_period/2  ) < 0.1 );      % Account for Erths rotation

% Calcuate UE Perspective
[ xyzU, visible, orientation ] = s.ue_perspective( [-98.89,-1],s.orbit_period/2 );

dD = 40000/360;  % km per degree
rE = 6371;       % Earth radius

assertTrue( abs( xyzU(3) - 6500 + rE ) < 20 )       % Height should roughly match
assertTrue( abs( xyzU(1) + 2*dD ) < 20 )            % X-positions should be 2 degrees west of the UE
assertTrue( abs( xyzU(2) - dD ) < 20 )              % Y-positions should be 1 degrees north of the UE

assertTrue( abs( orientation(1) +1 ) < 0.1 )        % Bank angle should be -1 degree
assertTrue( abs( orientation(2) -2 ) < 0.1 )        % Tilt angle should be 2 degree
assertTrue( abs( orientation(3) ) < 0.1 )           % Heading should be east

assertTrue( visible )
