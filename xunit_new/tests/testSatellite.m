function testSatellite 

s = qd_satellite([]); 
s = qd_satellite('gso',3,10);

assertTrue( all( abs( s(1).true_anomaly - [10 130 -110] ) < 1e-16) );

s = qd_satellite('walker-delta', qd_satellite.R_e + 100, 60, 3, 1  );

% Write TLE file
fid = fopen( 'astra1L.tle', 'w');
fprintf(fid,'ASTRA 1L                \n');
fprintf(fid,'1 31306U 07016A   20239.83454660  .00000100  00000-0  00000-0 0  9990\n');
fprintf(fid,'2 31306   0.0733 297.1425 0001796 267.1455  90.8519  1.00270503 23770\n');
fclose(fid);

s = qd_satellite.read_tle('astra1L.tle');

[ xyzI, xyzR, r, lat, lonI, lonR, pq, Omega, omega, v ] = s.orbit_predictor( 'utc+02' );
[ xyzI, xyzR, r, lat, lonI, lonR, pq, Omega, omega, v ] = s.orbit_predictor( 0 );
assertTrue( abs( lat ) < 0.1 );                 % GEO Satellite
assertTrue( abs( lonR-19.2 ) < 0.1 );           % Orbital slot of astra1L is 19.2 East
assertTrue( abs( r-qd_satellite.R_geo ) < 10 ); % GEO orbit height

[ xyzU, visible, orientation ] = s.ue_perspective([],'utc-03');
t = s.init_tracks([],'utc-06');
t = s.init_tracks([],200);
az = 90-atan2( t.initial_position(2), t.initial_position(1) )*180/pi;
el = atan2( t.initial_position(3), sqrt(sum(t.initial_position(1:2).^2)) )*180/pi;

assertTrue( abs(az-172.6 ) < 0.1 );             % Azimuth angle position from Berlin
assertTrue( abs(el-29.7 ) < 0.1 );              % Elevation angle from Berlin

delete('astra1L.tle')

s.visualize_lonlat('utc+02');
s.visualize_earth(0);
s.visualize_orbit;
close all

