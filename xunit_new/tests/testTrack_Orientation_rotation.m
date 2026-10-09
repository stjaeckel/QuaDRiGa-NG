function testTrack_Orientation_rotation
%%

t = qd_track('linear',0,0);
t.positions = [ 0 0.25 0.5 0.75 1 ; 0 0 0 0 0 ; 0 0 0 0 0];                 % Track
t.orientation = [ 0 90 180 270 0 ; 0 90 0 -90 0 ; 0 0 171 0 0]*pi/180;
t.interpolate_positions( 36 );
o = t.orientation*180/pi;

% Roll
ref = [ 0:10:170 , -170:10:0 ];
assertTrue(  all( abs( o(1,[1:18,20:end])-ref ) < 1e-13 ) ); % [deg]

o = o(:,[1:9,11:27,29:end]);    % Remove poles

% Pitch
ref = [ 0:10:80 , 80:-10:-80 , -80:10:0 ];
assertTrue(  all( abs( o(2,:)-ref ) < 1e-13 ) ); % [deg]

% Yaw
ref = [ zeros(1,9) , 19:19:170 , 171:-19:1 , zeros(1,9) ];
assertTrue(  all( abs( exp(1j*o(3,:)*pi/180) - exp(1j*ref*pi/180) ) < 1e-13 ) );       % [rad]
