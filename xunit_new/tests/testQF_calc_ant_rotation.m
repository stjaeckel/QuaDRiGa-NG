function testQF_calc_ant_rotation
%%

% Dummy
[ R, phiL, thetaL, gamma ] = qf.calc_ant_rotation( [], [], [], [0 0] );
assertTrue(all( abs( R(:) - [ 1;0;0 ; 0;1;0 ; 0;0;1 ] ) < 1e-10 ));
assertEqual( size(phiL,2) , 2 );
assertTrue(all( abs( phiL ) < 1e-10 ));
assertTrue(all( abs( thetaL ) < 1e-10 ));
assertTrue(all( abs( gamma ) < 1e-10 ));

% Rotate around z
[ R, phiL, thetaL, gamma ] = qf.calc_ant_rotation( pi/2, [], [], [0 0] );
assertTrue(all( abs( R(:) - [ 0;1;0 ; -1;0;0 ; 0;0;1 ] ) < 1e-10 ));
assertEqual( size(phiL,2) , 2 );
assertTrue(all( abs( phiL+pi/2 ) < 1e-10 ));
assertTrue(all( abs( thetaL ) < 1e-10 ));
assertTrue(all( abs( gamma ) < 1e-10 ));

[ R, phiL, thetaL, gamma ] = qf.calc_ant_rotation( [pi/2,-pi/2], [], [], [0 0] );
assertEqual( size(R,3) , 2 );
assertTrue(all( abs( phiL+[pi/2 -pi/2] ) < 1e-10 ));
assertTrue(all( abs( thetaL ) < 1e-10 ));
assertTrue(all( abs( gamma ) < 1e-10 ));

% Rotate around y
[ R, phiL, thetaL, gamma ] = qf.calc_ant_rotation( [], pi/4, [], [0 0] );
assertTrue(all( abs( R(:) - [ 1/sqrt(2);0;-1/sqrt(2) ; 0;1;0 ; 1/sqrt(2);0;1/sqrt(2) ] ) < 1e-10 ));
assertEqual( size(phiL,2) , 2 );
assertTrue(all( abs( phiL ) < 1e-10 ));
assertTrue(all( abs( thetaL-pi/4 ) < 1e-10 ));
assertTrue(all( abs( gamma ) < 1e-10 ));

[ R, phiL, thetaL, gamma ] = qf.calc_ant_rotation( [], [pi/4 -pi/2], [], [0 0] );
assertTrue(all( abs( phiL ) < 1e-10 ));
assertTrue(all( abs( thetaL-[pi/4 -pi/2] ) < 1e-10 ));
assertTrue(all( abs( gamma ) < 1e-10 ));

% Rotate around x
[ R, phiL, thetaL, gamma ] = qf.calc_ant_rotation( [], [], pi/4, [0 0] );
assertTrue(all( abs( R(:) - [ 1;0;0 ; 0;1/sqrt(2);1/sqrt(2) ; 0;-1/sqrt(2);1/sqrt(2) ] ) < 1e-10 ));
assertTrue(all( abs( phiL ) < 1e-10 ));
assertTrue(all( abs( thetaL ) < 1e-10 ));
assertTrue(all( abs( gamma-pi/4 ) < 1e-10 ));       % Positive

[ R, phiL, thetaL, gamma ] = qf.calc_ant_rotation( [], [], [pi/4 -pi/2], [0 0] );
assertTrue(all( abs( phiL ) < 1e-10 ));
assertTrue(all( abs( thetaL ) < 1e-10 ));
assertTrue(all( abs( gamma-[pi/4 -pi/2] ) < 1e-10 ));

% Combined y-z rotation
% This is equal to "z", "-x" rotation
[ R, phiL, thetaL, gamma ] = qf.calc_ant_rotation( pi/2, pi/4, [], [0 0] );
assertTrue(all( abs( R(:) - [ 0;1/sqrt(2);-1/sqrt(2) ; -1;0;0 ; 0;1/sqrt(2);1/sqrt(2) ] ) < 1e-10 ));
assertTrue(all( abs( phiL+pi/2 ) < 1e-10 ));
assertTrue(all( abs( thetaL ) < 1e-10 ));
assertTrue(all( abs( gamma+pi/4 ) < 1e-10 ));       % Negative

