function testTrack_interpolate
%%

% Test different interpolation algorithms
alg = {'linear','cubic'};
for ia = 1 : numel( alg )
    t = qd_track('linear',10,pi/2);
       [d,ti] = t.interpolate('distance',0.5,[],alg{ia});
    
    assertEqual( ti.no_snapshots, 21 );
    assertTrue( all( abs( ti.positions(2,:) - (0:0.5:10)  ) < 1e-6 ) );
    assertTrue( all( abs( d - (0:0.5:10)  ) < 1e-6 ) );
    assertTrue( all( abs( ti.orientation(3,:) - pi/2  ) < 1e-6 ) );
end

t = qd_track('linear',10,pi/2);                                         
t.name = 'bla';
[~,ti] = t.interpolate('distance',5,[],[],1);               % Update input track

assertTrue( qf.eqo(t,ti) );                                 % Same handle for output and input

t.segment_index = [1 2];                                    % Set segments
t.scenario = {'a','b'};                                     % Set scenarios
t.set_speed(2);                                             % Set speed

[d,ti] = t.interpolate('time',0.1);                         % Time-base interpolation 100 ms SR

assertFalse( qf.eqo(t,ti) );                                % Different handle for output and input
assertEqual( ti.no_snapshots, 51 );
assertTrue( all( abs( ti.positions(2,:) - (0:0.2:10)  ) < 1e-6 ) );
assertEqual( ti.segment_index,[1,26] );
assertEqual( ti.name,'bla' );

movement_profile = t.movement_profile;
t.movement_profile = [];
[d,ti] = t.interpolate('time',0.1,movement_profile,'cubic');      % Input movement profile

assertEqual( ti.movement_profile, movement_profile );       % Check if MP was assignet to new track
assertEqual( ti.no_snapshots, 51 );
assertTrue( all( abs( ti.positions(2,:) - (0:0.2:10)  ) < 1e-6 ) );
assertEqual( ti.segment_index,[1,26] );
assertEqual( ti.name,'bla' );

t = qd_track('circular',10,0);
t.set_speed(1.25);
[d,ti] = t.interpolate('time',1);

assertEqual( ti.no_snapshots, 9 );
assertTrue( ti.closed );

ang = angle( ti.positions(1,:) + 10/(2*pi) + 1j*ti.positions(2,:) );
assertTrue( all( abs( exp(1j*( 0:45:360 )*pi/180) - exp(1j*ang) ) < 1e-6 ))

t = qd_track('linear',0);
t.no_snapshots = 5;
t.orientation = [ 0,0,0,0,0 ; 0,0,0,0,0 ; 0,90,180,270,0 ]*pi/180;
t.movement_profile = [0 1 ; 1 5];

[d,ti] = t.interpolate('snapshot',1/16);

ang = ti.orientation(3,:);
assertTrue( all( abs( exp(1j*( 0:22.5:360 )*pi/180) - exp(1j*ang) ) < 1e-6 ))

%% Test for NaNs due to co-located points

t = qd_track('linear',10);
t.scenario = 'bla';
t.interpolate('distance',1,[],[],1);

t.positions(:,5) = t.positions(:,4);
t.positions(:,6) = t.positions(:,4);
t.positions(:,end+1) = t.positions(:,end);

t.interpolate('distance',0.45,[],[],1);

assertEqual( t.no_snapshots, 23 );

assertTrue(  all(~isnan( t.positions(:) )) );

t.positions(:,end+1) = t.positions(:,end);

t.segment_index(2) = t.no_snapshots;

t.interpolate('distance',0.45,[],[],1);

assertEqual( t.no_snapshots, 23 );

assertEqual( t.segment_index, [1 23] );






