function testBuilder_genFbsLbs
% This tests the FBS and LBS generation

N = 5;
L = 10;

tx_pos  = rand( 3,N )*100;
rx_pos  = rand( 3,N )*100;
taus    = rand( N,L )*99e-9 + 1e-9;
AoD     = rand( N,L )*2*pi - pi;
AoA     = rand( N,L )*2*pi - pi;
EoD     = 0.5*(rand( N,L )*pi - pi/2);
EoA     = 0.5*(rand( N,L )*pi - pi/2);
NumSubPaths = ones(1,L);
SubPathCPL = [];
PerClusterAS = [];

% Calculate angles between BS and MT
d_2d = hypot( tx_pos(1,:) - rx_pos(1,:), tx_pos(2,:) - rx_pos(2,:) );
d_2d( d_2d<1e-5 ) = 1e-5;
angles = zeros( 5,N );
angles(1,:) = atan2( rx_pos(2,:) - tx_pos(2,:) , rx_pos(1,:) - tx_pos(1,:) );   % Azimuth at BS
angles(2,:) = pi + angles(1,:);                                                 % Azimuth at MT
angles(3,:) = atan2( ( rx_pos(3,:) - tx_pos(3,:) ), d_2d );                     % Elevation at BS
angles(4,:) = -angles(3,:);                                                     % Elevation at MT
angles(5,:) = -atan2( ( rx_pos(3,:) + tx_pos(3,:) ), d_2d );                    % Ground Reflection Elevation at BS and MT
angles = angles.';


%% Pure NLOS, single-frequency, no sub-paths, bounce2

b = qd_builder;
b.tx_position = tx_pos;
b.rx_positions = rx_pos;
b.taus = taus;
b.AoD = AoD;
b.AoA = AoA;
b.EoD = EoD;
b.EoA = EoA;

b.gen_fbs_lbs;
b.gen_ssf_from_scatterers;

fbs_pos = b.fbs_pos;
lbs_pos = b.lbs_pos;
AoD_c = b.AoD;
AoA_c = b.AoA;
EoD_c = b.EoD;
EoA_c = b.EoA;

% Check size of outputs
assertEqual( size(fbs_pos), [3,L,N] )
assertEqual( size(lbs_pos), [3,L,N] )
assertEqual( size(AoD_c),   [N,L] )
assertEqual( size(AoA_c),   [N,L] )
assertEqual( size(EoD_c),   [N,L] )
assertEqual( size(EoA_c),   [N,L] )

% Check the path lengths
oL  = ones(1,L);
T   = permute( tx_pos , [1,3,2] );
R   = permute( rx_pos , [1,3,2] );
r   = -T(:,oL,:) + R(:,oL,:);
b   = -T(:,oL,:) + fbs_pos;
c   = -fbs_pos + lbs_pos;
a   = -lbs_pos + R(:,oL,:);
dx  = sqrt( sum(b.^2,1) ) + sqrt( sum(c.^2,1) ) + sqrt( sum(a.^2,1) ) - sqrt( sum(r.^2,1) );
dx  = permute( dx,[3,2,1] );
di  = taus*qd_simulation_parameters.speed_of_light;

% In "bounce2", all path lengths must match the given taus
assertTrue( all( abs( dx(:) - di(:) ) < 1e-8 ) );

% Arrival angles shoule be identical
assertTrue( all(abs(AoA_c(:) - AoA(:)) < 1e-8) )
assertTrue( all(abs(EoA_c(:) - EoA(:)) < 1e-8) )

% Some departure angles should match, but not all (due to single-bounce)
tmp = sum(abs(AoD_c(:) - AoD(:)) < 1e-8);
assertTrue( tmp > N*L/10 )
assertTrue( tmp < N*L )
assertEqual( sum(abs(EoD_c(:) - EoD(:)) < 1e-8), tmp  )

%% Pure NLOS, single-frequency, with sub-paths, bounce2
NumSubPaths = randi(19,1,L)+1;
ML = sum( NumSubPaths );
SubPathCPL = rand(4,ML);

b = qd_builder('Null');
b.scenpar.PerClusterAS_A = 5;
b.scenpar.PerClusterAS_D = 5;
b.scenpar.PerClusterES_A = 5;
b.scenpar.PerClusterES_D = 5;
b.tx_position = tx_pos;
b.rx_positions = rx_pos;
b.NumSubPaths = NumSubPaths;
b.subpath_coupling = SubPathCPL;
b.taus = taus;
b.AoD = AoD;
b.AoA = AoA;
b.EoD = EoD;
b.EoA = EoA;

b.gen_fbs_lbs;
b.gen_ssf_from_scatterers;

fbs_pos = b.fbs_pos;
lbs_pos = b.lbs_pos;
AoD_c = b.AoD;
AoA_c = b.AoA;
EoD_c = b.EoD;
EoA_c = b.EoA;

% Check size of outputs
assertEqual( size(fbs_pos), [3,ML,N] )
assertEqual( size(lbs_pos), [3,ML,N] )
assertEqual( size(AoD_c),   [N,L] )
assertEqual( size(AoA_c),   [N,L] )
assertEqual( size(EoD_c),   [N,L] )
assertEqual( size(EoA_c),   [N,L] )

% Arrival angles shoule be similar
assertTrue( all(abs(  angle(exp(1j*( AoA_c(:) - AoA(:) )))  )  < 0.5) )
assertTrue( all(abs(EoA_c(:) - EoA(:)) < 0.1) )

% Some departure angles should match, but not all (due to single-bounce)
tmp = sum(abs(AoD_c(:) - AoD(:)) < 0.5);
assertTrue( tmp > N*L/10 )
assertTrue( tmp < N*L )

tmp = sum(abs(EoD_c(:) - EoD(:)) < 0.5);
assertTrue( tmp > N*L/10 )
assertTrue( tmp < N*L )

%% LOS, single-frequency, with sub-paths, bounce2
NumSubPaths(1) = 1;
ML = sum( NumSubPaths );
SubPathCPL = rand(4,ML);

AoD(:,1) = angles(:,1);
AoA(:,1) = angles(:,2);
EoD(:,1) = angles(:,3);
EoA(:,1) = angles(:,4);
taus(:,1) = 0;

b = qd_builder('Null');
b.scenpar.PerClusterAS_A = 5;
b.scenpar.PerClusterAS_D = 5;
b.scenpar.PerClusterES_A = 5;
b.scenpar.PerClusterES_D = 5;
b.tx_position = tx_pos;
b.rx_positions = rx_pos;
b.NumSubPaths = NumSubPaths;
b.subpath_coupling = SubPathCPL;
b.taus = taus;
b.AoD = AoD;
b.AoA = AoA;
b.EoD = EoD;
b.EoA = EoA;

b.gen_fbs_lbs;
b.gen_ssf_from_scatterers;

fbs_pos = b.fbs_pos;
lbs_pos = b.lbs_pos;
AoD_c = b.AoD;
AoA_c = b.AoA;
EoD_c = b.EoD;
EoA_c = b.EoA;


% Check size of outputs
assertEqual( size(fbs_pos), [3,ML,N] )
assertEqual( size(lbs_pos), [3,ML,N] )
assertEqual( size(AoD_c),   [N,L] )
assertEqual( size(AoA_c),   [N,L] )
assertEqual( size(EoD_c),   [N,L] )
assertEqual( size(EoA_c),   [N,L] )

% Identical LBS and FBS positions
tmp = lbs_pos(:,1,:) - fbs_pos(:,1,:);
assertTrue( all(abs( tmp(:) ) < 1e-7) )

% Path lengths matche the LOS path lengths
T   = permute( tx_pos , [1,3,2] );
R   = permute( rx_pos , [1,3,2] );
r   = -T(:,1,:) + R(:,1,:);
b   = -T(:,1,:) + fbs_pos(:,1,:);
c   = -fbs_pos(:,1,:) + lbs_pos(:,1,:);
a   = -lbs_pos(:,1,:) + R(:,1,:);
dx  = sqrt( sum(b.^2,1) ) + sqrt( sum(c.^2,1) ) + sqrt( sum(a.^2,1) ) - sqrt( sum(r.^2,1) );
assertTrue( all( abs( dx(:) ) < 1e-7 ) )

%% LOS, multi-frequency, with sub-paths, bounce2
SubPathCPL = rand(4,ML,2);
SubPathCPL(:,:,3) = SubPathCPL(:,:,1);
SubPathCPL(:,:,4) = SubPathCPL(:,:,1);

% PerClusterAS = ones(4,4)*5;
% PerClusterAS(:,4)= 2;

b = qd_builder('Null');
b.simpar.center_frequency = [1:4]*1e9;
b.scenpar.PerClusterAS_A = 5;
b.scenpar.PerClusterAS_D = 5;
b.scenpar.PerClusterES_A = 5;
b.scenpar.PerClusterES_D = 5;
b.tx_position = tx_pos;
b.rx_positions = rx_pos;
b.NumSubPaths = NumSubPaths;
b.subpath_coupling = SubPathCPL;
b.taus = taus;
b.AoD = AoD;
b.AoA = AoA;
b.EoD = EoD;
b.EoA = EoA;

b.gen_fbs_lbs;
b.gen_ssf_from_scatterers;

fbs_pos = b.fbs_pos;
lbs_pos = b.lbs_pos;
AoD_c = b.AoD;
AoA_c = b.AoA;
EoD_c = b.EoD;
EoA_c = b.EoA;

% Check size of outputs
assertEqual( size(fbs_pos), [3,ML,N,4] )
assertEqual( size(lbs_pos), [3,ML,N,4] )
assertEqual( size(AoD_c),   [N,L] )
assertEqual( size(AoA_c),   [N,L] )
assertEqual( size(EoD_c),   [N,L] )
assertEqual( size(EoA_c),   [N,L] )

% LOS-scatterer should be the same
tmp = fbs_pos(:,1,:,1) - fbs_pos(:,1,:,2);
assertTrue( all( abs(tmp(:)) < 1e-12 ) );

% Most NLOS scatterer should be different
tmp = fbs_pos(:,2:end,:,1) - fbs_pos(:,2:end,:,2);
assertTrue( sum( abs(tmp(:)) < 1e-12  ) < N*L*3/2 )

% Identical clusters should be created if SubPathCPL and PerClusterAS are the same 
tmp = fbs_pos(:,:,:,1) - fbs_pos(:,:,:,3);
assertTrue( all( abs(tmp(:)) < 1e-12 ) );


%% LOS+GR, multi-frequency, with sub-paths, bounce2
NumSubPaths(2) = 1;
ML = sum( NumSubPaths );
SubPathCPL = rand(4,ML,2);

d_3d = sqrt( sum((rx_pos - tx_pos).^2,1) );
d_gf = sqrt( sum(([rx_pos(1:2,:);-rx_pos(3,:)] - tx_pos).^2,1) );
tau_gr = ( d_gf-d_3d ).' / qd_simulation_parameters.speed_of_light;

AoD(:,2) = AoD(:,1);
AoA(:,2) = AoA(:,1);
EoD(:,2) = angles(:,5);
EoA(:,2) = angles(:,5);
taus(:,2) = tau_gr;

b = qd_builder('Null');
b.simpar.center_frequency = [1,2]*1e9;
b.scenpar.PerClusterAS_A = 5;
b.scenpar.PerClusterAS_D = 5;
b.scenpar.PerClusterES_A = 5;
b.scenpar.PerClusterES_D = 5;
b.tx_position = tx_pos;
b.rx_positions = rx_pos;
b.NumSubPaths = NumSubPaths;
b.subpath_coupling = SubPathCPL;
b.taus = taus;
b.AoD = AoD;
b.AoA = AoA;
b.EoD = EoD;
b.EoA = EoA;

b.gen_fbs_lbs;
b.gen_ssf_from_scatterers;

fbs_pos = b.fbs_pos;
lbs_pos = b.lbs_pos;
AoD_c = b.AoD;
AoA_c = b.AoA;
EoD_c = b.EoD;
EoA_c = b.EoA;

% GR FBS and LBS must be identical
tmp = fbs_pos(:,2,:,:) - lbs_pos(:,2,:,:);
assertTrue( all( abs(tmp(:)) < 1e-12 ) );

% GR FBS and must be identical for all frequencies
tmp = fbs_pos(:,2,:,1) - fbs_pos(:,2,:,2);
assertTrue( all( abs(tmp(:)) < 1e-12 ) );

% Scatterer must be on the ground
tmp = fbs_pos(3,2,:,:);
assertTrue( all( abs(tmp(:)) < 1e-8 ) );


