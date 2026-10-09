function testBuilder_AzArr_Angle_Est
%%
% This test checks if the model generates the correct azimuth angles at the output

ant = qd_arrayant('omni');
ant.Fa(:) = 0;
ant.copy_element(1,2:36);
for n = 1 : 36
    ii = (n-1)*10 + (1:10) + 176;
    ii( ii > 361 ) = ii( ii > 361 )-361;
    %disp([ n ant.azimuth_grid( ii )*180/pi ]);
    ant.Fa( :,ii,n ) = 1;
end

b = qd_builder('Un');
b.tx_position = [0;0;25];
b.tx_array = qd_arrayant('omni');
b.rx_array = ant;
b.simpar.show_progress_bars = 0;
b.scenpar.SC_lambda = 0;
b.scenpar.PerClusterES_A = 0;
b.scenpar.PerClusterES_D = 0;
b.scenpar.PerClusterAS_A = 0;
b.scenpar.PerClusterAS_D = 0;
b.rx_track = qd_track([]);
b.rx_track.orientation = [0;0;pi/2];

max_dist = 300; % m
min_dist = 100;
no_user  = 10;
a = (2*rand(1,no_user)-1)*max_dist + 1j*(2*rand(1,no_user)-1)*max_dist;
while any( abs(a) < min_dist )
    ii = abs(a) < min_dist;
    a(ii) = (2*rand(1,sum(ii))-1)*max_dist + 1j*(2*rand(1,sum(ii))-1)*max_dist;
end

b.rx_positions = [ real(a) ; imag(a) ; ones(1,no_user)*1.5 ];
b.rx_positions(:,1) = [50;0;1.5];

gen_parameters(b);
c = get_channels(b);

AoAI = b.AoA*180/pi;
AoAO = zeros( no_user,b.NumClusters );
for n = 1:no_user
    for m = 1:b.NumClusters
        [~,ii] = max( abs( c(1,n).coeff(:,1,m) )  );
        ii = (ii-1)*10+90; % Angle estimate
        AoAO(n,m) = angle(exp( 1j*ii * pi/180 ))*180/pi;
    end
end
dA = angle(exp(1j*(AoAO - AoAI)*pi/180))*180/pi;

assertTrue( all( abs(dA(:))<6 ) );