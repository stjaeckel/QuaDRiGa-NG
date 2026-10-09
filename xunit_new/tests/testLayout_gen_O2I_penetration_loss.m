function testLayout_gen_O2I_penetration_loss
%%

s = qd_simulation_parameters;
s.center_frequency = [ 2e9, 28e9 ];                 % Two frequencies
s.show_progress_bars = 0;

l = qd_layout( s );
l.no_tx = 2;                                        % 2 BS
l.tx_position(2,2) = 100;                           % 100 m ISD

l.no_rx = 4;                                        % 10 Rx
l.randomize_rx_positions( 500,1.5,50,0 );           % 500 m radius, random heights
l.set_scenario('3GPP_38.901_UMi_LOS_O2I');          % Set scenario

t = qd_track('street',50,0);                        % New track with MT mobility
t.name = l.rx_name{1};                              % Keep Rx name
t.initial_position = [0,0,1.5]';                    % Start point under BS1

%t.interpolate_positions(10);

% The track onla has O2I segments
t.set_scenario({'3GPP_38.901_UMi_LOS_O2I','3GPP_38.901_UMi_NLOS_O2I'},[0.5,0.5],[]);

l.rx_track(1,1) = t;                                   % Assign track

o2i_loss = l.gen_o2i_loss('3GPP_38.901');                      % Generate o2i penetration loss
o2i_loss = l.gen_o2i_loss('mmMAGIC');                          % Generate o2i penetration loss

% Non-O2I scenario mix of one outdoor and O2I scenario at the same position
l.rx_track(1,2).scenario= {'3GPP_38.901_UMi_LOS';'3GPP_38.901_UMi_LOS_O2I'};     
l.rx_track(1,2).par.o2i_d3din(1,:,:) = 0;
l.rx_track(1,2).par.o2i_loss(1,:,:)  = 0;

% One MT has all outdoor scenarios
l.rx_track(1,3).scenario= {'3GPP_38.901_UMi_LOS';'3GPP_38.901_UMi_NLOS'}; 
l.rx_track(1,3).par = [];

assertEqual( size( o2i_loss{1} ), [ l.no_tx, l.rx_track(1,1).no_segments, numel(l.simpar.center_frequency)] )
assertEqual( l.rx_track(1,1).par.o2i_loss, o2i_loss{1} );

% Loss at 28 GHz must be larger than at 2 GHz
assertTrue( all(all( o2i_loss{1}(:,:,2) > o2i_loss{1}(:,:,1) )) );

% MT2, MS1 is outdoor
assertTrue( l.rx_track(1,2).par.o2i_loss(1,:,1) == 0 )
assertTrue( l.rx_track(1,2).par.o2i_loss(1,:,2) == 0 )
assertTrue( l.rx_track(1,2).par.o2i_d3din(1,:) == 0 )

% MT3 has no indoor users
assertTrue( isempty( l.rx_track(1,3).par ));

% Test if the segments are correctly split
l.rx_track(1,1).split_segment;

% Split track into subtracks
subtrack = l.rx_track(1,1).get_subtrack;

assertEqual( subtrack(1,1).par.o2i_loss, l.rx_track(1,1).par.o2i_loss(:,1,:) );
assertEqual( subtrack(1,2).par.o2i_loss, l.rx_track(1,1).par.o2i_loss(:,1:2,:) );
assertEqual( subtrack(1,end).par.o2i_loss, l.rx_track(1,1).par.o2i_loss(:,end-1:end,:) );

assertEqual( subtrack(1,1).par.o2i_d3din, l.rx_track(1,1).par.o2i_d3din(:,1) );
assertEqual( subtrack(1,2).par.o2i_d3din, l.rx_track(1,1).par.o2i_d3din(:,1:2) );
assertEqual( subtrack(1,end).par.o2i_d3din, l.rx_track(1,1).par.o2i_d3din(:,end-1:end) );

% Initialize builder objects
b = l.init_builder;
for n = 1:size( b,1 )
    for m = 1:size(b,2)
        scenpar_tmp = b(n,m).scenpar;
        scenpar_tmp.SC_lambda = 0;           % Disable SC for faster processing
        scenpar_tmp.SF_sigma = 0;            % No SF
        scenpar_tmp.XPR_mu = 100;            % No XPR
        scenpar_tmp.XPR_sigma = 0;           
        scenpar_tmp.PerClusterAS_A = 0;      % No cluster AS and DS
        scenpar_tmp.PerClusterAS_D = 0;
        scenpar_tmp.PerClusterES_A = 0;
        scenpar_tmp.PerClusterES_D = 0;
        scenpar_tmp.PerClusterDS = 0;
        b(n,m).scenpar = scenpar_tmp;
        b(n,m).plpar = [];                      % No PL, only O2I loss
    end
end

% Check if the parameter structs were correclty split
assertEqual( b(1,1).rx_track(1,1).par.o2i_loss, l.rx_track(1,1).par.o2i_loss(1,1,:) );
assertEqual( b(1,2).rx_track(1,1).par.o2i_loss, l.rx_track(1,1).par.o2i_loss(2,1,:) );

gen_parameters( b );                % Generate SSF parameters
b = split_multi_freq( b );              % Split the 2 frequencies

% Check if the parameter structs were correclty split
assertEqual( b(1,1).rx_track(1,1).par.o2i_loss, l.rx_track(1,1).par.o2i_loss(1,1,1) );
assertEqual( b(1,2).rx_track(1,1).par.o2i_loss, l.rx_track(1,1).par.o2i_loss(2,1,1) );
assertEqual( b(1,3).rx_track(1,1).par.o2i_loss, l.rx_track(1,1).par.o2i_loss(1,1,2) );
assertEqual( b(1,4).rx_track(1,1).par.o2i_loss, l.rx_track(1,1).par.o2i_loss(2,1,2) );

assertEqual( b(1,1).rx_track(1,1).par.o2i_d3din, l.rx_track(1,1).par.o2i_d3din(1,1) );
assertEqual( b(1,2).rx_track(1,1).par.o2i_d3din, l.rx_track(1,1).par.o2i_d3din(2,1) );
assertEqual( b(1,3).rx_track(1,1).par.o2i_d3din, l.rx_track(1,1).par.o2i_d3din(1,1) );
assertEqual( b(1,4).rx_track(1,1).par.o2i_d3din, l.rx_track(1,1).par.o2i_d3din(2,1) );

c = get_channels( b );                  % Get channel coefficients

d = merge( c ,[] ,l.simpar.show_progress_bars );       % Merge coefficients
d = qf.reshapeo( d, [ l.no_rx, l.no_tx, numel( l.simpar.center_frequency ) ] );  % Order

% Calculate actual PG from the output channel coefficients
pg = 10*log10( squeeze( sum( abs( d(1,1,1).coeff ).^2, 3 ) ) );    % BS1, F1

assertTrue( abs( pg(1,1)   + l.rx_track(1,1).par.o2i_loss(1,1,1) )   < 1e-4 );  % BS1, MT1, F1, Start
assertTrue( abs( pg(end,1) + l.rx_track(1,1).par.o2i_loss(1,end,1) ) < 1e-4 );  % BS1, MT1, F1, End

pg(:,2) = 10*log10( squeeze( sum( abs( d(1,1,2).coeff ).^2, 3 ) ) );    % BS1, F2
assertTrue( abs( pg(1,2)   + l.rx_track(1,1).par.o2i_loss(1,1,2) )   < 1e-4 );  % BS1, MT1, F2, Start
assertTrue( abs( pg(end,2) + l.rx_track(1,1).par.o2i_loss(1,end,2) ) < 1e-4 );  % BS1, MT1, F2, End

pg(:,3) = 10*log10( squeeze( sum( abs( d(1,2,1).coeff ).^2, 3 ) ) );    % BS2, F1
assertTrue( abs( pg(1,3)   + l.rx_track(1,1).par.o2i_loss(2,1,1) )   < 1e-4 );  % BS2, MT1, F2, Start
assertTrue( abs( pg(end,3) + l.rx_track(1,1).par.o2i_loss(2,end,1) ) < 1e-4 );  % BS2, MT1, F2, End

pg(:,4) = 10*log10( squeeze( sum( abs( d(1,2,2).coeff ).^2, 3 ) ) );    % BS2, F2
assertTrue( abs( pg(1,4)   + l.rx_track(1,1).par.o2i_loss(2,1,2) )   < 1e-4 );  % BS2, MT1, F2, Start
assertTrue( abs( pg(end,4) + l.rx_track(1,1).par.o2i_loss(2,end,2) ) < 1e-4 );  % BS2, MT1, F2, End

% Get the target path gain along the trajectory
pg_opt = d(1,1,1).par.pg';
pg_opt(:,2) = d(1,1,2).par.pg';
pg_opt(:,3) = d(1,2,1).par.pg';
pg_opt(:,4) = d(1,2,2).par.pg';

if 0 % Debugging plot
    figure(1)
    plot( pg_opt );
    hold on
    plot( pg, '--' );
    hold off
end

% Compare actual and target PG - the difference should be 0
pg_diff = pg - pg_opt;
segment_index = l.rx_track(1,1).segment_index;

tmp = pg_diff( 1:floor( segment_index(2)/2 )-1 , : );
assertTrue( all( abs(tmp(:)) < 1e-4 ) );

tmp = pg_diff( segment_index(end):end , : );
assertTrue( all( abs(tmp(:)) < 1e-4 ) );

% MT2, BS1 is outdoor
pgX = sum( squeeze(abs( d(2,1,1).coeff ).^2) );
assertTrue( abs(pgX - 1) < 1e-6 )

% MT2, BS2 is indoor
pgX = sum( squeeze(abs( d(2,2,1).coeff ).^2) );
o2i_loss = 10^(-0.1*l.rx_track(1,2).par.o2i_loss(2,1,1));
assertTrue( abs(pgX - o2i_loss) < 1e-6 )

% MT3 outdoor
for n = 1:2
    for m = 1:2
        pgX = sum( squeeze(abs( d(3,n,m).coeff ).^2) );
        assertTrue( abs(pgX - 1) < 1e-6 )
    end
end

% MT4 indoor
for n = 1:2
    for m = 1:2
        pgX = sum( squeeze(abs( d(4,n,m).coeff ).^2) );
        o2i_loss = 10^(-0.1*l.rx_track(1,4).par.o2i_loss(n,1,m));
        assertTrue( abs(pgX - o2i_loss) < 1e-6 )
    end
end

assertTrue( ~isempty( l.rx_track(1,1).par.o2i_loss ) );
assertTrue( ~isempty( l.rx_track(1,1).par.o2i_d3din ) );

