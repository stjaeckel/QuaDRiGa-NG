function testBuilder_add_sdc
%clear all

%nnn=85
%RandStream.setGlobalStream(RandStream('mt19937ar','seed',nnn));

% Create layout
l = qd_layout;
l.simpar.autocorrelation_function = 'Disable';
l.simpar.center_frequency(2) = 3.7e9;
l.simpar.show_progress_bars = 0;
l.tx_position(:,2) = [100;0;25];
l.no_rx = 5;
l.randomize_rx_positions(200,1.5,1.5,0,[],50);
l.set_scenario('3GPP_38.901_UMi',[],1,0,0,0,0);
l.set_scenario('3GPP_38.901_UMi_LOS_GR',[],2,0,0,0,0);
%l.set_scenario('Freespace');

%%
% Initialize builders
b = init_builder( l );

% Split into one builder per RX
bs = split_rx( b );

% Initialize all parameters
gen_parameters( b );

% Absolute positioning, Power relative to freespace,
add_sdc( b, [-10 10 100 ;10 10 0 ;10 10 10],-3,[],[],[],0.1,[5 1 25] );

sic = size( b );
for i_cb = 1 : prod( sic )
    [ i1,i2 ] = qf.qind2sub( sic, i_cb );
    
    if b(i1,i2).no_rx_positions > 0
        % Check if the GR is present
        if check_los( b(i1,i2) ) == 0 || check_los( b(i1,i2) ) == 1
            iL = 0;
        else
            iL = 1;
        end
        
        % Check the number of subpaths
        assertEqual( b(i1,i2).NumSubPaths( (2:4)+iL ), [5 1 25] );
        
        % All power values must sum up to 1
        assertTrue( all( abs( reshape( sum( b(i1,i2).pow,2 ),[],1 ) - 1 ) < 1e-13 ) );
        
        % Accuray limit of FBS/LBS placement
        assertTrue( all( reshape(abs(b(i1,i2).fbs_pos(:,7+iL,:)-10),1,[]) < 0.2 ) );   
        assertTrue( all( reshape(abs(b(i1,i2).lbs_pos(:,7+iL,:)-10),1,[]) < 0.2 ) );
        
        % Total path length
        d = sqrt( sum( abs( [-10;10;10]*ones(1,b(i1,i2).no_rx_positions) - b(i1,i2).tx_position  ).^2,1 )) + ...
            sqrt( sum( abs( -[-10;10;10]*ones(1,b(i1,i2).no_rx_positions) + b(i1,i2).rx_positions  ).^2,1 ));
        
        % Delays relative to LOS delay
        dtr = sqrt(sum(abs(b(i1,i2).rx_positions - b(i1,i2).tx_position).^2,1));
        taus = (d-dtr)./qd_simulation_parameters.speed_of_light;
        assertTrue( all( b(i1,i2).taus(:,2+iL) - taus' < 1e-18 ) );
        
        % Freespace gain relative to total path length
        g = [1;1]*20*log10(d) + 20*log10( b(i1,i2).simpar(1,1).center_frequency'/1e9 )*ones(1,b(i1,i2).no_rx_positions) + 32.45;
        assertTrue( all(all( abs( 10*log10(reshape( b(i1,i2).gain(:,2+iL,:),[],2)') + 3 + g ) < 1e-12 )) );
    end
end

gen_parameters(bs,0);

% Add SDC to empty builder, power relative to freespce
add_sdc( bs(1,1), [10;10;10],0,[],[],[],0 );
assertEqual( bs(1,1).NumClusters, 2 )
assertEqual( bs(1,1).NumSubPaths, [1,1] )
assertTrue( all(  abs(bs(1,1).pow(:) - [0;1;0;1]) < 1e-13 ) );      % Zero-power-LOS
d3d = sqrt(sum(( bs(1,1).rx_positions - bs(1,1).tx_position ).^2));
path_lenght = (d3d / qd_simulation_parameters.speed_of_light + bs(1,1).taus(2) )*...
    qd_simulation_parameters.speed_of_light;
path_loss = 20*log10(path_lenght) + 32.45 + 20*log10( bs(1,1).simpar(1,1).center_frequency/1e9);
assertTrue( all(abs( 10*log10( squeeze(bs(1,1).gain(1,2,:)) ) + path_loss(:) ) < 1e-12) )

% Power relative to PL
add_sdc( bs(1,2), [10;10;10],0,[],'pathloss',[],0 );
d3d = sqrt(sum(( bs(1,2).rx_positions - bs(1,2).tx_position ).^2));
path_lenght = (d3d / qd_simulation_parameters.speed_of_light + bs(1,2).taus(2) )*...
    qd_simulation_parameters.speed_of_light;
hBS = bs(1,2).tx_position(3);
hUE = bs(1,2).rx_positions(3);
d2d = sqrt( path_lenght.^2 - ( hBS-hUE ).^2 );
path_loss = bs(1,2).get_pl( [ d2d; 0; hUE ], [], [ 0;0; hBS ] );  % PG
assertTrue( all(abs( 10*log10( squeeze(bs(1,2).gain(1,2,:)) ) + path_loss(:) ) < 1e-12) )

% Absolute power
add_sdc( bs(1,3), [10;10;10],-119,[],'absolute',[],0.1 );
assertEqual( bs(1,3).NumSubPaths, [1,20] )
assertTrue( all(abs( 10*log10( squeeze(bs(1,3).gain(1,2,:)) ) + 119 ) < 1e-12) )

% Position relative to RX
add_sdc( bs(1,4), [10;10;10],0,'rx_abs',[],[],0 );
assertTrue( all( abs( bs(1,4).fbs_pos(:,2,1,1) - bs(1,4).rx_positions - [10;10;10] ) < 1e-8 ) );

% Position relative to RX heading direction
add_sdc( bs(1,5), [10;0;0],0,'rx_heading',[],[],0 );
alpha = bs(1,5).rx_track.orientation(3,1);
assertTrue( all( abs( bs(1,5).fbs_pos(:,2,1,1) - bs(1,5).rx_positions - [cos(alpha)*10;sin(alpha)*10;0] ) < 1e-8 ) );

% Position relative to RX orientation
bs(1,6).rx_track(1,1).orientation = [ rand ; pi/4 ; pi/2  ];
add_sdc( bs(1,6), [],0,'rx_full',[],[],0,1,10,0,0 );
assertTrue( all( abs( bs(1,6).fbs_pos(:,2,1,1) - bs(1,6).rx_positions - [0;10;-10]/sqrt(2) ) < 1e-8 ) );

% Position relative to TX position
add_sdc( bs(1,1), [],0,'tx_abs',[],[],0,1,10,0,pi/2 );
assertEqual( bs(1,1).NumClusters, 3 );
assertEqual( bs(1,1).NumSubPaths, [1,1,1] );
assertTrue( all( abs( bs(1,1).fbs_pos(:,2,1,1) - bs(1,1).tx_position - [0;0;10] ) < 1e-8 ) );

% Position relative to TX heading direction
bs(1,7).tx_track(1,1).orientation = [ pi/4 ; 0 ; -pi/2  ];
add_sdc( bs(1,7), [0;10;0],0,'tx_heading',[],[],0 );
assertTrue( all( abs( bs(1,7).fbs_pos(:,2,1,end) - bs(1,7).tx_position - [10;0;0] ) < 1e-8 ) );

% Position relative to TX orientation
add_sdc( bs(1,7), [],0,'tx_full',[],[],0,1,10,pi/2,0 );
assertTrue( all( abs( bs(1,7).fbs_pos(:,2,1,end) - bs(1,7).tx_position - [10;0;10]/sqrt(2) ) < 1e-8 ) );

% Position relative to LOS path
add_sdc( bs(1,8), [10;5;0],0,'los_rx',[],[],0 );
assertTrue( abs( sum(( bs(1,8).fbs_pos(:,2,1,end) - bs(1,8).rx_positions ).^2) - 125 ) < 1e-8  );
d3d = sqrt(sum(( bs(1,8).rx_positions - bs(1,8).tx_position ).^2));
dpath = sqrt( sum( ( bs(1,8).fbs_pos(:,2,1,end) - bs(1,8).tx_position ).^2 ) ) + sqrt(125);
assertTrue( abs( dpath - d3d - bs(1,8).taus(2) * qd_simulation_parameters.speed_of_light  ) < 1e-8 )
add_sdc( bs(1,8), [],0,'los_tx',[],[],0,[],sqrt(200),0,pi/4 );
assertTrue( abs( sum(( bs(1,8).fbs_pos(:,2,1,end) - bs(1,8).tx_position ).^2) - 200 ) < 1e-8  );
dpath = sqrt( sum( ( bs(1,8).fbs_pos(:,2,1,end) - bs(1,8).rx_positions ).^2 ) ) + sqrt(200);

% Test XPR
gen_parameters(bs,0)
ang = bs(1,9).get_angles*pi/180;
bs(1,9).tx_track(1,1).orientation = [0;0;ang(1)];    % Pint BS towards RX
bs(1,9).rx_track(1,1).orientation = [0;0;ang(2)];    % Pint RX towards BS

add_sdc( bs(1,9), [10 10 ;5 -5 ;0 0],0,'los_rx',[],[10,16],0.1,[10 25] );
add_paths( bs(1,9), 'Freespace' );

[ xprL, xprC, cprL, cprC ] = qf.calc_xpr( bs(1,9).xprmat(:,1:end,1,1) );

xL = qd_builder.call_private_fcn( 'clst_avg', xprL, bs(1,9).NumSubPaths );
% xC = qd_builder.call_private_fcn( 'clst_avg', xprC, bs(1,9).NumSubPaths );
% cL = qd_builder.call_private_fcn( 'clst_avg', cprL, bs(1,9).NumSubPaths );
% cC = qd_builder.call_private_fcn( 'clst_avg', cprC, bs(1,9).NumSubPaths );

assertTrue( all( abs( xL(2:end) - 10.^(0.1*[10,16]) ) < 1e-8 ) )

bx = split_multi_freq( bs(1,9) );
bx = bx(1,1);

bx.tx_array = qd_arrayant('xpol');
bx.rx_array = qd_arrayant('xpol');

c = get_channels( bx );

xpr_L = qf.calc_xpr( reshape(c.coeff(:,:,2),[],1) );
assertTrue( abs( xpr_L - 10 ) < 0.2 )

xpr_L = qf.calc_xpr( reshape(c.coeff(:,:,3),[],1) );
assertTrue( abs( xpr_L - 40 ) < 0.5 )
% 
% bx.tx_array = qd_arrayant('lhcp-rhcp-dipole');
% bx.rx_array = qd_arrayant('lhcp-rhcp-dipole');
% 
% c = get_channels( bx );
% 
% xpr_C = qf.calc_xpr( reshape(c.coeff(:,:,3),[],1) );
% assertTrue( abs( xpr_C - 10 ) < 2 )


