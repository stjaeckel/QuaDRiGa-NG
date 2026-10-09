function testBuilder_get_channels
%%

%Disable some warnings
warning('off','QuaDRiGa:qd_builder:check_dual_mobility:no_rx_antenna');
warning('off','QuaDRiGa:qd_builder:check_dual_mobility:no_tx_antenna');

% Setup basic test scenario
b = qd_builder('3GPP_38.901_UMi_LOS');
b.rx_array = qd_arrayant('ula4');
b.tx_array = qd_arrayant('ula8');
b.name = 'Scen_Tx1';
b.simpar.show_progress_bars = false;
b.tx_position = [0;0;25];
b.rx_positions = [10,0,1.5;50,10,1.5]';
b.rx_track = qd_track.generate('linear',1,0);
b.rx_track.initial_position = b.rx_positions(:,1);
b.rx_track.name = 'Rx1';
b.rx_track(1,2) = qd_track.generate('circular',1);  % Closed track
b.rx_track(1,2).initial_position = b.rx_positions(:,2);
b.rx_track(1,2).interpolate_positions(5);
b.rx_track(1,2).name = 'Rx2';

set_rand_state(1);
gen_parameters(b); % Parameters are independent of the simulation options
rsA = get_rand_state; % Rand-state should not change at any time in the future

% Test differnt options
for n = 1 : 3
    switch n % Set options
        case 3
            b.simpar.use_3GPP_baseline = 1;
            gen_parameters(b,4); 
    end
    
    c = b.get_channels;
    rsB = get_rand_state;
    if n<3
        assertTrue( all( rsA-rsB == 0 ) );
    end
    
    assertEqual( numel(c) , 2 )
    assertEqual( c(1,1).name , 'Scen_Tx1_Rx1' )
    assertEqual( c(1,1).no_rxant , 4 );
    assertEqual( c(1,1).no_txant , 8 );
    assertEqual( c(1,1).no_snap , 2 );
    assertEqual( c(1,1).tx_position , b.tx_position(:,1) );
    assertEqual( c(1,1).rx_position , b.rx_track(1,1).positions_abs );
    assertFalse( any(isnan( c(1,1).coeff(:) )) );
    assertFalse( any(isnan( c(1,1).delay(:) )) );
    
    assertEqual( c(1,2).name , 'Scen_Tx1_Rx2' )
    assertEqual( c(1,2).rx_position , b.rx_track(1,2).positions_abs );
    assertFalse( any(isnan( c(1,2).coeff(:) )) );
    assertFalse( any(isnan( c(1,2).delay(:) )) );
    
    % Test features of different options
    switch n % Set options
        case 1
            assertEqual( c(1,1).individual_delays , true );     % Spherical waves
            assertEqual( c(1,2).coeff(:,:,:,1) , c(1,2).coeff(:,:,:,end) ); % Closed circutar track
            A = c(1,1).coeff(:,:,1,1); % Diagonal structure
            for aa = -6:6
                B = diag(A,n);
                assertTrue( all( abs(B-B(1)) < 1e-9 ) );
            end
        case 2
            assertEqual( c(1,1).coeff(:) , d(1,1).coeff(:) );   % Results in case 1 and 2 must be identical
            assertEqual( c(1,1).delay(:) , d(1,1).delay(:) );   % Results in case 1 and 2 must be identical
        case 3
            assertEqual( c(1,1).individual_delays , false );    % Planar waves
            
            assertTrue( c(1,1).no_path < d(1,1).no_path);       % Different subpath splitting
            %assertTrue(  sum( abs( c(1,1).coeff(:) - d(1,1).coeff(:)  ) ~= 0 ) > 500  ); % Different results
            
            % In low precision, the elements of the matrix are identical
            % Broadside of the arrays are aligned
            A = c(1,1).coeff(:,:,1,1);
            assertTrue( all( abs(A(:)-A(1,1)) < 1e-9 ) );
            A = c(1,1).coeff(:,:,1,2);
            assertTrue( all( abs(A(:)-A(1,1)) < 1e-9 ) );
    end
    
    d = c; % Copy handle for comaprisons in next step
end
