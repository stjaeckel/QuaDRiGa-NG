function testBuilder_LinearMovement
%% The Rx is moving away from the Tx
% The phase is tracked

b = qd_builder('LOSonly');
b.tx_position = [0;0;0];
b.rx_array = qd_arrayant('omni');
b.tx_array = qd_arrayant('omni');
b.name = 'Tx1';
b.simpar.show_progress_bars = false;
b.simpar.use_absolute_delays = true;
b.simpar.sample_density = 8;
b.rx_track = qd_track('linear', b.simpar.wavelength ,0);
b.rx_track(1,1).initial_position = [30;0;0];
b.rx_track(1,1).interpolate_positions( b.simpar.samples_per_meter );

for n = 1 : 4
    switch n
        case 2
            b.rx_track = qd_track('linear', b.simpar.wavelength ,pi/2);
            b.rx_track.initial_position = [0;30;0];
            b.rx_track(1,1).interpolate_positions( b.simpar.samples_per_meter );
        case 3
            b.rx_track = qd_track('linear', b.simpar.wavelength ,3*pi/4);
            b.rx_track.initial_position = [-30;30;0];
            b.rx_track(1,1).interpolate_positions( b.simpar.samples_per_meter );
        case 4
            b.rx_track = qd_track('linear', 0.9999999*b.simpar.wavelength ,-3*pi/4);
            b.rx_track.initial_position = [-30;-30;0];
            b.rx_track(1,1).interpolate_positions( b.simpar.samples_per_meter );
    end
    b.rx_positions = [];
    check_dual_mobility( b );
    gen_parameters(b);
    
    c = b.get_channels;
    
    ang = unwrap(angle( squeeze(c.coeff) ));
    ang = ang - ang(1);
    
    assertTrue(  all( abs(  ang(5)  + pi/2  ) < 0.01 ) );
    assertTrue(  all( abs( ang(9) +pi ) < 0.01 ) );
    assertTrue(  all( abs( ang(13) +3*pi/2 ) < 0.01 ) );
    assertTrue(  all( abs( ang(17) +2*pi ) < 0.01 ) );
end