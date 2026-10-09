function testLayout_Closed_Chan_Generation
%%

set_rand_state( 1 );

t = qd_track.generate('circular',2);    % Closed track, 4 m length
t.initial_position = [50,0,1.5]';
t.set_speed(2);
t.scenario = 'unit_test_conf';

l = qd_layout;
l.simpar.show_progress_bars = 0;
l.simpar.center_frequency = 500e6;
l.simpar.use_absolute_delays = 1;
l.simpar.sample_density = 2.1;

for n = 1:6
    switch n
        case 1 % No segments, Single Antenna, Linear Interpolation
            l.rx_track = t.copy;
            l.tx_array = qd_arrayant('omni');
            l.rx_array = qd_arrayant('omni');
            l.update_rate = 10e-3;
            c = l.get_channels;
            
        case 2 % No segments, Multi-Antenna, Linear Interpolation
            l.rx_track = t.copy;
            l.rx_array = qd_arrayant('ula2');
            l.tx_array = qd_arrayant('ula4');
            l.update_rate = 10e-3;
            c = l.get_channels;
            
        case 3 % No segments, Single Antenna, Cubic Interpolation
            l.rx_track = t.copy;
            l.tx_array = qd_arrayant('omni');
            l.rx_array = qd_arrayant('omni');
            l.update_rate = 10e-3;
            c = l.get_channels([],[],'cubic');
            
        case 4 % No segments, Multi-Antenna, Cubic Interpolation
            l.rx_track = t.copy;
            l.rx_array = qd_arrayant('ula2');
            l.tx_array = qd_arrayant('ula4');
            l.update_rate = 10e-3;
            c = l.get_channels([],[],'cubic');
            
        case 5 % Segments, Single Antenna, Linear Interpolation
            l.rx_track = t.copy;
            l.rx_track.segment_index = [1 50];
            l.tx_array = qd_arrayant('omni');
            l.rx_array = qd_arrayant('omni');
            l.update_rate = 10e-3;
            c = l.get_channels;
            
        case 6 % Tx mobility, segments
            l.rx_track = t.copy;
            l.tx_track = t.copy;
            l.tx_position(1,1) = 55;
            l.rx_track.segment_index = [1 50];
            l.update_rate = 10e-3;
            c = l.get_channels;
    end
    
    assertEqual( c.no_snap, 101 );
    assertTrue( sum( abs( c.rx_position(:,1) - c.rx_position(:,end) ).^2 ) < 1e-5 );
    assertTrue( sum( abs( c.tx_position(:,1) - c.tx_position(:,end) ).^2 ) < 1e-5 );
    if n < 5
        x = c.coeff(:,:,:,1) - c.coeff(:,:,:,end);
        assertTrue( all( abs( x(:) ) < 1e-11 ) );
        x = c.delay(:,:,:,1) - c.delay(:,:,:,end);
        assertTrue( all( abs( x(:) ) < 1e-11 ) );
        assertTrue( abs( c.par.pg(1) - c.par.pg(end) ) < 1e-5 );
    end
    assertTrue( ~isempty( c.par ) );
    assertEqual( numel( c.par.pg ) , c.no_snap );
    
end
