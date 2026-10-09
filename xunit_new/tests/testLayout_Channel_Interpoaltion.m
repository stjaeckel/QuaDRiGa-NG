function testLayout_Channel_Interpoaltion
%%
%set_rand_state(1);

a = qd_arrayant;
a.copy_element(1,2);
a.element_position(1,:) = [-0.5,0.5];
a.combine_pattern;

l = qd_layout;
l.simpar.show_progress_bars = 0;
l.simpar.sample_density = 1.1;

l.no_rx = 1;
l.randomize_rx_positions( 10, 1.5, 1.5, 0.2 );
l.set_scenario( 'WINNER_UMa_C2_LOS' );
l.rx_array = a;

l.no_tx = 1;
l.tx_position(3,:) = 25;
l.tx_array = a;

[c,cb] = l.get_channels;

set_speed(l.rx_track, 3/3.6 );

for o = 1:2
    dist = interpolate( l.rx_track, 'time', 10e-3 );
    if o == 2
        dist = reshape( dist, 1,1,[] );
        dist = dist([1 1],[1 1],:);
    end

    for n = 0:1
        cb(1,1).simpar(1,1).use_3GPP_baseline = n;
        cb(1,1).gen_parameters;
        c = cb(1,1).get_channels;

        for m = 1:2
            switch m
                case 1
                    ci = interpolate(c, dist,'linear' );
                case 2
                    ci = interpolate(c, dist,'cubic' );
            end

            % There are more snapshots in the interpolated version
            assertTrue( ci.no_snap > c.no_snap )

            % Tx-pos
            assertTrue( all(  abs( ci.tx_position - l.tx_position ) < 1e-13 ) );

            % First Rx-Pos
            assertTrue( all(  abs( ci.rx_position(:,1) - l.rx_track(1,1).initial_position ) < 1e-13 ) );
            assertTrue( all(  abs( ci.rx_position(:,1) - c.rx_position(:,1) ) < 1e-13 ) );

            % First coeff is the same
            X = ci.coeff(:,:,:,1) - c.coeff(:,:,:,1);
            assertTrue( all( abs(X(:)) < 1e-13 ) );

            % First delay is the same
            if ci.individual_delays
                X = ci.delay(:,:,:,1) - c.delay(:,:,:,1);
                assertTrue( all( abs(X(:)) < 1e-13 ) );
            else
                X = ci.delay(:,1) - c.delay(:,1);
                assertTrue( all( abs(X(:)) < 1e-13 ) );
            end

        end
    end
end
