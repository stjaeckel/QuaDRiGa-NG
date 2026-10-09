function testLayout_GetChannels
%%
l = qd_layout;
l.simpar.show_progress_bars = 0;
l.simpar.center_frequency = 500e6;
l.simpar.use_absolute_delays = 1;
l.no_tx = 2;
l.no_rx = 2;
l.randomize_rx_positions( 100, 1.5, 1.5, 1.999 );
l.rx_track(1,1).segment_index = [1 8 16 ];
l.rx_track(1,1).set_speed(1);

l.rx_track(1,2) = qd_track.generate('circular',4);  % Closed track
l.rx_track(1,2).initial_position = [50,0,1.5]';
l.rx_track(1,2).interpolate_positions(l.simpar.samples_per_meter);
l.rx_track(1,2).name = 'Rx2';
l.rx_track(1,2).set_speed(2);

l.set_scenario('Freespace');
l.rx_track(1,1).scenario{2,2} = 'LOSonly';
l.update_rate = 10e-3;

c = l.get_channels;

assertEqual( size(c),[2,2] );

assertEqual( c(1,1).no_snap,201 ); % Interpolarion to default 10 m sampling time
assertEqual( c(1,2).no_snap,201 );
assertEqual( c(2,1).no_snap,201 ); 
assertEqual( c(2,2).no_snap,201 ); 

% Path loss after interpolation and merging of the freespace coefficients
d3d = sqrt(sum(( c(1,1).tx_position(:,ones(1,201)) - c(1,1).rx_position ).^2));
PL_free = 20*log10( d3d ) + 32.45 + 20*log10( l.simpar.center_frequency/1e9 );
assertTrue( all( abs( c(1,1).par.pg + PL_free ) < 1e-4 ) );

PG_coeff = 10*log10(abs(reshape(c(1,1).coeff,1,[])).^2);
assertTrue( all( abs( PG_coeff + PL_free ) < 1e-4 ) );

assertTrue( all( abs( reshape( c(1,1).delay , 1 , [] ) - d3d / l.simpar.speed_of_light ) < 1e-12 ) );


