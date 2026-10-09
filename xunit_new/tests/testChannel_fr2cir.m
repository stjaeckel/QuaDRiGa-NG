function testChannel_fr2cir

a = qd_arrayant;
a.copy_element(1,2);
a.element_position(1,:) = [-0.5,0.5];

l = qd_layout;
l.simpar.use_3GPP_baseline = 1;
l.simpar.samples_per_meter = 1;
l.simpar.use_absolute_delays = 1;
l.simpar.show_progress_bars = 0;

l.no_rx = 1;
l.randomize_rx_positions( 10, 1.5, 1.5, 1.005);
l.track(1,1).set_scenario( 'WINNER_UMa_C2_LOS' )
l.rx_array = a;
l.rx_position(1) = l.rx_position(1) + 100;

l.no_tx = 1;
l.tx_position(3,:) = 25;
l.tx_array = a;

b = l.init_builder;

b.sf = 1;
b.xpr = Inf;
b.scenpar.PerClusterAS_A = 0;
b.scenpar.PerClusterES_A = 0;
b.scenpar.PerClusterAS_D = 0;
b.scenpar.PerClusterES_D = 0;

b.gen_parameters;

c = b.get_channels;

b.taus = (0:7)*1e-7;
b.pow  = [ 0.79 ones(1,7)*0.03 ];
b.pin(:) = 0;

c = b.get_channels;

BW = 100e6;
N  = 200;

Y = c.fr( BW, N );
d = qd_channel.fr2cir( Y, BW, 8, [],[],[],[],[], 0 );

if 0
    s = 2;
    
    delay_axis = 1e9*(1:200)/BW;
    dst = sqrt(sum(( b.rx_track.positions_abs - b.tx_position ).^2));
    delays_bld = 1e9* ( dst(s)/b.simpar.speed_of_light + b.taus );
    
    P = squeeze(mean(mean(abs(d.coeff(:,:,:,s)).^2,1),2));
    D = d.delay(:,s)*1e9;
    
    plot( delay_axis, 10*log10( d.par.PDPo(:,s) ) )
    hold on
    plot( delay_axis, 10*log10( d.par.PDPr(:,s) ),'r' )
    plot( delay_axis, 10*log10( d.par.PDPn(:,s) ),'--k' )
    plot( delays_bld, 10*log10(b.gain) , '+m','Markersize',8 )
    plot( D , 10*log10(P) , 'ob','Markersize',8 )
    hold off
    
    grid on
    xlabel('Delay (ns)')
    ylabel('Gain (dB)')
    title(s)
end

assertTrue( all( 10*log10(d.par.snr_est) > 30 ) );

% Delay resolution < 2 ns
[dd,ii] = sort( d.delay,1 );

assertTrue( all( c.delay(:) - dd(:) < 2e-9 ) );

% Coefficient resolution
dc = cat(4, d.coeff(:,:,ii(:,1),1), d.coeff(:,:,ii(:,2),2));

SE = abs( c.coeff(:) - dc(:) ).^2 ./ abs(c.coeff(:)).^2;
assertTrue( all( 10*log10(SE) < -12 ) );

