function testBuilder_gen_cdl_model

% Generate
b = qd_builder.gen_cdl_model( 'nr-cdl-a', 2.6e9, 1, 1, 30 );
c = b.get_channels;

pow = reshape( 10*log10(mean(abs(c.coeff).^2,4)), 1,[] );
assertTrue( isinf( pow(1)))

% Test power values
tmp = pow(2:end) - [ -13.4 , 0, -2.2, -4, -6, -8.2, -9.9, -10.5, -7.5, -15.9, -6.6, -16.7, -12.4,...
    -15.2, -10.8, -11.3, -12.7, -16.2, -18.3, -18.9, -16.6, -19.9, -29.7 ];
assertTrue( all( abs(tmp)<1e-12 ));

% Test if DS is 30 ns
assertTrue( abs( qf.calc_delay_spread( c.delay(:,1)', 10.^(0.1*pow) ) - 30e-9 ) < 1e-20 )

% Test departure angles
a = qd_arrayant('testarray',5);
a = a.sub_array(1:29);  % CDL models do not support circular polarization

b = qd_builder.gen_cdl_model( 'nr-cdl-a', 2.6e9, 0, 1, 30, [],[],[],[],[],0 );
b.tx_array = a;
c = b.get_channels;

[ az, el ] = qf.calc_angles( permute(c.coeff,[2,1,3]), a ,1,[],[],1,0 );      % Calculate angles
az = az * 180/pi;
zd = 90 - el * 180/pi;

tmp = az(2:end) - [ -178.1, -4.2, -4.2, -4.2, 90.2, 90.2, 90.2, 121.5, -81.7, 158.4, -83, 134.8, -153,...
    -172, -129.9, -136, 165.4, 148.4, 132.7, -118.6, -154.1, 126.5, -56.2 ];
assertTrue( all( abs(tmp) < 1 ));

tmp = zd(2:end) - [ 50.2, 93.2, 93.2, 93.2, 122, 122, 122, 150.2, 55.2, 26.4, 126.4, 171.6, 151.4, 157.2,...
    47.2, 40.4, 43.3, 161.8, 10.8, 16.7, 171.7, 22.7, 144.9 ];
assertTrue( all( abs(tmp) < 1 ));

% Test the polarizazion
co = sum(abs(c.coeff(1,1:28,:)).^2,2);
cross = abs(c.coeff(1,29,:)).^2;
tmp = 10*log10(co(:)./cross(:))-10;
assertTrue( all( abs(tmp(2:end)) < .05 ));

% Test arrival angles
b = qd_builder.gen_cdl_model( 'nr-cdl-a', 2.6e9, 0, 1, 30, [],[],[],[],[],0 );
b.rx_array = a;
c = b.get_channels;

[ az, el ] = qf.calc_angles( c.coeff, a ,1,[],[],1,0 );       % Calculate angles
az = az * 180/pi;
zd = 90 - el * 180/pi;

tmp = az(2:end) - [ 51.3, -152.7, -152.7, -152.7, 76.6, 76.6, 76.6, -1.8, -41.9, 94.2, 51.9, -115.9, ...
    26.6, 76.6, -7, -23, -47.2, 110.4, 144.5, 155.3, 102, -151.8, 55.2 ];
assertTrue( all( abs(tmp) < 1 ));

tmp = zd(2:end) - [ 125.4, 91.3, 91.3, 91.3, 94, 94, 94, 47.1, 56, 30.1, 58.8, 26, 49.2, 143.1, 117.4,...
    122.7, 123.2, 32.6, 27.2, 15.2, 146, 150.7, 156.1 ];
assertTrue( all( abs(tmp) < 1 ));

% Test the polarizazion
co = sum(abs(c.coeff(1:28,1,:)).^2,1);
cross = abs(c.coeff(29,1,:)).^2;
tmp = 10*log10(co(:)./cross(:))-10;
assertTrue( all( abs(tmp(2:end)) < .05 ));

% Set departure angle spreads
b = qd_builder.gen_cdl_model( 'nr-cdl-a', 2.6e9, 0, 1, 30,[], 10,[],5,[],0 );
assertTrue( abs(b.asD - 10)<1e-6 );
assertTrue( abs(b.esD - 5)<1e-6 );
b.tx_array = a;
c = b.get_channels;

[ az, el ] = qf.calc_angles( permute(c.coeff,[2,1,3]), a ,1,[],[],1,0 );  
pow = sum(abs(c.coeff(:,1:28,:)).^2,2);
assertTrue(  abs(qf.calc_angular_spreads( az(:)', pow(:)' )*180/pi - 10) < 0.1 );
assertTrue(  abs(qf.calc_angular_spreads( el(:)', pow(:)' )*180/pi - 5) < 0.1 );

% Set arrival angle spreads
b = qd_builder.gen_cdl_model( 'nr-cdl-a', 2.6e9, 0, 1, 30,[],[],9,[],7,0 );
assertTrue( abs(b.asA - 9)<1e-6 );
assertTrue( abs(b.esA - 7)<1e-6 );
b.rx_array = a;
c = b.get_channels;

[ az, el ] = qf.calc_angles( c.coeff, a, 1,[],[],1,0 ); 
pow = sum(abs(c.coeff(1:28,:,:)).^2,1);
assertTrue(  abs(qf.calc_angular_spreads( az(:)', pow(:)' )*180/pi - 9) < 0.1 );
assertTrue(  abs(qf.calc_angular_spreads( el(:)', pow(:)' )*180/pi - 7) < 0.1 );

% Test LOS models
b = qd_builder.gen_cdl_model( 'nr-cdl-e', 2.6e9, 0, 1, 45, -3, [],[],[],[], 0.1 );
b.tx_array = a;
c = b.get_channels;
kf = sum(abs(c.coeff(:,1:28,1)).^2,2) ./ sum(sum(abs(c.coeff(:,1:28,2:end)).^2,2),3);
assertTrue( abs(kf - 0.5)<0.01 );

% Test TDL NLOS model
b = qd_builder.gen_cdl_model( 'nr-tdl-e', 2.6e9, 2, 5, 45 );
c = b.get_channels;

w = c.no_snap;                  % Doppler analysis windows size (100 ms)
update_rate = 5/c.no_snap;     % Channel update rate
Doppler_axis = -( (0:w-1)/(w-1)-0.5)/update_rate;     % The Doppler axis in Hz

H = permute( c.coeff, [3,4,1,2] );
G = ifft(H,[],2);                                   % 2D IFFT
G = fftshift( G,2);
Doppler = 10*log10(sum(abs(G(:,:)).^2,1));
[~,ii] = max(Doppler);

fD_max = 2/b.simpar(1,1).wavelength;

assertTrue( abs(Doppler_axis(ii)/fD_max + 0.7) < 0.01 );
%plot( Doppler_axis, 10*log10(sum(abs(G(:,:)).^2,1)),'k' )

pow = reshape( 10*log10(mean(abs(c.coeff).^2,4)), 1,[] );
assertTrue( abs( qf.calc_delay_spread( c.delay(:,1)', 10.^(0.1*pow) ) - 45e-9 ) < 1e-20 )



