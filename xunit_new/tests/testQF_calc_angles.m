function testQF_calc_angles
%%

% Build an antenna with good angular resolution
[ theta, phi ] = qf.pack_sphere( 30 );
N = numel( theta );
a = qd_arrayant('custom',20,20,0.05);               % Main beam opening and front-back ratio
a.set_grid( (-180:10:180)*pi/180, (-90:10:90)*pi/180  );
a.element_position(1) = 0.2;                     % Distance from phase-center
a.copy_element(1,2:N+1);
for n = 1:N                                         % Create sub-elements
    a.rotate_pattern( theta(n)*180/pi,'y',n,1);
    a.rotate_pattern( phi(n)*180/pi,'z',n,1);
end
a.center_frequency = 299792458/0.125;

P = sum( abs(a.Fa(:,:,1:N)).^2,3 );
a.Fa(:,:,1:N) = a.Fa(:,:,1:N) ./ sqrt(P(:,:,ones(1,N)));
a.Fb(:,:,N+1)=1;
a.Fa(:,:,N+1)=0;

a1 = a.copy;

a.combine_pattern;

%% Create channel coefficients
b = qd_builder('3GPP_38.901_UMa_NLOS');
b.scenpar.SC_lambda = 0;
b.scenpar.XPR_mu  = 6;
b.scenpar.XPR_sigma  = 0;
b.scenpar.NumSubPaths = 1;
b.simpar.center_frequency = a.center_frequency;
b.simpar.use_3GPP_baseline = 1;
b.simpar.show_progress_bars = 0;
b.tx_array = a;
b.rx_array = qd_arrayant('omni');
b.tx_position = [0;0;25];
b.rx_positions = rand(3,10)*100;
b.rx_positions(1,:) = b.rx_positions(1,:) +20;
b.rx_positions(3,:) = b.rx_positions(3,:)*0.1;
b.gen_parameters;


b.EoD( b.EoD > 1 ) = 1;
b.EoD( b.EoD < -1 ) = -1;

c = b.get_channels;

for n = 1 : b.no_rx_positions
    b.tx_array(1,n) = a1;
end

d = b.get_channels;

%%
element_position = a1.element_position;
for n = 1 : b.no_rx_positions
    
    % Test if the estimated angles are correct
    [az,el,J] = qf.calc_angles( permute( c(1,n).coeff, [2,1,3] ), a,1,[],[],1,0 );
    
    assertTrue( all( abs( exp(1j*b.AoD(n,:)) - exp(1j*az) ) < 0.1 ) )
    assertTrue( all( abs( b.EoD(n,:) - el )*180/pi < 5 ) )
    
    % Test if the same results can be achieved with the antenna that includes the element positions
    if n == 1
        [az1,el1] = qf.calc_angles( permute( c(1,n).coeff, [2,1,3] ), a1,1,[],[],1,0 );
        assertTrue( all( abs( az-az1  ) < 0.002 ) );
        assertTrue( all( abs( el-el1  ) < 0.002 ) );
        assertEqual( element_position, a1.element_position );
        
        [az2,el2] = qf.calc_angles( permute( d(1,n).coeff, [2,1,3] ), a,1,[],[],1,0 );
        assertTrue( all( abs( exp(1j*b.AoD(n,:)) - exp(1j*az2) )*180/pi < 5 ) )
        assertTrue( all( abs( b.EoD(n,:) - el2 )*180/pi < 5 ) )
    end
    
    J = permute( J , [3,1,2] );
    J = mean( abs(J(2:end,1).^2 ) ./ abs(J(2:end,2).^2 ) );
    
    assertTrue( J > 3.7 & J < 4.3);
end

