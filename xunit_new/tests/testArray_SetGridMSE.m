function testArray_SetGridMSE
%%

[ theta, phi ] = qf.pack_sphere( 16 );        

N = numel( theta );
a = qd_arrayant('custom',20,20,0.05);               % Main beam opening and front-back ratio
a.set_grid( (-180:10:180)*pi/180, (-90:10:90)*pi/180  );

assertEqual( a.no_az, 37 );
assertEqual( a.no_el, 19 );

a.element_position(1) = 0.2;                     % Distance from phase-center
a.copy_element(1,2:N+1);
for n = 1:N                                         % Create sub-elements
    a.rotate_pattern( theta(n)*180/pi,'y',n,1);
    a.rotate_pattern( phi(n)*180/pi,'z',n,1);
end
a.center_frequency = 299792458/0.125;
a.combine_pattern;
P = sum( abs(a.Fa(:,:,1:N)).^2,3 );
a.Fa(:,:,1:N) = a.Fa(:,:,1:N) ./ sqrt(P(:,:,ones(1,N)));
a.Fb(:,:,N+1)=1;
a.Fa(:,:,N+1)=0;

a(1,2) = a(1,1);

mse = set_grid( a, (-180:15:180)*pi/180, (-90:15:90)*pi/180  );

assertEqual( a.no_az, 25 );
assertEqual( a.no_el, 13 );

assertTrue( qf.eqo(a(1,1),a(1,2)));
assertEqual( size(mse), [1,2] )
assertEqual( mse(1), mse(2) )

