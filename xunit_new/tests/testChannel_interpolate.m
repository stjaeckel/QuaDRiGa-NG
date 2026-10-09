function testChannel_interpolate
%%
dist  = 0:10;
distI = 0:0.1:10;

P = zeros(3,11);                           % Positions
P(1,:) = 10 + dist;
G = randn( 2,2,10 ) + 1j*randn( 2,2,10 );  % Coefficients
D = rand(10,1) * 5 * 1e-7;                 % Delays
S = 1 + dist*0.01;                         % Scale from 1 to 1.1

D = D*S;
G = G(:,:,:,ones(1,11));
for n = 1:11
   G(:,:,:,n) = G(:,:,:,n) * S(n);
end

c = qd_channel( G,D );
c.rx_position = P;

for n = 1:4
    switch n
        case 1
            ci = c.interpolate( distI );
        case 2
            ci = c.interpolate( reshape( [1;1;1;1] * distI,2,2,[] ) );
        case 3
            ci = c.interpolate( distI ,'cubic' );
        case 4
            ci = c.interpolate( reshape( [1;1;1;1] * distI,2,2,[] ) ,'cubic' );
    end
    assertEqual( size( ci.rx_position ) , [3 101] )
    assertEqual( size( ci.coeff ) , [2,2,10,101] )
    assertEqual( size( ci.delay ) , [10,101] )
    assertTrue( all( abs( ci.rx_position(1,:) - 10 - distI ) < 1e-5 ) )
    assertTrue( all( abs( reshape( ci.coeff( :,:,:,1:10:end ) , [] , 1 ) - G(:) ) < 1e-6 ) )
    assertTrue( all( abs( reshape( ci.delay( :,1:10:end ) , [] , 1 ) - D(:) ) < 1e-6 ) )
end

c.individual_delays = true;


for n = 1:4
    switch n
        case 1
            ci = c.interpolate( distI );
        case 2
            ci = c.interpolate( reshape( [1;1;1;1] * distI,2,2,[] ) );
        case 3
            ci = c.interpolate( distI ,'cubic' );
        case 4
            ci = c.interpolate( reshape( [1;1;1;1] * distI,2,2,[] ) ,'cubic' );
    end
    assertEqual( size( ci.rx_position ) , [3 101] )
    assertEqual( size( ci.coeff ) , [2,2,10,101] )
    assertEqual( size( ci.delay ) , [2,2,10,101] )
    assertTrue( all( abs( ci.rx_position(1,:) - 10 - distI ) < 1e-5 ) )
    assertTrue( all( abs( reshape( ci.coeff( :,:,:,1:10:end ) , [] , 1 ) - G(:) ) < 1e-6 ) )
    assertTrue( all( abs( reshape( ci.delay( 1,1,:,1:10:end ) , [] , 1 ) - D(:) ) < 1e-6 ) )
end

