function testChannel_fr

G = randn( 2,2,10,3 ) + 1j*randn( 2,2,10,3 );
D = rand(10,3) * 5 * 1e-7;
c = qd_channel(G, D);

H = c.fr( 20e6, 20 );
assertEqual( size(H), [2,2,20,3] )
H = c.fr( 20e6, 20, 2 );
assertEqual( size(H), [2,2,20] )
H = c.fr( 20e6, [0,0.5,1], 1 );
assertEqual( size(H), [2,2,3] )

c.individual_delays = 1;
H = c.fr( 20e6, [0,0.5,1], 3 );
assertEqual( size(H), [2,2,3] )
