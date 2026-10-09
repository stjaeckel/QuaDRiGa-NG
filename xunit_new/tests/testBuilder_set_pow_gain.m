function testBuilder_set_pow_gain

b = qd_builder('3GPP_38.901_UMi_LOS');

assertTrue( isempty( b.pow ) );

b.simpar.center_frequency(2) = 22e9;

assertTrue( isempty( b.pow ) );

b.tx_position = [0;0;1.5]*[1 1 1];
b.rx_positions = [100 200 300;0 0 0 ;1.5 1.5 1.5];

assertTrue( isempty( b.pow ) );

T = [0.5,0.4,0.6,0.3];

b.pow = T;

% Size of gain and pow must match
assertEqual( size( b.gain ), [3,4,2] );
assertEqual( size( b.pow ), [3,4,2] );

% Pow must be normalized to 1
X = sum(b.pow,2);   
assertTrue( all(abs(X(:)-1) < 1e-13) );

% Relative differences must match
X = b.pow - repmat(T/sum(T),[b.no_rx_positions,1,b.no_freq]);
assertTrue( all(abs(X(:)) < 1e-13) );

% Check if gain is normalized to PG
PL = b.get_pl;
X = 10*log10(sum(b.gain,2)) + reshape(PL',3,1,2);
assertTrue( all(abs(X(:)) < 1e-13) );

% Set a shadow fading value
SF = 10*rand( 2,3 );
b.sf = SF;

% Check if gain is normalized to PG --> SF did not change the gains
PL = b.get_pl;
X = 10*log10(sum(b.gain,2)) + reshape(PL',3,1,2);
assertTrue( all(abs(X(:)) < 1e-13) );

% Powers should be scaled down by SF
X = sum(b.pow,2) .* reshape(SF',3,1,2);
assertTrue( all(abs(X(:)-1) < 1e-13) );

% Set power again --> should include SF in the gain
b.pow = T;
X = sum(b.pow,2);   
assertTrue( all(abs(X(:)-1) < 1e-13) );

% Check if gain contains SF and PL
X = 10*log10(sum(b.gain,2)) + reshape(PL',3,1,2) - 10*log10(reshape(SF',3,1,2));
assertTrue( all(abs(X(:)) < 1e-13) );

