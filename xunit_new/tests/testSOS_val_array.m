function testSOS_val_array

set_rand_state( 1 );
ca = rand( 3,50 ) * 50;
cb = rand( 3,50 ) * 50;

% Generators with the same tables and different phases, mixed evaluation modes and distributions
x = qd_sos( 'Comb256', 'Normal', 10 );
x(1,2) = qd_sos( 'Comb256', 'Normal', 10 );
x(1,3) = qd_sos( 'Comb256', 'Uniform', 10 );
x(1,4) = qd_sos( 'Comb256', 'Normal', 10 );
init( x(1,4), 0 );                                  % Different maps
x(1,5) = qd_sos( 'Comb128', 'Normal', 10 );         % Other frequency table
x(1,6) = qd_sos( 'Comb256', 'Normal', 15 );         % Other decorrelation distance
x(1,7) = qd_sos( 'Comb256', 'Normal', 10 );
y = x(1,7);
y.dist_decorr = [10, 20];                            % Two decorrelation distances

for cc = 1 : 2
    if cc == 1
        s = val( x, ca );
    else
        s = val( x, ca, cb );
    end
    assertEqual( size(s), [7,50] );
    for n = 1 : 7
        if cc == 1
            r = val( x(1,n), ca );
        else
            r = val( x(1,n), ca, cb );
        end
        assertElementsAlmostEqual( s(n,:), r, 'absolute', 1e-12 );
    end
end

% Uniform values lie between 0 and 1
assertTrue( all( s(3,:) >= 0 & s(3,:) <= 1 ) );

% Same map: the value does not change when the devices swap their positions
s1 = val( x(1,1), ca, cb );
s2 = val( x(1,1), cb, ca );
assertElementsAlmostEqual( s1, s2, 'absolute', 1e-12 );

% Positions with fewer rows are zero-padded
s1 = val( x(1,1), ca(1:2,:) );
s2 = val( x(1,1), [ ca(1:2,:); zeros(1,50) ] );
assertElementsAlmostEqual( s1, s2, 'absolute', 1e-12 );

% Empty position list
assertEqual( size( val( x, zeros(3,0) ) ), [7,0] );

% Old preset names are mapped to the presets of Quadriga-Lib
y = qd_sos( 'Exp300' );
assertEqual( y.name, 'Exp256' );
assertEqual( y.no_coefficients, 256 );

% The phases follow the state of the MATLAB random number generator
set_rand_state( 2 );
a = qd_sos;
set_rand_state( 2 );
b = qd_sos;
assertEqual( a.sos_phase, b.sos_phase );

end
