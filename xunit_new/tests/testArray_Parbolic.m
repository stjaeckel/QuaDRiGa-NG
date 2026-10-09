function testArray_Parbolic

% Test if the correct gain is obtained from the paatterns
r = [ 0.5  10   0.5 10 ];
f = [   1   1   100 40 ];
for n = 1:numel(r)
    a = qd_arrayant.generate('parabolic', r(n) , f(n)*1e9 , -50 , 1 );
    
    % Ideal gain formula
    Gt = 10*log10(( 2 * pi * r(n)/( qd_simulation_parameters.speed_of_light/a.center_frequency )).^2);
    
    % Gain calculated from the pattern
    Gp = a.calc_gain;
    
    % Compare
    assertTrue( abs( Gt - Gp ) < 0.1 );
    
    % Test for circular polarization
    if n == 1 
        a = qd_arrayant.generate('parabolic', r(n) , f(n)*1e9 , -50 , 3 );
        Gp = a.calc_gain;
        assertTrue( abs( Gt - Gp ) < 0.1 );
    end
end

% Test pattern rotations
a = qd_arrayant.generate('parabolic', 0.5 , 2e9 , -50 , 3 );

% Test if there is a partial grid
d = diff(a.azimuth_grid);
assertTrue( all( d-d(1) < 1e-12 ) );
assertTrue( sum(d) < pi );

% Test if there is a partial grid
d = diff(a.elevation_grid);
assertTrue( all( d-d(1) < 1e-12 ) );
assertTrue( sum(d) < pi );

% Store data
Fa = a.Fa;
Fb = a.Fb;
az = a.azimuth_grid;
el = a.elevation_grid;
G = a.calc_gain;

a.rotate_pattern(-80,'y');
assertTrue( abs( a.calc_gain - G ) < 0.1 );

% Test if there is a full azimuth grid
d = diff(a.azimuth_grid);
assertTrue( all( d-d(1) < 1e-12 ) );
assertTrue( sum(d) - 2*pi < 1e-12 );

a.rotate_pattern(90,'z');
assertTrue( abs( a.calc_gain - G ) < 0.1 );

a.rotate_pattern(-80,'x');
assertTrue( abs( a.calc_gain - G ) < 0.1 );

a.rotate_pattern(270,'z');
assertTrue( abs( a.calc_gain - G ) < 0.1 );

% Test if there is a partial grid
d = diff(a.azimuth_grid);
assertTrue( all( d-d(1) < 1e-12 ) );
assertTrue( sum(d) < pi );

% Test if there is a partial grid
d = diff(a.elevation_grid);
assertTrue( all( d-d(1) < 1e-12 ) );
assertTrue( sum(d) < pi );

% Restore original grid
a.set_grid( az, el );

% Compare patterns
assertTrue( all( abs( a.Fa(:) - Fa(:) ) < 0.6 ) );
assertTrue( all( abs( a.Fb(:) - Fb(:) ) < 0.6 ) );

