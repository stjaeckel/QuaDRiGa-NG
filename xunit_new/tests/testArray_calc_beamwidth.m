function testArray_calc_beamwidth

% Custom antenna with 20 deg azimuth and 30 deg elevation FWHM
a = qd_arrayant('custom',20,30,0);
[ bw_az, bw_el, az, el ] = a.calc_beamwidth;
assertEqual( size(bw_az), [1,1] );
assertTrue( abs( bw_az - 20 ) < 0.1 );
assertTrue( abs( bw_el - 30 ) < 0.1 );
assertTrue( abs( az ) < 0.01 );
assertTrue( abs( el ) < 0.01 );

% A larger threshold must give a wider beam
[ bw_az10, bw_el10 ] = a.calc_beamwidth( [], 10 );
assertTrue( bw_az10 > bw_az + 5 );
assertTrue( bw_el10 > bw_el + 5 );

% Pointing angle must follow the rotation of the pattern
a.rotate_pattern( 30, 'z' );
[ bw_az, bw_el, az, el ] = a.calc_beamwidth;
assertTrue( abs( bw_az - 20 ) < 0.1 );
assertTrue( abs( bw_el - 30 ) < 0.1 );
assertTrue( abs( az - 30 ) < 0.1 );
assertTrue( abs( el ) < 0.1 );

% Element selection
a = qd_arrayant('custom',20,30,0);
b = qd_arrayant('custom',40,60,0);
a.append_array( b );
[ bw_az, bw_el ] = a.calc_beamwidth;
assertEqual( size(bw_az), [2,1] );
assertEqual( size(bw_el), [2,1] );
assertTrue( all( abs( bw_az - [20;40] ) < 0.2 ) );
assertTrue( all( abs( bw_el - [30;60] ) < 0.2 ) );

[ bw_az, bw_el ] = a.calc_beamwidth( [2,1,2] );
assertTrue( all( abs( bw_az - [40;20;40] ) < 0.2 ) );
assertTrue( all( abs( bw_el - [60;30;60] ) < 0.2 ) );

% Invalid element index
try
    a.calc_beamwidth( 3 );
    error('moxunit:exceptionNotRaised', 'Expected an error!');
catch ME
    if strcmp( ME.identifier, 'moxunit:exceptionNotRaised' )
        error('moxunit:exceptionNotRaised', 'Expected an error!');
    end
end

end
