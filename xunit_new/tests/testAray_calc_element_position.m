function testAray_calc_element_position
% Generate omni
%%

a = qd_arrayant('dipole');
a.center_frequency = 3.7e9;
a.element_position = [ 0.7 ; -0.3 ; 0.18 ];

pos = a.element_position;
Fa = a.Fa;

a.combine_pattern;

Fa2 = a.Fa;

elp = a.calc_element_position(0);

assertTrue( all( abs( elp-pos ) < 1e-3 ))
assertTrue( all(abs(a.Fa(:)-Fa(:)) < 0.1) )

a.combine_pattern;

assertTrue( all(abs(a.Fa(:)-Fa2(:)) < 1e-12) )


