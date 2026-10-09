function testArray_Calc_Gain

a = qd_arrayant('dipole');
[x,y] = a.calc_gain;
assertTrue( x - 10*log10(1.5) < 1e-3 )
assertTrue( y - 10*log10(1.5) < 1e-3 )
