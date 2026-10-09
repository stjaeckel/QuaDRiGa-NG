function testArray_Rotate_Combine

% Rotate - Combine
a = qd_arrayant('patch',90,90,0.01);
a.copy_element(1,2);
a.rotate_pattern(90,'x',2);
a.rotate_pattern(-90,'y')
a.coupling = [1;1j];
a.combine_pattern(2.6e9);

% Combine - Rotate
b = qd_arrayant('patch',90,90,0.01);
b.copy_element(1,2);
b.rotate_pattern(90,'x',2);
b.coupling = [1;1j];
b.combine_pattern(2.6e9)
b.rotate_pattern(-90,'y')

errA = a.Fa(2:end-1,2:end-1) - b.Fa(2:end-1,2:end-1);
errB = a.Fb(2:end-1,2:end-1) - b.Fb(2:end-1,2:end-1);
errB(90,180) = 0;  % Pole

assertTrue( max(abs(errA(:))) < 1e-3 )
assertTrue( max(abs(errB(:))) < 1e-3 )
