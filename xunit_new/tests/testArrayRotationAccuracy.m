function testArrayRotationAccuracy

a = qd_arrayant('custom',10,10,0.1);
a.rotate_pattern(-33.5,'y');
a.rotate_pattern(90,'z');

b = qd_arrayant('custom',10,10,0.1);
b.rotate_pattern(90,'z');
b.rotate_pattern(33.5,'x');

x = a.Fa(2:end-1,:) - b.Fa(2:end-1,:);  
assertTrue( all( abs(x(:)) < 1e-12 ) );

x = a.Fb(2:end-1,:) - b.Fb(2:end-1,:);  
assertTrue( all( abs(x(:)) < 1e-12 ) );
