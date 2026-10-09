function testSOS_copy
%%

x = qd_sos;
x(2) = qd_sos;
x(3) = x(1);

assertTrue( ~isequal( x(1),x(2) ) )
assertTrue( isequal( x(1),x(3) ) )

y = copy(x);

assertTrue( all( size(y) == [1,3] ) )
assertTrue( ~isequal( y(1),y(2) ) )
assertTrue( isequal( y(1),y(3) ) )
