function testArrayAppend

a = qd_arrayant('ula2');
b = qd_arrayant('ula4');
a.append_array( b );
assertTrue( a.no_elements == 6 );