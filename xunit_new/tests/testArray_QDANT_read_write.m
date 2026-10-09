function testArray_QDANT_read_write
%%

a = qd_arrayant('ula2');
a.center_frequency = 1.9e9;
a.Fa = randn( size( a.Fa ) ) + 1j*randn( size( a.Fa ) );
a.Fb = randn( size( a.Fa ) ) + 1j*randn( size( a.Fa ) );
a.coupling = randn( size( a.coupling ) ) + 1j*randn( size( a.coupling ) );

a(1,2) = qd_arrayant('xpol');
a(1,2) = qd_arrayant('ula2');
a(1,2).Fa = randn( size( a(1,2).Fa ) ) + 1j*randn( size( a(1,2).Fa ) );

xml_write(a,'tst.qdant')

b = qd_arrayant.xml_read('tst.qdant');

assertEqual( size(b), size(a) );

for n = 1:2
    assertEqual( a(1,n).name, b(1,n).name );
    assertTrue( abs( a(1,n).center_frequency - b(1,n).center_frequency ) < 1e-5);
    assertTrue( all( abs( a(1,n).Fa(:)-b(1,n).Fa(:) ) < 1e-2 ) )
    assertTrue( all( abs( a(1,n).Fb(:)-b(1,n).Fb(:) ) < 1e-2 ) )
    assertTrue( all( abs( a(1,n).coupling(:)-b(1,n).coupling(:) ) < 1e-2 ) )
end

delete('tst.qdant')
