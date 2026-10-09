function testArray_Subarray

a = qd_arrayant( 'ula4' );
b = a.sub_array( [2,4] );

assertEqual( b.no_elements , 2 );
assertEqual( b.element_position , a.element_position(:,[2,4]) );
assertEqual( b.coupling , eye(2) );