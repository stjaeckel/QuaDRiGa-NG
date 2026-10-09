function testArraySetFunctions
% Thest the get and set-Interface

a = qd_arrayant('omni');
a.no_elements = 10;
assertEqual(  size(a.element_position) , [3 10] );
assertEqual(  size(a.Fa)  ,  [181 361 10] );
assertEqual(  size(a.Fb)  ,  [181 361 10] );
assertEqual(  size(a.coupling)  ,  [10 10] );

a.no_elements = 2;
assertEqual(  size(a.element_position)  ,  [3 2] );
assertEqual(  size(a.Fa)  ,  [181 361 2] );
assertEqual(  size(a.Fb)  ,  [181 361 2] );
assertEqual(  size(a.coupling)  ,  [2 2] );

a.set_grid([1 0 1],0,0)
assertEqual(  size(a.Fa)  ,  [1 3 2] );
assertEqual(  size(a.Fb)  ,  [1 3 2] );
assertEqual(  a.no_az  , 3   );
assertEqual(  a.no_el  , 1  );

a.set_grid((-180:180)*pi/180,(-90:90)*pi/180,0)
assertEqual(  size(a.element_position)  ,  [3 2] );
assertEqual(  size(a.Fa)  ,  [181 361 2] );
assertEqual(  size(a.Fb)  ,  [181 361 2] );
assertEqual(  size(a.coupling)  ,  [2 2] );

