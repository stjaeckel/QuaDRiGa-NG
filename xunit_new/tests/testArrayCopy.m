function testArrayCopy
% Test the copy function

a = qd_arrayant;
a.name = 'Bla1';
a.no_elements = 2;
a.Fb = rand( size(a.Fb));
a.Fa = rand( size(a.Fa));
a.coupling = rand( size(a.coupling));

a(1,2) = qd_arrayant;
a(1,2).name = 'Bla2';
a(1,2).Fb = rand( size(a(2).Fb));
a(1,2).Fa = rand( size(a(2).Fa));
a(1,2).coupling = rand( size(a(2).coupling));

a(1,3) = a(1,1);

set_grid( a, (-180:10:180)*pi/180 , (-90:10:90)*pi/180 );

b = copy( a );

names = {'name','no_elements','elevation_grid','azimuth_grid',...
  'element_position','Fa','Fb','coupling'};

assertEqual( numel(b)  , 3  );
assertTrue( ~isequal( b(1),b(2) ) )
assertTrue( isequal( b(1),b(3) ) )


for n=numel(a)
    for m = 1:numel(names)
        assertEqual(  a(n).( names{m} ) , b(n).( names{m} ) );
    end
end
