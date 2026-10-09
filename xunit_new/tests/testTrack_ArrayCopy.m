function testTrack_ArrayCopy
%%
a = qd_track;
a.name = 'Bla1';
a.no_snapshots = 3;

a(1,1,2) = qd_track('circular');
a(1,1,2).name = 'Bla2';

a(1,1,3) = a(1,1,1);

tst = true(1,1,3);
tst(1,1,2) = false;

assertEqual( qf.eqo( a(1,1,1), a )  , tst  );

b = copy(a);

assertEqual( numel(b)  , 3  );
assertTrue( ~isequal( b(1),b(2) ) )
assertTrue( isequal( b(1),b(3) ) )

assertEqual( qf.eqo( b(1,1,1), b )  , tst  );
assertEqual( qf.eqo( b(1,1,1), a )  , false(1,1,3)  );

names = {'name','initial_position','no_snapshots','positions',...
  'movement_profile','no_segments','segment_index','scenario','closed'};

for n=numel(a)
    for m = 1:numel(names)
        assertEqual(  a(n).( names{m} ) , b(n).( names{m} ) );
    end
end

