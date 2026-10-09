function testLayout_PowerMap
%%
l = qd_layout.generate('regular',7,50);
l.no_tx = 2;
l.simpar.show_progress_bars = 0;

[map , posx , posy ] = l.power_map( 'BERLIN_UMa_LOS' , 'quick' , 50 );

assertEqual(  numel(map)  ,  2 );
assertEqual(  size(map{1})  ,  [10,10,1,3]  );

v = reshape( map{1} ,[],1);

assertTrue( isnumeric(v) && isreal(v) );
assertTrue( all(~isnan(v)) );
assertTrue( min(v)>0 );
assertTrue( max(v)<Inf );

p = qd_builder('BERLIN_UMa_LOS');
[map2 , posx , posy ] = power_map( l , p , 'sf' , 50 );

assertEqual(  numel(map2)  ,  2 );
assertEqual(  size(map2{1})  ,  [10,10,1,3]  );
v =  reshape( map2{1} ,[],1);
assertTrue( isnumeric(v) && isreal(v) );
assertTrue( all(~isnan(v)) );
assertTrue( min(v)>0 );
assertTrue( max(v)<Inf );

[map3 , posx , posy ] = power_map( l , {'BERLIN_UMa_LOS','Ul'} , 'phase' , 50 );

assertEqual(  numel(map3)  ,  2 );
assertEqual(  size(map3{1})  ,  [10,10,1,3]  );
v =  reshape( map3{1} ,[],1);
assertTrue( isnumeric(v) && ~isreal(v) );
assertTrue( all(~isnan(v)) );

assertTrue( all( reshape( abs( map{1} - abs(map3{1}).^2  ),1,[] ) < 1e-12 ));

p = qd_builder('BERLIN_UMa_LOS');
p(1,2) = qd_builder('BERLIN_UMa_LOS');

[map4 , posx , posy ] = power_map( l , p , 'detailed' , 50 );
assertEqual(  numel(map4)  ,  2 );
assertEqual(  size(map4{1})  ,  [10,10,1,3]  );
v =  reshape( map4{1} ,[],1);
assertTrue( isnumeric(v) && isreal(v) );
assertTrue( all(~isnan(v)) );
assertTrue( min(v)>0 );
assertTrue( max(v)<Inf );

