function testMesh_read_obj

% Write a small OBJ file: two walls, each made of two triangles
fn_obj = 'test_mesh.obj';
fn_mtl = 'test_mesh.mtl';

fid = fopen( fn_mtl,'w' );
fprintf( fid,'newmtl itu_concrete\nKd 0.5 0.6 0.7\nnewmtl itu_glass\nKd 1 0 0\n' );
fclose( fid );

fid = fopen( fn_obj,'w' );
fprintf( fid,'mtllib %s\n',fn_mtl );
fprintf( fid,'o wall1\nv 5 -5 -5\nv 5 5 -5\nv 5 5 5\nv 5 -5 5\nusemtl itu_concrete\nf 1 2 3\nf 1 3 4\n' );
fprintf( fid,'o wall2\nv 8 -5 -5\nv 8 5 -5\nv 8 5 5\nv 8 -5 5\nusemtl itu_glass\nf 5 6 7\nf 5 7 8\n' );
fclose( fid );

m = qd_mesh;
m.read_obj( fn_obj );

delete( fn_obj );
delete( fn_mtl );

% The first object and the first material are the defaults of the empty mesh (without faces)
assertEqual( m.no_face, 4 );
assertEqual( m.no_vert, 8 );
assertEqual( m.no_obj, 3 );
assertEqual( m.no_mtl, 3 );
assertEqual( m.obj_name(2:3), {'wall1','wall2'} );
assertEqual( m.mtl_name(2:3), {'itu_concrete','itu_glass'} );
assertEqual( double( m.obj_index ), [2,2,3,3] );
assertEqual( double( m.mtl_index ), [2,2,3,3] );

% Diffuse color from the MTL file
assertElementsAlmostEqual( m.mtl_color(:,2:3), [0.5,1 ; 0.6,0 ; 0.7,0], 'absolute', 1e-6 );

% Electric properties from the ITU material table (Rec. ITU-R P.2040)
assertElementsAlmostEqual( m.mtl_prop(:,2), [5.24; 0; 0.0462; 0.7822; 0], 'absolute', 1e-6 );
assertElementsAlmostEqual( m.mtl_prop(:,3), [6.31; 0; 0.0036; 1.3394; 0], 'absolute', 1e-6 );

% Ray 1 passes both walls, ray 2 misses them, ray 3 ends between the walls
orig = [0;0;0];
dest = [ 10,10,6 ; 1,20,1 ; 0.5,0.1,0.5 ];
[ islos, no_trans, fbs, sbs, iFBS, iSBS ] = m.intersect_mesh( orig, dest );

assertEqual( islos(:)', [false,true,false] );
assertEqual( double( no_trans(:)' ), [2,0,1] );
assertElementsAlmostEqual( fbs(1,:), [5,0.5,0.25], 'absolute', 1e-5 );
assertElementsAlmostEqual( sbs(1,:), [8,0.8,0.4], 'absolute', 1e-5 );
assertElementsAlmostEqual( fbs(3,:), [5,5/6,5/12], 'absolute', 1e-5 );
assertTrue( all( double( m.obj_index( iFBS([1,3]) ) ) == 2 ) );
assertEqual( double( m.obj_index( iSBS(1) ) ), 3 );
assertEqual( double( iFBS(2) ), 0 );

end
