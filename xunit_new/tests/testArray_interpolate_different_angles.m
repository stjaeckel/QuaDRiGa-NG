function testArray_interpolate_different_angles
%%

a = qd_arrayant('patch'); 

a.Fa = a.Fa ./ max(a.Fa(:));

a.copy_element(1,2:3);
a.Fa(:,:,2) = 2*a.Fa(:,:,2);
a.Fa(:,:,3) = 3*a.Fa(:,:,3);

a.Fb = -1j*a.Fa;

a.element_position(1,:) = [ 1 , 2 , 3 ];

% Interpolater broadside
[V,H,D] = a.interpolate( 0,0 );
assertTrue( all( abs( V - [1;2;3] ) < 1e-14 ))
assertTrue( all( abs( H + 1j*[1;2;3] ) < 1e-14 ))
assertTrue( all( abs( D + [1;2;3] ) < 1e-14 ))

% Select subset
[V,H,D] = a.interpolate( 0,0,2:3 );
assertTrue( all( abs( V - [2;3] ) < 1e-14 ))
assertTrue( all( abs( H + 1j*[2;3] ) < 1e-14 ))
assertTrue( all( abs( D + [2;3] ) < 1e-14 ))

% Change orientation z
[V,H,D] = a.interpolate( pi/2,0,[],[0;0;pi/2] );
assertTrue( all( abs( V - [1;2;3] ) < 1e-14 ))
assertTrue( all( abs( H + 1j*[1;2;3] ) < 1e-14 ))
assertTrue( all( abs( D + [1;2;3] ) < 1e-14 ))

% Change orientation y
[V,H,D] = a.interpolate( 0,pi/2,[],[0;pi/2;0] );
assertTrue( all( abs( V - [1;2;3] ) < 1e-14 ))
assertTrue( all( abs( H + 1j*[1;2;3] ) < 1e-14 ))
assertTrue( all( abs( D + [1;2;3] ) < 1e-14 ))

% Change orientation x (swap polarization) and z
[V,H,D] = a.interpolate( pi,0,[],[pi;0;pi] );
assertTrue( all( abs( V + [1;2;3] ) < 1e-14 ))
assertTrue( all( abs( H - 1j*[1;2;3] ) < 1e-14 ))
assertTrue( all( abs( D + [1;2;3] ) < 1e-14 ))

% Different angles for 1st and 3rd element
[V,H,D] = a.interpolate( [0;pi/2],[0;0],[1,3] );
assertTrue( all( abs( V - [1;0.75] ) < 1e-8 ))
assertTrue( all( abs( H + 1j*[1;0.75] ) < 1e-8 ))
assertTrue( all( abs( D + [1;0] ) < 1e-14 ))

% Different angles for 1st and 3rd element
[V,H,D] = a.interpolate( [0;pi/2],[0;0],[1,3],[],pi ); % Force linear interpolation
assertTrue( all( abs( V - [1;0.75] ) < 1e-8 ))
assertTrue( all( abs( H + 1j*[1;0.75] ) < 1e-8 ))
assertTrue( all( abs( D + [1;0] ) < 1e-14 ))

% Different orientations for different angles
orientation = permute( [0 0;0 0;pi/2 -pi/2], [1,3,2] );
[V,~,D] = a.interpolate( [0 pi/2;pi/2 -pi/2],[0 0;0 0],[1,2], orientation);
assertTrue( all(abs(V(2,:)-2)<1e-8))
assertTrue( all(abs(V(1,1,:)-0.25)<1e-8))
assertTrue( all(abs(V(1,2,:))<1e-2))
assertTrue( all(abs(D(:)+[0;2;-1;2])<1e-14))


