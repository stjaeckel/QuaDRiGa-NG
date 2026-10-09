function testArray_ImportPattern
%%

a = qd_arrayant('patch');
a.append_array( qd_arrayant('xpol') );

a.set_grid((-180:3:180)*pi/180,(-90:3:90)*pi/180)

ii = [61:121,2:60];
azimuth_grid = (0:3:359)*pi/180;
elevation_grid = a.elevation_grid;
fVi = a.Fa(:,ii,:);
fHi = a.Fb(:,ii,:);

b = qd_arrayant.import_pattern( fVi, fHi , azimuth_grid , elevation_grid );

assertTrue( all( abs( a.Fa(:) - b.Fa(:) ) < 1e-14 ) );
assertTrue( all( abs( a.Fb(:) - b.Fb(:) ) < 1e-14 ) );


end

