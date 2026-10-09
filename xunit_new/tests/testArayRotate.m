function testArayRotate
% Rotte antennas around axis

a = qd_arrayant('dipole');
fp = a.Fa;
a.set_grid( (-180:5:180)*pi/180 , (-90:5:90)*pi/180 );

a.rotate_pattern(-45,'y');
a.rotate_pattern(180,'z');
a.rotate_pattern(90,'x');
a.rotate_pattern(45,'z');
a.rotate_pattern(-90,'y');

a.set_grid( (-180:1:180)*pi/180 , (-90:1:90)*pi/180 );
fpo = a.Fa;

%imagesc(fp-fpo);colorbar

assertTrue( all(abs( fp(:)-fpo(:) ) < 0.01) );