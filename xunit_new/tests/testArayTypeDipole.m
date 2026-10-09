function testArayTypeDipole
% Generate dipole

a = qd_arrayant('dipole');

elevation_grid = (-90:90)*pi/180;
azimuth_grid = (-180:180)*pi/180;

% short dipole
[~, theta_grid] = meshgrid(azimuth_grid, elevation_grid);

E_theta = cos((1 - 1e-6) * theta_grid);
E_phi = zeros(size(E_theta));
% calculate radiation power pattern
P = E_theta.^2 + E_phi.^2;
% normalize by max value
P_max = max(max(P));
P = P ./ P_max;
% the gain of a short dipole is 1.76 dBi
gain_lin = sum(sum((cos(theta_grid))))/sum(sum(P.*cos(theta_grid)));
% gain_dbi = 10*log10(gain_lin)

E_theta = E_theta .* sqrt(gain_lin./P_max);

assertTrue( all( abs( a.Fa(:) - E_theta(:) ) <1e-6  ) );
assertTrue(  all( abs( a.Fb(:) ) < 1e-6  ) );