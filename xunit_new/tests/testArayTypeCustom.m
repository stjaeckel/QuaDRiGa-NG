function testArayTypeCustom
% Generate custom

ar = qd_arrayant('custom',89.99,89.99,0.1);

elevation = (-90:90)*pi/180;
azimuth = (-180:180)*pi/180;

% short dipole
[~, theta_grid] = meshgrid(azimuth, elevation);

C = exp(-1.314630546162659*azimuth.^2);
D = cos(elevation).^2;

P = zeros(181,361);
for a = 1:181
    for b = 1:361
        P(a,b) = D(a) * C(b);
    end
end
P = 0.1 + (1-0.1)*P;

% normalize by max value
P_max = max(max(P));
P = P ./ P_max;
% the gain of a half-wave dipole is 2.15 dBi
gain_lin = sum(sum((cos(theta_grid))))/sum(sum(P.*cos(theta_grid)));
% gain_dbi = 10*log10(gain_lin)

E_theta = sqrt(P .* gain_lin);

assertTrue( all( abs( ar.Fa(:) - E_theta(:) ) < 1e-2 ) )
assertTrue(  all( abs( ar.Fb(:) ) < 1e-6  ) );
