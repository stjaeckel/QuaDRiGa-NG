function testQF_slerp
%%
x = 0:120:360;
phi = x*pi/180;
xi = 0:22.5:360;

[ phiI, thetaI, pI ]  = qf.slerp ( x, [phi ; phi+pi/2]', 0, xi );

