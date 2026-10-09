function testSOS_default_approx
%%

x = qd_sos('Exp300',[],20);
[Ro,Do] = x.acf_approx;

assertTrue( all( abs(Ro(399,399,:) - 1) < 1e-6 ) );

assertTrue(all(abs(Ro(399,399:598,1) - x.acf)<0.1))
assertTrue(all(abs(Ro(399,399:598,2) - x.acf)<0.1))
assertTrue(all(abs(Ro(399,399:598,3) - x.acf)<0.1))

assertTrue(all(abs(Ro(399:598,399,1)' - x.acf)<0.1))
assertTrue(all(abs(Ro(399:598,399,2)' - x.acf)<0.1))
assertTrue(all(abs(Ro(399:598,399,3)' - x.acf)<0.1))