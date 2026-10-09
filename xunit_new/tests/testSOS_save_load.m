function testSOS_save_load
%%
% Delete all mat files
fn = dir('*.mat');
for n = 1:numel( fn )
    delete(fn(n).name);
end

x = qd_sos;
x.save('sos_test.mat');

y = qd_sos.load( 'sos_test.mat' );

assertTrue( strcmp(y.name,'sos_test.mat') )
assertTrue( strcmp(y.distribution,'Normal') )
assertTrue( y.dist_decorr == 10 )
assertTrue( y.dimensions == 3 )
assertTrue( y.no_coefficients == 300 )
assertTrue( y.dist_decorr - 10 < 1e-5 )
assertTrue( sum( abs(x.sos_freq(:) - y.sos_freq(:)) < 1e-5 ) == 900 )
assertTrue( sum( abs(x.acf(:) - y.acf(:)) < 1e-5 ) == 200 )
assertTrue( sum( abs(x.dist(:) - y.dist(:)) < 1e-5 ) == 200 )

delete('sos_test.mat')
