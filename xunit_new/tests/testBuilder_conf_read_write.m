function testBuilder_conf_read_write
%%
b = qd_builder('LOSonly');

% Fill scenario table with random variables
scenpar = b.scenpar;
names = fieldnames( scenpar );
for n = 1:numel( names )
    if strcmp( names{n}, 'NumClusters' ) || strcmp( names{n}, 'NumSubPaths' )
        scenpar.(names{n}) = 2 + randi(20);
    elseif strcmp( names{n}, 'GR_enabled' )
        scenpar.(names{n}) = 1;
    elseif strcmp( names{n}, 'GR_epsilon' )  
        scenpar.(names{n}) = rand+1j*rand;
    elseif strcmp( names{n}, 'SubpathMethod' )   
        scenpar.(names{n}) = 'bla';
    else
        scenpar.(names{n}) = rand;
    end
end
b.scenpar = scenpar;

% Check if the XCORR matris is correct
assertEqual( b.lsp_xcorr(1,2) , b.scenpar.ds_kf )
assertEqual( b.lsp_xcorr(1,3) , b.scenpar.ds_sf );
assertEqual( b.lsp_xcorr(1,4) , b.scenpar.asD_ds );
assertEqual( b.lsp_xcorr(1,5) , b.scenpar.asA_ds');
assertEqual( b.lsp_xcorr(1,6) , b.scenpar.esD_ds');
assertEqual( b.lsp_xcorr(1,7) , b.scenpar.esA_ds');
assertEqual( b.lsp_xcorr(1,8) , b.scenpar.xpr_ds');
assertEqual( b.lsp_xcorr(2,3) , b.scenpar.sf_kf');
assertEqual( b.lsp_xcorr(2,4) , b.scenpar.asD_kf');
assertEqual( b.lsp_xcorr(2,5) , b.scenpar.asA_kf');
assertEqual( b.lsp_xcorr(2,6) , b.scenpar.esD_kf');
assertEqual( b.lsp_xcorr(2,7) , b.scenpar.esA_kf');
assertEqual( b.lsp_xcorr(2,8) , b.scenpar.xpr_kf');
assertEqual( b.lsp_xcorr(3,4) , b.scenpar.asD_sf');
assertEqual( b.lsp_xcorr(3,5) , b.scenpar.asA_sf');
assertEqual( b.lsp_xcorr(3,6) , b.scenpar.esD_sf');
assertEqual( b.lsp_xcorr(3,7) , b.scenpar.esA_sf');
assertEqual( b.lsp_xcorr(3,8) , b.scenpar.xpr_sf');
assertEqual( b.lsp_xcorr(4,5) , b.scenpar.asD_asA');
assertEqual( b.lsp_xcorr(4,6) , b.scenpar.esD_asD');
assertEqual( b.lsp_xcorr(4,7) , b.scenpar.esA_asD');
assertEqual( b.lsp_xcorr(4,8) , b.scenpar.xpr_asd');
assertEqual( b.lsp_xcorr(5,6) , b.scenpar.esD_asA');
assertEqual( b.lsp_xcorr(5,7) , b.scenpar.esA_asA');
assertEqual( b.lsp_xcorr(5,8) , b.scenpar.xpr_asa');
assertEqual( b.lsp_xcorr(6,7) , b.scenpar.esD_esA');
assertEqual( b.lsp_xcorr(6,8) , b.scenpar.xpr_esd');
assertEqual( b.lsp_xcorr(7,8) , b.scenpar.xpr_esa');

assertTrue( all( abs( reshape( tril( b.lsp_xcorr )' - triu( b.lsp_xcorr ),[],1) ) < 1e-12 ) );


% Write a conf file
b.write_conf_file('test.conf');

% Read the conf file
c = qd_builder('test');

% Compare results
s = c.scenpar;

for n = 1:numel( names )
    if strcmp( names{n}, 'SubpathMethod' )  
        assertTrue( strcmp( scenpar.(names{n}), s.(names{n})  ) );
    else
        assertTrue( abs( scenpar.(names{n}) - s.(names{n}) ) < 1e-5 );
    end
end

assertTrue( all( abs( b.lsp_xcorr(:) - c.lsp_xcorr(:) ) < 5e-6 ) )

delete('test.conf')

% Check the xcorr-set function
b.lsp_xcorr = eye(8);
assertTrue( all(abs( reshape( b.lsp_xcorr-eye(8),[],1 ) )<1e-12) )

