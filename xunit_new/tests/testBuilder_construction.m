function testBuilder_construction
%%
% Delete all conf files
fn = dir('testXYZ.conf');
for n = 1:numel( fn )
    delete(fn(n).name);
end
fid = fopen('testXYZ.conf','w');
fprintf(fid,'NumClusters = 66\n');
fclose(fid);

b = qd_builder;
b = qd_builder([]);

b = qd_builder('Freespace');
assertEqual( b.scenario,'Freespace' );
assertEqual( b.scenpar.NumClusters,1 );

b = qd_builder('Ul');
assertEqual( b.scenpar.NumClusters,15 );

% Advanced get-function
b.lsp_xcorr_chk;
b.lsp_vals;
b.lsp_xcorr;

% Set funtion
b.scenario = 'Un';
assertEqual( b.scenpar.NumClusters,25 );

b.scenario = 'testXYZ';
assertEqual( b.scenpar.NumClusters,66 );

tmp = b.scenpar;
tmp.NumClusters = 55;
b.scenpar = tmp;
assertEqual( b.scenpar.NumClusters,55 );
assertEqual( b.scenario,'Custom' );

delete('testXYZ.conf');