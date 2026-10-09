function test_all_config_files
%%

set_rand_state( 1 );

show = 1;

list = qd_builder.supported_scenarios(0);

l = qd_layout;
l.simpar.show_progress_bars = 0;
l.simpar.center_frequency = 2.6e9;
l.randomize_rx_positions(100,0,0,0.03);
l.tx_position(3) = 25;

l = l(1,1);

for n = 1 : numel(list)
    l.track(1,1).scenario = list{n};
    
    if show
        fprintf( [num2str(n),' - ',list{n},' ... '] );
    end
        
    l.track(1,1).par = [];
    [c,p] = l.get_channels;
    p = p(1,1);
    
    assertEqual( list{n} , p.scenario );
    
    n_clusters = p.scenpar.NumClusters;
    
    if p.scenpar.GR_enabled == 1; % Ground reflection
        nLOS = 2;
    else
        nLOS = 1;
    end
    
    switch p.scenpar.SubpathMethod
        case 'legacy'
            if p.scenpar.PerClusterDS > 0       % Sub-Path splitting
                if p.scenpar.SC_lambda > 0      % SSF spatial consistency
                    n_clusters = max( (n_clusters-nLOS),0) *3 + nLOS;
                else
                    nNLOS = n_clusters-nLOS;
                    nSPLIT = min( nNLOS,2 );
                    n_clusters = nLOS + (nNLOS-nSPLIT) + nSPLIT*3;
                end
            end
        case 'mmMAGIC'
            n_clusters = max( (n_clusters-nLOS),0) *p.scenpar.NumSubPaths + nLOS;
        otherwise
            error('wrong Subpath Method');
            
    end
    assertEqual( n_clusters , c.no_path );
    
    tmp = reshape( c.coeff,[],1 );
    if ~(   all( isnumeric( tmp ) ) &&...
            all( ~isnan( tmp ) )  &&...
            all( abs(tmp) > 1e-60 )  )
        error( [ list{n},' : Incorrect coefficients.' ] );
    end
    
    tmp = reshape( c.delay,[],1 );
    if ~(   all( isnumeric( tmp ) ) &&...
            all( ~isnan( tmp ) ) )
        error( [ list{n},' : Incorrect delays.' ] );
    end
    
    
    if show
        fprintf( ['OK\n'] );
    end
end

