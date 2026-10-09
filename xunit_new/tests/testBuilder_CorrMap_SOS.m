function testBuilder_CorrMap_SOS
%%
set_rand_state( 1 );

b = qd_builder('BERLIN_UMa_LOS');
b.simpar.show_progress_bars = 0;
b.simpar.use_3GPP_baseline = 1;
b.tx_position = [0;0;25];
b.rx_positions = [0;0;1.5];

% We want to test the exponential decay of the correlation at these
% distances
vv = [4,5,10,15,20,35,50];

% Set the correlation distances (the is for the map resolution)
b.scenpar.DS_lambda = vv(1);
b.scenpar.KF_lambda = vv(2);
b.scenpar.SF_lambda = vv(3);
b.scenpar.AS_D_lambda = vv(4);
b.scenpar.AS_A_lambda = vv(5);
b.scenpar.ES_D_lambda = vv(6);
b.scenpar.ES_A_lambda = vv(7);

% Set some values for the cross-polarization
M = eye(8);
M(1,2) = 0.3;
M(1,5) = 0.1;
M(2,4) = -0.4;
b.lsp_xcorr = M;

b.init_sos;

% Start Values
st = [1,50,100,150,200];

runs = 15;

cx = zeros(numel(st),numel(vv),runs);
cy = cx;
crc = zeros(5,numel(st),runs);

for o=1:runs
    b.init_sos(1);
    for n = 1:numel(st)
        val_x = b.get_lsp_map( (-20:20)*20 , [st(n),st(n)+vv] ) ;
        val_y = b.get_lsp_map( [st(n),st(n)+vv] , (-20:20)*20  );
        for m = 1:7
            tmp = qf.xcorrcoeff( val_x(1,:,1,m)', val_x(2:end,:,1,m)' );
            cx(n,m,o) = tmp(m);
            tmp = qf.xcorrcoeff( val_y(:,1,1,m) , val_y(:,2:end,1,m) );
            cy(n,m,o) = tmp(m);
        end
        
        T1 = [ val_x(:,:,:,1)' , val_y(:,:,:,1) ];
        T2 = [ val_x(:,:,:,2)' , val_y(:,:,:,2) ];
        T3 = [ val_x(:,:,:,3)' , val_y(:,:,:,3) ];
        T4 = [ val_x(:,:,:,4)' , val_y(:,:,:,4) ];
        T5 = [ val_x(:,:,:,5)' , val_y(:,:,:,5) ];
               
        crc(1,n,o) = qf.xcorrcoeff( T1(:), T2(:) );
        crc(2,n,o) = qf.xcorrcoeff( T1(:), T3(:) );
        crc(3,n,o) = qf.xcorrcoeff( T1(:), T5(:) );
        crc(4,n,o) = qf.xcorrcoeff( T2(:), T4(:) );
        crc(5,n,o) = qf.xcorrcoeff( T2(:), T5(:) );
    end
end

Cx  = mean(mean(cx,3),1);
Cy  = mean(mean(cy,3),1);
Crc = mean( mean( crc , 2 ), 3 ).';

% Check the outputs
assertTrue( all( abs( Cx-exp(-1) )<0.2 ) );
assertTrue( all( abs( Cy-exp(-1) )<0.2 ) );
assertTrue( all( abs(Crc - [0.3,0,0.1,-0.4,0]) < 0.1 ) );