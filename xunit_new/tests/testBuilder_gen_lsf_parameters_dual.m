function testBuilder_gen_lsf_parameters_dual
%%
% If tx and rx positions are swapped, LSPs should be identical!

b = qd_builder('3GPP_3D_UMa_LOS');
b.simpar.show_progress_bars = 0;

b.lsp_xcorr = eye(8);       % Disable cross-correlation
b.scenpar.AS_D_lambda = b.scenpar.AS_A_lambda;
b.scenpar.ES_D_mu_A = 0;

b.rx_positions = [0,10,10;0,50,0;200,200,1]';
b.tx_position = b.rx_positions(:,[2 1 3]);
gen_parameters(b);

assertTrue( b.dual_mobility == true );

% DS, KF, SF and XPR should be identical
assertTrue( abs( b.ds(1) - b.ds(2) ) < 1e-12 )
assertTrue( abs( b.kf(1) - b.kf(2) ) < 1e-7 )
assertTrue( abs( b.sf(1) - b.sf(2) ) < 1e-7 )
assertTrue( abs( b.xpr(1) - b.xpr(2) ) < 1e-7 )

% ASD, ASA, ESD, ESA should be different
assertTrue( abs( b.asD(1) - b.asD(2) ) > 1e-7 )
assertTrue( abs( b.asA(1) - b.asA(2) ) > 1e-7 )
assertTrue( abs( b.esD(1) - b.esD(2) ) > 1e-7 )
assertTrue( abs( b.esA(1) - b.esA(2) ) > 1e-7 )

% ASD shoud become ASA 
asD = (log10(b.asD)-b.scenpar.AS_D_mu)./b.scenpar.AS_D_sigma;
asA = (log10(b.asA)-b.scenpar.AS_A_mu)./b.scenpar.AS_A_sigma;
assertTrue( abs( asD(1) - asA(2) ) < 1e-7 )
assertTrue( abs( asD(2) - asA(1) ) < 1e-7 )
assertTrue( abs( asD(3) - asA(3) ) < 1e-5 )

% ESD shoud become ESA 
esD = (log10(b.esD)-b.scenpar.ES_D_mu)./b.scenpar.ES_D_sigma;
esA = (log10(b.esA)-b.scenpar.ES_A_mu)./b.scenpar.ES_A_sigma;
assertTrue( abs( esD(1) - esA(2) ) < 1e-7 )
assertTrue( abs( esD(2) - esA(1) ) < 1e-7 )
assertTrue( abs( esD(3) - esA(3) ) < 1e-5 )

% Different positions should be uncorrelated
assertTrue( abs( b.ds(1) - b.ds(3) ) > 1e-12 )
assertTrue( abs( b.kf(1) - b.kf(3) ) > 1e-7 )
assertTrue( abs( b.sf(1) - b.sf(3) ) > 1e-7 )
assertTrue( abs( b.asD(1) - b.asD(3) ) > 1e-7 )
assertTrue( abs( b.asA(1) - b.asA(3) ) > 1e-7 )
assertTrue( abs( b.esD(1) - b.esD(3) ) > 1e-7 )
assertTrue( abs( b.esA(1) - b.esA(3) ) > 1e-7 )
assertTrue( abs( b.xpr(1) - b.xpr(3) ) > 1e-7 )

