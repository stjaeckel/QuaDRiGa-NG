function testBuilder_get_pl_dual
%%

b = qd_builder('3GPP_3D_UMa_LOS');
t = qd_track('linear',100);
t.interpolate_positions(0.1);
b.simpar.show_progress_bars = 0;

% Parralel tracks should have same PL
b.rx_track = t.copy;
b.rx_track.initial_position(3)=1;
b.tx_track = t.copy;
b.tx_track.initial_position(1) = 50;
pl = b.get_pl(t,[],b.tx_track);
assertEqual( numel(pl), t.no_snapshots );
assertTrue( all( abs(pl-pl(1)) < 1e-7 ) );
assertEqual( b.dual_mobility, -1 );     % PL doesnt care for dual mobility

% Reversed tracks should have matching SF
b.tx_track = t.copy;
b.tx_track.positions = b.tx_track.positions(:,end:-1:1);

b.gen_parameters;

[ sf,kf ] = get_sf_profile( b, b.rx_track, b.tx_track );

all( abs( sf - sf(end:-1:1) ) < 1e-7 );
all( abs( kf - kf(end:-1:1) ) < 1e-7 );

