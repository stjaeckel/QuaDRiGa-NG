

b = qd_builder('Freespace');
b.simpar.show_progress_bars = 0;
b.simpar.center_frequency = 2.7e9;
b.tx_position = [0;0;10];
b.rx_positions = [100;0;1.5];

bx = copy(b);

b.gen_parameters;
b.add_sdc([50;50;3],-3,[],[],[],[],1);

rt_struct = struct([]);
rt_struct(1,1).tx_pos = b.tx_position;
rt_struct(1,1).rx_pos = b.rx_positions;
rt_struct(1,1).frequency = b.simpar.center_frequency;
rt_struct(1,1).pow = b.gain;
rt_struct(1,1).delay = norm(b.tx_position - b.rx_positions)./qd_simulation_parameters.speed_of_light + b.taus;
rt_struct(1,1).aod = b.AoD;
rt_struct(1,1).eod = b.EoD;
rt_struct(1,1).xprmat = reshape(b.xprmat,2,2,[]);

bx = bx.import_rt_paths( rt_struct );

assertTrue(  all(abs(exp(1j*bx.AoA) - exp(1j*b.AoA)) < 1e-14) )

assertTrue(  all(abs( bx.EoA - b.EoA ) < 1e-14) )

