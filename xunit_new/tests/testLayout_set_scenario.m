function testLayout_set_scenario
%%

s = qd_simulation_parameters;
s.center_frequency = [ 2e9, 28e9 ];                 % Two frequencies
s.show_progress_bars = 0;

l = qd_layout( s );
l.no_tx = 2;                                        % 2 BS
l.tx_position(2,2) = 100;                           % 100 m ISD

l.no_rx = 4;                                        % 10 Rx
l.randomize_rx_positions( 500,1.5,50,0 );           % 500 m radius, random heights

t = qd_track('street',50,0);                        % New track with MT mobility
t.name = l.rx_name{1};                              % Keep Rx name
t.initial_position = [0,0,1.5]';                    % Start point under BS1

scen = {'3GPP_3D_UMi','3GPP_3D_UMa','3GPP_38.901_UMi','3GPP_38.901_UMa','3GPP_38.901_RMa',...
    '3GPP_38.901_Indoor_Mixed_Office','3GPP_38.901_Indoor_Open_Office','mmMAGIC_UMi','mmMAGIC_Indoor','LOSonly'};

for n = 1:numel( scen )
    l.set_scenario(scen{n},[],[],0.8);
end
