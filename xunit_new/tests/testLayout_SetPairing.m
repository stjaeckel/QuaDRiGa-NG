function testLayout_SetPairing
%%

l = qd_layout;
l.simpar.center_frequency = 2e9;
l.simpar.show_progress_bars = 0;
l.no_tx = 2;
l.tx_position(2,:) = [0,200];
l.no_rx = 50;
l.randomize_rx_positions(100,1.5,1.5,0);
l.set_scenario('Freespace');

% Set pairing to only one link
l.pairing = [1;1];
assertTrue( l.no_links == 1 );

% Reset pairing and check if all links are active
l.set_pairing;      
assertEqual( l.pairing, [ ones(1,l.no_rx) , 2*ones(1,l.no_rx) ; 1:l.no_rx , 1:l.no_rx ] );

% The freespace PL formula
d3d = sqrt(sum(( l.rx_position - l.tx_position(:,1) * ones(1,l.no_rx) ).^2));
dt  = mean(d3d);

% Use poweer level to inlcude only users below the average power
P_tres = 20 * log10(dt) + 32.45 + 20*log10(2);
l.set_pairing( 'power', -P_tres, [],[],0 );

X = find( d3d<dt );
assertEqual( l.pairing, [ ones(1,numel(X)) ; X] );


warning('off','QuaDRiGa:qd_layout:init_builder:no_tx')
c = l.get_channels( [],0 );

% There should only be as many channels as there are entries in X
assertEqual( numel(c), numel(X) ); % One additional user from BS2

for n = 2:numel(X)
    Pc = 10*log10( abs(c(1,n).coeff).^2 );
    Pl = 20 * log10(d3d(X(n))) + 32.45 + 20*log10(2);
    assertTrue( abs( Pc + Pl ) < 0.01 );
end

