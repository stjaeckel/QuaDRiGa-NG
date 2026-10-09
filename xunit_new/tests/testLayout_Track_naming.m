function testLayout_Track_naming
%%

l = qd_layout;
l.simpar.show_progress_bars = 0;
l.no_rx = 2;
assertFalse( strcmp( l.rx_name{1} , l.rx_name{2} ) );
l.randomize_rx_positions;
l.set_scenario('LOSonly');

% This should not work, but it does because we don't chnage anything in layout, only track handles
try
    % Throws error in octave but still assigns name
    % Throws no error in matlab
    l.rx_track(1,2).name = l.rx_track(1,1).name;
catch err
    assertEqual( err.identifier, 'QuaDRiGa:qd_layout:WrongInput' );
end
    
assertFalse( l.has_unique_track_names );

try % Check if error is trown if we try to generate channels
    c = l.get_channels;
catch err
    assertEqual( err.identifier, 'QuaDRiGa:qd_layout:has_unique_track_names' );
end


