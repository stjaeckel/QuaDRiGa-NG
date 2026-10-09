function testBuilder_pol_rotation_dual
%%

for reverse = 0 : 1
    reverse = logical( reverse );
    
    %
    ap = qd_arrayant('omni');                                           % Probe antenna
    ap.copy_element(1,2);                                               % Swap polarization
    ap.Fa(:,:,2) = 0;
    ap.Fb(:,:,2) = 1;
    
    at = qd_arrayant('custom',90,45,0);                                 % Test antenna
    gain = at.calc_gain;                                                % Calc gain
    
    % Don't change array grid! It must match the sample grid to avoid interpoaltion artefacts.
    at.set_grid( (-180:22.5:180)*pi/180 , (-90:22.5:90)*pi/180 );
    at.rotate_pattern(90,'z');                                          % Beam faces north
    amp = sqrt(10^(0.1*gain));                                          % Calc amplitude in main direction
    
    b = qd_builder('LOSonly');                                          % LOS builder
    b.simpar.show_progress_bars = 0;                                    % Disable progress bar
    b.simpar.use_absolute_delays = 1;                                   % Use absolute delays
    
    ang = 0:22.5:360;                                                   % Sample angles
    pos = 100 * b.simpar.wavelength * exp(1j*ang*pi/180);
    pos = [ real(pos); imag(pos); zeros(1,numel(ang)) ];                % Probe positions
    
    if reverse
        b.tx_position = pos;                                            % Probe positions
        b.tx_array = ap;                                                % Probe antenna = Tx
        b.rx_positions  = [0;0;0];                                      % Test position
        b.rx_array = at;                                                % Test antenna
    else
        b.rx_positions = pos;                                           % Probe positions
        b.rx_array = ap;                                                % Probe antenna = Rx
        b.tx_position  = [0;0;0];                                       % Test position
        b.tx_array = at;                                                % Test antenna
    end
    
    b.gen_parameters;                                               % Create SSF parameters
    
    % Test 1 : Test the 90° opening in az-direction
    c1 = b.get_channels;                                                % Get channels (regular)
    h1 = cat(2,c1.coeff);
    h1 = reshape( h1,2,[] );
    c2 = b.get_los_channels;                                            % Get channels (LOS)
    h2 = reshape( c2.coeff,2,[] );
    
    assertTrue( all( abs( h1(:)-h2(:)) < 1e-12 ) );                     % Test for identical results
    assertTrue( all( abs( angle(h1(1,:)) ) < 1e-12 ) );                 % Phase mest be zero (distance = 100 lambda)
    assertTrue( all( abs( h1(2,:)) < 1e-12 ) );                         % H-pol must be 0
    assertTrue( all( abs( h1(1,5) - amp ) < 0.05 ) );                   % Gain must match
    assertTrue( all( abs( h1(1,[3,7]) - amp/sqrt(2) ) < 0.05 ) );       % 3dB opening must match
    
    % Test 2 : Rotate antenna py 90° around y-axis and z-axis
    if reverse
        b.rx_array = at.copy;                                           % Assign test antenna
        b.rx_array.rotate_pattern(90,'y');
        b.rx_array.rotate_pattern(90,'z');
    else
        b.tx_array = at.copy;                                           % Assign test antenna
        b.tx_array.rotate_pattern(90,'y');
        b.tx_array.rotate_pattern(90,'z');
    end
    b.check_dual_mobility;                                              % Validate input variables
    
    c1 = b.get_channels;                                                % Get channels (regular)
    h1 = cat(2,c1.coeff);
    h1 = reshape( h1,2,[] );
    c2 = b.get_los_channels;                                            % Get channels (LOS)
    h2 = reshape( c2.coeff,2,[] );
    
    assertTrue( all( abs( h1(:)-h2(:)) < 1e-12 ) );                     % Test for identical results
    assertTrue( all( abs( h1(1,:)) < 1e-12 ) );                         % V-pol must be 0
    assertTrue( all( abs( imag(h1(2,:)) ) < 1e-12 ) );                   % Phase mest be zero or pi (distance = 100 lambda)
    assertTrue( all( abs( h1(2,9) + amp ) < 0.05 ) );                   % Gain must match (negative phase due to pol-rotation)
    assertTrue( all( abs( h1(2,[8,10]) + amp/sqrt(2) ) < 0.05 ) );      % 3dB opening must match
    
    % Test 3 : Same as test 2, but using track instead
    if reverse
        b.rx_array = at.copy;                                           % Assign test antenna
        b.rx_track(1,1).orientation = [ 0 ; -pi/2 ; pi/2 ];             % Orientation (pitch and yaw)
    else
        b.tx_array = at.copy;                                           % Assign test antenna
        b.tx_track(1,1).orientation = [ 0 ; -pi/2 ; pi/2 ];             % Orientation (pitch and yaw)
    end
    b.check_dual_mobility;                                              % Validate input variables
    
    c2 = b.get_channels;                                                % Get channels (regular)
    h2 = cat(2,c2.coeff);
    h2 = reshape( h2,2,[] );
    c3 = b.get_los_channels;                                            % Get channels (LOS)
    h3 = reshape( c3.coeff,2,[] );
    
    assertTrue( all( abs(h2(:)-h3(:)) < 1e-12 ) );                      % Test for identical results
    assertTrue( all( abs(h1(:)-h3(:)) < 1e-12 ) );                      % Test for identical results
    
    % Test 4 : Rotate antenna by 45° around y-axis
    if reverse
        b.rx_array = at.copy;                                           % Assign test antenna
        b.rx_track(1,1).orientation = [0;0;0];                          % Remove track rotation
        b.rx_array.rotate_pattern(45,'y');
    else
        b.tx_array = at.copy;                                           % Assign test antenna
        b.tx_track(1,1).orientation = [0;0;0];                          % Remove track rotation
        b.tx_array.rotate_pattern(45,'y');
    end
    b.check_dual_mobility;                                              % Validate input variables
    
    c1 = b.get_channels;                                                % Get channels (regular)
    h1 = cat(2,c1.coeff);
    h1 = reshape( h1,2,[] );
    c2 = b.get_los_channels;                                            % Get channels (LOS)
    h2 = reshape( c2.coeff,2,[] );
    
    assertTrue( all( abs( h1(:)-h2(:)) < 1e-12 ) );                     % Test for identical results
    assertTrue( all( abs( imag(h1(:)) ) < 1e-12 ) );                    % Phase mest be zero or pi (distance = 100 lambda)
    assertTrue( all( abs( h1(1,5) - amp/sqrt(2) ) < 0.05 ) );           % V-Pol has positive amplitude
    assertTrue( all( abs( h1(2,5) + amp/sqrt(2) ) < 0.05 ) );           % H-Pol has negative amplitude
    
    % Test 5 : Same as test 4, but using track instead
    if reverse
        b.rx_array = at.copy;                                           % Assign test antenna
        b.rx_track(1,1).orientation = [ 0;-pi/4;0 ];                    % Pitch (must be negative)
    else
        b.tx_array = at.copy;                                           % Assign test antenna
        b.tx_track(1,1).orientation = [ 0;-pi/4;0 ];                    % Pitch (must be negative)
    end
    b.check_dual_mobility;                                              % Validate input variables
    
    c2 = b.get_channels;                                                % Get channels (regular)
    h2 = cat(2,c2.coeff);
    h2 = reshape( h2,2,[] );
    c3 = b.get_los_channels;                                            % Get channels (LOS)
    h3 = reshape( c3.coeff,2,[] );
    
    assertTrue( all( abs(h2(:)-h3(:)) < 1e-12 ) );                      % Test for identical results
    assertTrue( all( abs(h1(:)-h3(:)) < 1e-12 ) );                      % Test for identical results
    
    % Test 6 : Rotate antenna by -45° around y-axis
    if reverse
        b.rx_array = at.copy;                                           % Assign test antenna
        b.rx_track(1,1).orientation = [0;0;0];                          % Remove track rotation
        b.rx_array.rotate_pattern(-45,'y');
    else
        b.tx_array = at.copy;                                           % Assign test antenna
        b.tx_track(1,1).orientation = [0;0;0];                          % Remove track rotation
        b.tx_array.rotate_pattern(-45,'y');
    end
    b.check_dual_mobility;                                              % Validate input variables
    
    c1 = b.get_channels;                                                % Get channels (regular)
    h1 = cat(2,c1.coeff);
    h1 = reshape( h1,2,[] );
    c2 = b.get_los_channels;                                            % Get channels (LOS)
    h2 = reshape( c2.coeff,2,[] );
    
    assertTrue( all( abs( h1(:)-h2(:)) < 1e-12 ) );                     % Test for identical results
    assertTrue( all( abs( imag(h1(:)) ) < 1e-12 ) );                    % Phase mest be zero or pi (distance = 100 lambda)
    assertTrue( all( abs( h1(1,5) - amp/sqrt(2) ) < 0.05 ) );           % V-Pol has positive amplitude
    assertTrue( all( abs( h1(2,5) - amp/sqrt(2) ) < 0.05 ) );           % H-Pol has positive amplitude
    
    % Test 7 : Same as test 6, but using track instead
    if reverse
        b.rx_array = at.copy;                                           % Assign test antenna
        b.rx_track(1,1).orientation = [ 0;pi/4;0 ];                     % Rotation around y-axis
    else
        b.tx_array = at.copy;                                           % Assign test antenna
        b.tx_track(1,1).orientation = [ 0;pi/4;0 ];                     % Rotation around y-axis
    end
    b.check_dual_mobility;                                              % Validate input variables
    
    c2 = b.get_channels;                                                % Get channels (regular)
    h2 = cat(2,c2.coeff);
    h2 = reshape( h2,2,[] );
    c3 = b.get_los_channels;                                            % Get channels (LOS)
    h3 = reshape( c3.coeff,2,[] );
    
    assertTrue( all( abs(h2(:)-h3(:)) < 1e-12 ) );                      % Test for identical results
    assertTrue( all( abs(h1(:)-h3(:)) < 1e-12 ) );                      % Test for identical results
    
    % Test 8 : Rotate antenna by -45° around x-axis
    if reverse
        b.rx_array = at.copy;                                           % Assign test antenna
        b.rx_track(1,1).orientation = [0;0;0];                          % Remove track rotation
        b.rx_array.rotate_pattern(-90,'z');                             % Align with broadside
        b.rx_array.rotate_pattern(45,'x');
    else
        b.tx_array = at.copy;                                           % Assign test antenna
        b.tx_track(1,1).orientation = [0;0;0];                          % Remove track rotation
        b.tx_array.rotate_pattern(-90,'z');                             % Align with broadside
        b.tx_array.rotate_pattern(45,'x');
    end
    b.check_dual_mobility;                                              % Validate input variables
    
    c1 = b.get_channels;                                                % Get channels (regular)
    h1 = cat(2,c1.coeff);
    h1 = reshape( h1,2,[] );
    c2 = b.get_los_channels;                                            % Get channels (LOS)
    h2 = reshape( c2.coeff,2,[] );
    
    assertTrue( all( abs( h1(:)-h2(:)) < 1e-12 ) );                     % Test for identical results
    assertTrue( all( abs( imag(h1(:)) ) < 1e-12 ) );                    % Phase mest be zero or pi (distance = 100 lambda)
    assertTrue( all( abs( h1(1,[1,17]) - amp/sqrt(2) ) < 0.05 ) );      % V-Pol has positive amplitude
    assertTrue( all( abs( h1(2,[1,17]) + amp/sqrt(2) ) < 0.05 ) );      % H-Pol has negative amplitude
    
    % Test 9 : Same as test 8, but using track instead
    if reverse
        b.rx_array = at.copy;                                           % Assign test antenna
        b.rx_array.rotate_pattern(-90,'z');                             % Align with broadside
        b.rx_track(1,1).orientation = [ pi/4; 0 ; 0 ];                  % Roll
    else
        b.tx_array = at.copy;                                           % Assign test antenna
        b.tx_array.rotate_pattern(-90,'z');                             % Align with broadside
        b.tx_track(1,1).orientation = [ pi/4; 0 ; 0 ];                  % Roll
    end
    b.check_dual_mobility;                                              % Validate input variables
    
    c2 = b.get_channels;                                                % Get channels (regular)
    h2 = cat(2,c2.coeff);
    h2 = reshape( h2,2,[] );
    c3 = b.get_los_channels;                                            % Get channels (LOS)
    h3 = reshape( c3.coeff,2,[] );
    
    assertTrue( all( abs(h2(:)-h3(:)) < 1e-12 ) );                      % Test for identical results
    assertTrue( all( abs(h1(:)-h3(:)) < 1e-12 ) );                      % Test for identical results
    
end
