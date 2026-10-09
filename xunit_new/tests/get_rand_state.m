function state = get_rand_state

v = version;
if ~isempty( strfind( v,'R20' ) )
    % MATLAB
    stream = RandStream.getGlobalStream;
    state = stream.State;
else 
    % Octave
    state = rand('state');
end

end

