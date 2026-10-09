function set_rand_state( seed )
%SET_RAND_STATE Helper function to set the rand state on MATLAB and Octave

if ~exist( 'seed','var' )
    seed = 1;
end

v = version;
if ~isempty( strfind( v,'R20' ) )
    % MATLAB
    RandStream.setGlobalStream(RandStream('mt19937ar','seed',seed));
else 
    % Octave
    rand('seed',seed);
end

end

