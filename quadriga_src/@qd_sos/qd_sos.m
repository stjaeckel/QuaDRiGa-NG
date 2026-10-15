classdef qd_sos < handle
%QD_SOS QuaDRiGa sum-of-sinuoids autocorrelation model
%
% DESCRIPTION
%
% This class implements the sum-of-sinusoids method for generating spatially correlated random
% variables. An autocorrelation function (ACF), i.e. a description of the correlation vs. distance,
% needs to be provided. This function is approximated by a fourier series. The coefficients of the
% series can the be used to generate spatially correlated random variables. The computations are
% done by Quadriga-Lib (quadriga_lib.sos_* functions).
%
% Calling qd_sos without options creates a new SOS object with default settings. The default ACF
% is defined by a combination of the Gauss profile (ranging from 0 to the decorrelation distance)
% and an exponential decaying profile (ranging from the decorrelation distance to a maximum
% distance of 5 time the decorrelation distance. It is approximated by 256 sinusoids in 3D
% coordinates. There are 3 ACF types that can be selected by acf_type. The 3D SOS frequencies are
% precomputed for 128, 256, 512, and 1024 sinusoids.
%
%    * Exponential ACF (Exp128, Exp256, Exp512, Exp1024 )
%    * Gaussian ACF (Gauss128, Gauss256, Gauss512, Gauss1024 )
%    * Combined Gaussian and Exponential ACF (Comb128, Comb256, Comb512, Comb1024)
%
%    The names of previous versions (150, 300, 500 and 1000 sinusoids) are mapped to 128, 256, 512
%    and 1024 sinusoids.
%
%    The distribution function can be either ('Normal' or 'Uniform').
%
% Input (constructor):
%   acf_type
%   String describing the shape of the autocorrelation function and the number of sinusoids,
%   Default: 'Comb256'
%
%   distribution
%   Distribution of the random variables ('Normal' or 'Uniform'), Default: 'Normal'
%
%   dist_decorr
%   Decorrelation distance in [m], Default: 10 m
%
%
% QuaDRiGa Copyright (C) 2011-2019
% Fraunhofer-Gesellschaft zur Foerderung der angewandten Forschung e.V. acting on behalf of its
% Fraunhofer Heinrich Hertz Institute, Einsteinufer 37, 10587 Berlin, Germany
% All rights reserved.
%
% e-mail: quadriga@hhi.fraunhofer.de
%
% This file is part of QuaDRiGa.
%
% The Quadriga software is provided by Fraunhofer on behalf of the copyright holders and
% contributors "AS IS" and WITHOUT ANY EXPRESS OR IMPLIED WARRANTIES, including but not limited to
% the implied warranties of merchantability and fitness for a particular purpose.
%
% You can redistribute it and/or modify QuaDRiGa under the terms of the Software License for
% The QuaDRiGa Channel Model. You should have received a copy of the Software License for The
% QuaDRiGa Channel Model along with QuaDRiGa. If not, see <http://quadriga-channel-model.de/>.

properties
    % Name of the preset
    name = 'Custom';

    % Distribution of random variables (Normal or Uniform)
    distribution = 'Normal';
end

properties(Dependent)
    % Decorrelation distance in [m]
    dist_decorr

    % The number of dimensions (1, 2 or 3)
    dimensions

    % Number of sampling frequencies
    no_coefficients
end

properties
    % Vector of sample points for the ACF in [m]
    dist = [];

    % Desired aotocorrelation function
    acf = [];

    % Sinusoid coefficients (frequencies), [ no_coefficients x dimensions ]
    sos_freq = [];

    % Phase offsets, [ no_coefficients x 2 ]
    %   Column 1 is used for the first device, column 2 for the second device in dual-mobility
    %   scenarios. If both columns are identical, both devices are on the same map.
    sos_phase = [];
end

properties(Access=private)
    Pdist_decorr = [];
end

properties(Hidden)
    OctEq = false; % For qf.eq_octave
end

methods
    % Constructor
    function h_sos = qd_sos( acf_type, distribution, dist_decorr )
        if nargin == 0
            acf_type = 'Comb256';
        elseif isempty( acf_type )
            return
        end

        if ~ischar( acf_type )
            error('QuaDRiGa:qd_sos:type','ACF type must be a string.');
        end

        % Map the names of previous versions to the presets of Quadriga-Lib
        tok = regexp( acf_type, '^(Exp|Gauss|Comb)(150|300|500|1000)$', 'tokens', 'once' );
        if ~isempty( tok )
            new_no = {'128','256','512','1024'};
            acf_type = [ tok{1}, new_no{ strcmp( tok{2}, {'150','300','500','1000'} ) } ];
        end

        % The builder creates thousands of generators, so the preset tables are read only once
        persistent presets
        if isempty( presets )
            presets = struct();
        end
        if ~isfield( presets, acf_type )
            [ freq, acf, dist ] = quadriga_lib.sos_preset( acf_type, 10 );
            presets.( acf_type ) = { freq, acf, dist };
        end
        h_sos.sos_freq = presets.( acf_type ){1};
        h_sos.acf = presets.( acf_type ){2};
        h_sos.dist = presets.( acf_type ){3};
        h_sos.name = acf_type;
        h_sos.Pdist_decorr = 10;
        h_sos.init;

        if exist( 'distribution','var' ) && ~isempty( distribution )
            h_sos.distribution = distribution;
        end

        if exist( 'dist_decorr','var' ) && ~isempty( dist_decorr ) && any( dist_decorr ~= 10 )
            h_sos.dist_decorr = dist_decorr;
        end
    end

    % Get functions
    function out = get.dist_decorr(h_sos)
        out = h_sos.Pdist_decorr;
    end
    function out = get.dimensions(h_sos)
        out = size( h_sos.sos_freq,2 );
    end
    function out = get.no_coefficients(h_sos)
        out = size( h_sos.sos_freq,1 );
    end

    % Set functions
    function set.dist_decorr(h_sos,val)
        if ~( ( all(size(val) == [1 1]) || all(size(val) == [1 2]) ) && isnumeric(val) && isreal(val) && all(val > 0) )
            error('QuaDRiGa:qd_sos:wrongInputValue','??? "decorr_dist" must be scalar and > 0')
        end
        if isempty( h_sos.dist )
            error('QuaDRiGa:qd_sos:noACFgiven','??? No ACF specified')
        end

        if numel( val ) == 2 && val(1) == val(2)
            val = val(1);
        end

        % Adjust distance vector
        maxD_old = h_sos.dist(end);
        h_sos.dist = val(1) ./ h_sos.Pdist_decorr(1) .* h_sos.dist;

        % Adjust the frequencies
        maxD_new = max(h_sos.dist);
        h_sos.sos_freq = h_sos.sos_freq * maxD_old / maxD_new;

        % Save new decorr distance
        h_sos.Pdist_decorr = val;
    end

    function set.distribution(h_sos,value)
        supported_types = {'Normal','Uniform'};
        if ~( ischar(value) && any( strcmpi(value,supported_types)) )
            error('QuaDRiGa:qd_sos:wrongDistributionType',...
                'Distribution not supported; supported types are: Normal, Uniform');
        end
        h_sos.distribution = value;
    end

    function init( h_sos, use_same )
        %INIT Initializes the random phases
        %
        %   The phases are drawn by quadriga_lib.sos_init_phases. Its seed is taken from the MATLAB
        %   random number generator, so the phases follow the state set by "rng".
        %
        % Input:
        %   use_same
        %   If set to 1 (default), identical phases are used for the dual-mobility option, assuming
        %   that devices are on the same radio-map. If set to 0, different values are used.

        if ~exist( 'use_same','var' )
            use_same = true;
        end
        for n = 1 : numel( h_sos )
            seed = floor( rand * 2^52 );
            if use_same
                phase = quadriga_lib.sos_init_phases( h_sos(n).no_coefficients, 1, seed );
                h_sos(n).sos_phase = [ phase, phase ];
            else
                h_sos(n).sos_phase = quadriga_lib.sos_init_phases( h_sos(n).no_coefficients, 2, seed );
            end
        end
    end

    function s = val( h_sos, ca, cb )
        %VAL Returns correlated values at given coordinates
        %
        % Input:
        %   ca
        %   Coordinates for the first mobile device in [m] given as [3 x N] matrix. The rows
        %   correspond to the x,y and z coordinate.
        %
        %   cb
        %   Coordinates for the corresponding second mobile device in [m] given as [3 x N] matrix.
        %   The rows correspond to the x,y and z coordinate. This variable must either be empty or
        %   have the same size as "ca".
        %
        % Output:
        %   s
        %   Random spatially correlated values for each object in the array, [ numel(h_sos) x N ]

        if ~exist( 'cb','var' ) || isempty( cb )
            cb = [];
        elseif any( size(cb) ~= size(ca) )
            error('QuaDRiGa:qd_sos:val','Size of "cb" mut match the size of "ca".' );
        end

        n_sos = numel( h_sos );
        no_val = size( ca,2 );
        s = zeros( n_sos, no_val );
        if no_val == 0
            return
        end

        [ mode, uniform, scale_b ] = eval_mode( h_sos, ~isempty( cb ) );
        if n_sos == 1
            s = sos_eval( h_sos, mode, uniform, scale_b, ca, cb ).';
            return
        end
        
        % Objects that share the frequencies, the ACF, the distribution and the evaluation mode are
        % evaluated in one call. This is the case for the generator arrays of the builder.
        same_tables = false( 1, n_sos );
        for n = 1 : n_sos
            same_tables(n) = isequal( h_sos(n).sos_freq, h_sos(1).sos_freq ) && ...
                isequal( h_sos(n).acf, h_sos(1).acf ) && isequal( h_sos(n).dist, h_sos(1).dist );
        end
        done = false( 1, n_sos );
        for n = 1 : n_sos
            if done(n)
                continue
            end
            if same_tables(n)
                ii = find( same_tables & ~done & mode == mode(n) & uniform == uniform(n) & scale_b == scale_b(n) );
            else
                ii = n;
            end
            s( ii,: ) = sos_eval( h_sos(ii), mode(n), uniform(n), scale_b(n), ca, cb ).';
            done( ii ) = true;
        end
    end

    function remove_bias( h_sos, ca, cb )
        %REMOVE_BIAS Removes the bias from the SOS generator for the given set of positions
        %
        % Description:
        %   One realisation of a spatially correlated random field has a non-zero mean over a finite
        %   set of positions. This method adds an offset to the phases of the SOS generator that
        %   forces the mean to be zero for a set of given positions. The offset is calculated in
        %   closed form by quadriga_lib.sos_remove_bias. Note that this removes variance from the
        %   random variables when the positions lie within a few decorrelation distances.
        %
        % Input:
        %   ca
        %   Coordinates for the first mobile device in [m] given as [3 x N] matrix. The rows
        %   correspond to the x,y and z coordinate.
        %
        %   cb
        %   Coordinates for the corresponding second mobile device in [m] given as [3 x N] matrix.
        %   The rows correspond to the x,y and z coordinate. This variable must either be empty or
        %   have the same size as "ca".

        if ~exist( 'cb','var' ) || isempty( cb )
            cb = [];
        elseif any( size(cb) ~= size(ca) )
            error('QuaDRiGa:qd_sos:remove_bias','Size of "cb" mut match the size of "ca".' );
        end
        if size( ca,2 ) == 0
            return
        end

        for n = 1 : numel( h_sos )
            h = h_sos(n);
            ph = h.sos_phase;
            use_same = all( ph(:,1) == ph(:,2) );
            if numel( h.Pdist_decorr ) == 2
                scale_b = h.Pdist_decorr(1) / h.Pdist_decorr(2);
            else
                scale_b = 1;
            end

            if use_same && scale_b == 1 && isempty( cb )
                phase = quadriga_lib.sos_remove_bias( h.sos_freq, ph(:,1), ca.' );
                h.sos_phase = [ phase, phase ];
            elseif use_same && scale_b == 1
                phase = quadriga_lib.sos_remove_bias( h.sos_freq, ph(:,1), ca.', cb.', [], h.acf, h.dist );
                h.sos_phase = [ phase, phase ];
            else
                % Each device on its own map; without second positions, the second phase column is
                % adjusted for the positions of the first device
                if isempty( cb )
                    pos_b = ca.';
                else
                    pos_b = cb.';
                end
                [ phase, phase_b ] = quadriga_lib.sos_remove_bias( h.sos_freq, ph(:,1), ca.', pos_b, ph(:,2), [], [], scale_b );
                h.sos_phase = [ phase, phase_b ];
            end
        end
    end

    function val = acfi( h_sos, dist )
        %ACFI Linear interpolation of the ACF
        %
        % Input:
        %   dist
        %   Vector containing the distances in [m] for which the ACF should be interpolated.
        %
        % Output:
        %   val
        %   Interpolated ACF at the given distances

        if numel( h_sos ) > 1
            error('QuaDRiGa:qd_sos:acfi','acfi not definded for object arrays.');
        end
        val = quadriga_lib.sos_acfi( h_sos.acf, h_sos.dist, dist(:) );
        val = reshape( val, size(dist) );
    end

    function out = copy( obj, in )
        %COPY Creates a copy of the SOS object
        %
        % While the standard copy command creates new physical objects for each element of the
        % object array (in case obj is an array of object handles), copy checks whether there are
        % object handles pointing to the same object and keeps this information.

        if ~exist( 'in','var' )
            sic = size( obj );
            prc = false( sic );
            out = obj;
            for n = 1 : prod( sic )
                if ~prc( n )
                    [ i1,i2,i3,i4 ] = qf.qind2sub( sic, n );
                    out( i1,i2,i3,i4 ) = qd_sos( [] );
                    copy( out(i1,i2,i3,i4), obj(i1,i2,i3,i4) );
                    prc( i1,i2,i3,i4 ) = true;
                    m = qf.eqo( obj(i1,i2,i3,i4), obj );
                    m(i1,i2,i3,i4) = false;
                    if any( m(:) )
                        out( m ) = out( i1,i2,i3,i4 );
                        prc( m ) = true;
                    end
                end
            end
            if numel( obj ) == 1
                out = out(1,1);
            end
        else
            prop = {'name','distribution','Pdist_decorr','dist','acf','sos_freq','sos_phase'};
            out = [];
            for n = 1 : numel(prop)
                obj.( prop{n} ) = in.( prop{n} );
            end
        end
    end
end

methods(Static)
    function [ s, h_sos ] = randn( dist_decorr, ca, cb, acf_type )
        %RANDN Generates spatially correlated normal distributed random numbers
        %
        % Input:
        %   dist_decorr
        %   Vector of decorrelation distances  [1 x M] or [ M x 1 ]
        %
        %   ca
        %   Coordinates for the first mobile device in [m] given as [3 x N] matrix.
        %
        %   cb
        %   Coordinates for the corresponding second mobile device in [m] given as [3 x N] matrix.
        %   This variable must either be empty or have the same size as "ca".
        %
        %   acf_type
        %   String describing the shape of the autocorrelation function and the number of
        %   sinusoids, Default: 'Comb256'
        %
        % Output:
        %   s
        %   Random spatially correlated numbers [ M x N ]
        %
        %   h_sos
        %   A handle to the used qd_sos object

        if ~exist( 'cb','var' )
            cb = [];
        end
        if ~exist( 'acf_type','var' ) || isempty( acf_type )
            acf_type = 'Comb256';
        end
        h_sos = qd_sos.create_array( dist_decorr, acf_type, 'Normal' );
        s = val( h_sos, ca, cb );
    end

    function [ s, h_sos ] = rand( dist_decorr, ca, cb, acf_type )
        %RAND Generates spatially correlated uniformly distributed random numbers
        %
        % Input:
        %   dist_decorr
        %   Vector of decorrelation distances  [1 x M] or [ M x 1 ]
        %
        %   ca
        %   Coordinates for the first mobile device in [m] given as [3 x N] matrix.
        %
        %   cb
        %   Coordinates for the corresponding second mobile device in [m] given as [3 x N] matrix.
        %   This variable must either be empty or have the same size as "ca".
        %
        %   acf_type
        %   String describing the shape of the autocorrelation function and the number of
        %   sinusoids, Default: 'Comb256'
        %
        % Output:
        %   s
        %   Random spatially correlated numbers [ M x N ]
        %
        %   h_sos
        %   A handle to the used qd_sos object

        if ~exist( 'cb','var' )
            cb = [];
        end
        if ~exist( 'acf_type','var' ) || isempty( acf_type )
            acf_type = 'Comb256';
        end
        h_sos = qd_sos.create_array( dist_decorr, acf_type, 'Uniform' );
        s = val( h_sos, ca, cb );
    end

    function [ s, h_sos ] = randi( dist_decorr, ca, imax, cb, acf_type )
        %RANDI Generates spatially correlated random integers between 1 and imax
        %
        % Input:
        %   dist_decorr
        %   Vector of decorrelation distances  [1 x M] or [ M x 1 ]
        %
        %   ca
        %   Coordinates for the first mobile device in [m] given as [3 x N] matrix.
        %
        %   imax
        %   The maximum value [ 1 x 1 ]
        %
        %   cb
        %   Coordinates for the corresponding second mobile device in [m] given as [3 x N] matrix.
        %   This variable must either be empty or have the same size as "ca".
        %
        %   acf_type
        %   String describing the shape of the autocorrelation function and the number of
        %   sinusoids, Default: 'Comb256'
        %
        % Output:
        %   s
        %   Random spatially correlated numbers [ M x N ]
        %
        %   h_sos
        %   A handle to the used qd_sos object

        if ~exist( 'cb','var' )
            cb = [];
        end
        if ~exist( 'acf_type','var' ) || isempty( acf_type )
            acf_type = 'Comb256';
        end
        h_sos = qd_sos.create_array( dist_decorr, acf_type, 'Uniform' );
        s = val( h_sos, ca, cb );
        s = min( max( ceil( s * imax ), 1 ), imax );
    end
end

methods(Static, Access=private)
    function h_sos = create_array( dist_decorr, acf_type, distribution )
        % One generator per decorrelation distance; a decorrelation distance of 0 is replaced by 0.1 m
        h_sos = qd_sos([]);
        for n = 1 : numel( dist_decorr )
            if dist_decorr(n) == 0
                dist_decorr(n) = 0.1;
            end
            h_sos(1,n) = qd_sos( acf_type, distribution, dist_decorr(n) );
        end
    end
end
end

function [ mode, uniform, scale_b ] = eval_mode( h_sos, has_cb )
% Evaluation mode of each object: 1 = one device, 2 = two devices on the same map, 3 = two devices
% on different maps. Two decorrelation distances always use different maps.
n_sos = numel( h_sos );
mode = ones( 1, n_sos );
uniform = false( 1, n_sos );
scale_b = ones( 1, n_sos );
for n = 1 : n_sos
    uniform(n) = strcmpi( h_sos(n).distribution, 'Uniform' );
    d = h_sos(n).dist_decorr;
    if numel( d ) == 2
        scale_b(n) = d(1) / d(2);
    end
    if has_cb
        ph = h_sos(n).sos_phase;
        if scale_b(n) == 1 && all( abs( ph(:,1) - ph(:,2) ) < 1e-6 )
            mode(n) = 2;
        else
            mode(n) = 3;
        end
    end
end
end

function v = sos_eval( h_sos, mode, uniform, scale_b, ca, cb )
% Evaluates objects that share the frequencies and the ACF in one call, output [ N x numel(h_sos) ]
n_sos = numel( h_sos );
phase = zeros( h_sos(1).no_coefficients, n_sos );
for n = 1 : n_sos
    phase(:,n) = h_sos(n).sos_phase(:,1);
end
switch mode
    case 1
        v = quadriga_lib.sos_val( h_sos(1).sos_freq, phase, ca.', [], [], [], [], 1, uniform );
    case 2
        v = quadriga_lib.sos_val( h_sos(1).sos_freq, phase, ca.', cb.', [], h_sos(1).acf, h_sos(1).dist, 1, uniform );
    case 3
        phase_b = zeros( h_sos(1).no_coefficients, n_sos );
        for n = 1 : n_sos
            phase_b(:,n) = h_sos(n).sos_phase(:,2);
        end
        v = quadriga_lib.sos_val( h_sos(1).sos_freq, phase, ca.', cb.', phase_b, [], [], scale_b, uniform );
end
end
