function write_conf_file( h_builder, fn, write_defaults, ref, ref_sec, valid_range )
%WRITE_CONF_FILE Writes configuration files from qd_builder objects
%
% Calling object:
%   Single object
%
% Description:
%   Writes the scenario parameters 'qd_builder.scenpar' and the path-loss parameters
%   'qd_builder.plpar' to a configuration file. The file is written by
%   'quadriga_lib.conf_file_write' and can be read back with the same values.
%
% Input:
%   fn
%   String containing the filename
%
%   write_defaults
%   If set to true, all parameters are written to the file. Otherwise, only the parameters that
%   differ from their default value are written. Default: false
%
%   ref
%   String containing the literature reference of the parameters
%
%   ref_sec
%   Reference for each section; Cell-Array containing up to 10 strings. The sections are:
%   (1) Large scale distributions, (2) Large-Scale fading decorrelation distances, (3)
%   Inter-parameter correlations, (4) Model parameters, (5) Decorrelation distance for the
%   small-scale fading spatial consistency, (6) Parameters for the ground reflection extension, (7)
%   Absolute time of arrival offset, (8) Path-loss model, (9) EoD offset, (10) Optional modelling
%   components
%
%   valid_range
%   A 4x2 matrix containing the validity range of the parameters (to be written as text into the
%   conf file). The first row is for the frequency in GHz, the second row contains the distance in
%   m, the third row the allowed BS antenna height and the fourth row contains the allowed MT
%   height in m. The first column is for the lower range, the second for the upper range. Entering
%   NaN into a row prevents the writing of that row to the file.
%
%
% QuaDRiGa Copyright (C) 2011-2021
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

if numel( h_builder ) > 1
    error('QuaDRiGa:qd_builder:ObjectArray','??? "write_conf_file" is only defined for scalar objects.')
else
    h_builder = h_builder(1,1); % workaround for octave
end

if ~exist( 'fn','var' ) || isempty( fn )
    error('QuaDRiGa:qd_builder:write_conf_file',...
        'You did not specify a filename.');
end

if ~exist( 'write_defaults','var' ) || isempty( write_defaults )
    write_defaults = false;
else
    write_defaults = logical( write_defaults(1) );
end

if ~exist( 'ref','var' ) || isempty( ref )
    ref = '';
end

if ~exist( 'ref_sec','var' ) || isempty( ref_sec )
    ref_sec = {};
end

if ~exist( 'valid_range','var' ) || isempty( valid_range )
    valid_range = [];
end

quadriga_lib.conf_file_write( fn, h_builder.scenpar, h_builder.plpar, {}, ref, ref_sec, ...
    valid_range, write_defaults );

end
