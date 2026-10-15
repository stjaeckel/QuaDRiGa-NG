function [ scenarios, file_names, file_dir ] = supported_scenarios( parse_shortnames )
%SUPPORTED_SCENARIOS Returns a list of supported scenarios
%
% Calling object:
%   None (static method)
%
% Description:
%   The configuration files are searched in the current working directory first and then in the
%   folder 'config' of QuaDRiGa. A file in the current directory replaces a file of the same name
%   in the config folder. The list is created by 'quadriga_lib.conf_file_list'.
%
% Input:
%   parse_shortnames
%   This optional parameter can disable (0) the shortname parsing. This is faster.
%   By default shortname parsing is enabled (1)
%
% Output:
%   scenarios
%   A cell array containing the scenario names. Those can be used in 'qd_track.scenario'
%
%   file_names
%   The names of the configuration files for each scenario
%
%   file_dir
%   The directory where each file was found. You can place configuration file also in you current
%   working directory
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

% Parse input variables.
if nargin > 0
    parse_shortnames = logical(parse_shortnames(1));
else
    parse_shortnames = true;
end

% The config folder is next to the class folder "@qd_builder"
config_folder = fullfile( fileparts( fileparts( mfilename('fullpath') ) ), 'config' );

[ scenarios, file_names, file_dir ] = quadriga_lib.conf_file_list( { '', config_folder }, parse_shortnames );

scenarios  = scenarios.';
file_names = file_names.';
file_dir   = file_dir.';

end
