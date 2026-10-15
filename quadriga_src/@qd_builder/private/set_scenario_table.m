function file_name = set_scenario_table( h_builder, scenario, check )
%SET_SCENARIO_TABLE Parse conf files and load values into MATLAB
%
%   The conf file is read and validated by 'quadriga_lib.conf_file_read'. A struct of scenario
%   parameters is validated and completed with default values by 'quadriga_lib.scenpar_validate'.
%   A field 'GR_epsilon' with value 0 is accepted and removed; any other value is an error. The
%   ground permittivity is set by 'qd_builder.gr_epsilon_r'.
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

if numel( h_builder ) > 1
    error('QuaDRiGa:qd_builder:ObjectArray','??? "set_scenario_table" is only defined for scalar objects.')
else
    h_builder = h_builder(1,1); % workaround for octave
end

% We can optionally disable the validity check for speedup.
if nargin == 3 && numel( check ) == 1
    check = logical( check );
else
    check = true;
end

if ischar(scenario)
    
    [ sup_scenarios , file_names , file_dirs ] = qd_builder.supported_scenarios;
    
    ind = find( strcmp( scenario, sup_scenarios ), 1 );
    if isempty( ind )
        str = ['??? Scenario "',scenario,'" not found; supported are: '];
        no = numel(sup_scenarios);
        for n = 1:no
            str = [str,sup_scenarios{n}];
            if n<no
                str = [str,', '];
            end
        end
        error('QuaDRiGa:qd_builder:WrongInput',str);
    end
    
    file_name = file_names{ind};
    [ scen, plpar ] = quadriga_lib.conf_file_read( fullfile( file_dirs{ind}, file_name ), check );
    
    h_builder.plpar = plpar;
    h_builder.Pscenpar = scen;
    h_builder.Pscenario = scenario;
    
elseif isstruct( scenario )
    
    if isfield( scenario, 'GR_epsilon' )
        if any( scenario.GR_epsilon(:) ~= 0 )
            error('QuaDRiGa:qd_builder:WrongInput',...
                ['??? "scenpar.GR_epsilon" is not supported. ',...
                'Set the ground permittivity by "qd_builder.gr_epsilon_r".']);
        end
        scenario = rmfield( scenario, 'GR_epsilon' );
    end
    
    h_builder.Pscenpar = quadriga_lib.scenpar_validate( scenario, check );
    h_builder.Pscenario = 'Custom';
    
else
    error('QuaDRiGa:qd_builder:WrongInput','Scenpar has wrong format');
end

end
