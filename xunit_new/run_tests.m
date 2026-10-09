function run_tests
tests = dir('tests/test*.m');

% Set paths:
current_path = pwd;
if isempty(strfind(current_path,'xunit_new')) %#ok
    error('Must be in xunit_new folder');
end

tmp = which('qd_simulation_parameters');
if isempty(strfind(tmp,'quadriga_src')) %#ok
    addpath([current_path(1:end-9),'quadriga_src']);
end

tmp = which('quadriga_lib.version');
if isempty(strfind(tmp,['+quadriga_lib',filesep,'version.mex'])) %#ok
    addpath([current_path(1:end-9),'quadriga_lib']);
end

tmp = which('MOxUnitTestSuite');
if isempty(tmp)
    current_dir = pwd;
    cd('../external/MOxUnit-master/MOxUnit');
    moxunit_set_path();
    cd(current_dir);
end

% Get version number
r = which('qd_simulation_parameters.m');
ri = regexp(r,'/quadriga_src/@', 'once');

clc
disp(['Quadriga v',qd_simulation_parameters.version])
disp(['quadriga-lib v',quadriga_lib.version])
disp(r(1:ri-1))

if true % ~isempty(  strfind( r(1:ri-1), '_v2.6.') ) %#ok
    exclude_list = {'testChannel_fr2cir','bls'};
else
    exclude_list = {};
end

warning('off','all')
addpath( 'tests' );

N = numel( tests );

test_suite=MOxUnitTestSuite();
Ne = 0;
for n = 1 : N
    subFunctionName = tests(n).name(1:end-2);
    if any( any(strcmp( subFunctionName, exclude_list )) )
        disp(['Exclude test: ',subFunctionName]);
        Ne = Ne + 1;
    else
        test_case = MOxUnitFunctionHandleTestCase(subFunctionName,...
            'run_tests', str2func( subFunctionName ));
        test_suite=addTest(test_suite, test_case);
    end
end

disp(['Running ',num2str(N-Ne),' tests'])

tic
disp(run(test_suite));
toc

disp(' ');
disp('Testing config files:');
test_all_config_files;
toc

rmpath( 'tests' );
warning('on','all')
