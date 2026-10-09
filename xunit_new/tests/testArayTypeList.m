function testArayTypeList
% Generate all array antennas on the list

types =  {'omni', 'dipole', 'half-wave-dipole', ...
    'custom', '3gpp', '3gpp-macro', '3gpp-3d', '3gpp-mmw', ...
    'parametric', 'rhcp-dipole', 'lhcp-dipole', 'lhcp-rhcp-dipole', ...
    'xpol', 'ula2', 'ula4', 'ula8', 'patch', 'multi', 'vehicular'};

for n = 1:numel(types)
    try
        a = qd_arrayant(types{n});
        assertEqual(  a.name  ,  types{n} );
        assertEqual(  numel(a.azimuth_grid)  , 361 );
    catch
        disp(types{n})
        assertTrue( false );
    end
end
