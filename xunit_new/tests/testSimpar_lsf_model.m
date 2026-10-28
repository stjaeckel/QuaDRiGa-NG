function testSimpar_lsf_model

s = qd_simulation_parameters;
assertEqual( s.lsf_model, '38.901' );
c0 = s.checksum;

s.lsf_model = 'legacy';
assertEqual( s.lsf_model, 'legacy' );
assertFalse( s.checksum == c0 );

c = copy( s );
assertEqual( c.lsf_model, 'legacy' );

s.lsf_model = '38.901';
assertEqual( s.checksum, c0 );
assertEqual( c.lsf_model, 'legacy' );

assertExceptionThrown( @() set_model( s, 'Legacy' ), 'QuaDRiGa:qd_simulation_parameters:lsf_model' );
assertExceptionThrown( @() set_model( s, 1 ), 'QuaDRiGa:qd_simulation_parameters:lsf_model' );
assertEqual( s.lsf_model, '38.901' );

end

function set_model( s, value )
s.lsf_model = value;
end
