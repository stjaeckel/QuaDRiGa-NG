function testArayTypeOmni
% Generate omni

a = qd_arrayant('omni');
assertEqual(  a.Fa  ,  ones(181,361) );
assertEqual(  a.Fb  ,   zeros(181,361) );