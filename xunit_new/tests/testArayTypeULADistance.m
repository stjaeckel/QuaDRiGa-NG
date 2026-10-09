function testArayTypeULADistance
% Check the distances between ULA-Elements

types = {'ula2','ula4','ula8'};
for n=1:numel(types)
    a = qd_arrayant(types{n});
    for n = 1:a.no_elements
        l=[];
        for m=1:3
            l(m,:) = a.element_position(m,setdiff(1:a.no_elements,n)) - a.element_position(m, n ) ;
        end
        assertTrue( all( sqrt(sum(l.^2,1))>0.09 ) );
    end
end
