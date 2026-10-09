function testSOS_val_NumericStability

for n = 1:100
    
    pR1 = rand( 3,100 )*100;
    pR2 = rand( 3,20  )*100;
    pR2(:,[1,5]) = pR1(:,[1,50]);
    
    pT1 = rand( 3,100 )*100;
    pT2 = rand( 3,20  )*100;
    pT2(:,[1,5]) = pT1(:,[1,50]);
    
    pT1(:,2) = pR1(:,3);
    pT1(:,3) = pR1(:,2);
    
    sos = qd_sos;

    s1 = val( sos, pR1 );
    s2 = val( sos, pR2 );
    
    assertTrue( abs( s1(1) - s2(1) ) < 1e-14 )
    assertTrue( abs( s1(50) - s2(5) ) < 1e-14 )
    
    sos = qd_sos('Gauss150', 'Normal' , 2.3 );
    
    s1 = val( sos, pR1, pT1 );
    s2 = val( sos, pR2, pT2 );
    
    assertTrue( abs( s1(1) - s2(1) ) < 1e-14 )
    assertTrue( abs( s1(50) - s2(5) ) < 1e-14 )
    assertTrue( abs( s1(2) - s1(3) ) < 1e-14 )
       
 
end







