function testQF_acdf

data = rand(10000,1);

[ Sh , bins, Sc ] = qf.acdf( data  );

assertEqual( size(bins) , [1,201] )
assertTrue( abs(bins(1)) < 1e-3 )
assertTrue( abs(bins(end)-1) < 1e-3 )

assertEqual( size(Sh) , [201,1] )
assertTrue( all(all(abs( (Sh - ((0:200)'/200)) ) < 0.1)) )
assertEqual( Sh, Sc )

data = rand(2,10000,4);
[ Sh , bins, Sc , mu , sig ] = qf.acdf( data, (0:100)/100, 2,3  );

assertEqual( size(bins) , [1,101] )

assertEqual( size(Sh) , [101,2,4] )
assertEqual( size(Sc) , [101,2] )
assertTrue( all(all(abs( (Sc - ((0:100)'/100)*[1 1]) ) < 0.1)) )
assertTrue( all(all(abs(mu-(1:9)'/10*[1 1]) < 0.05)) )

data = rand(10000,5);
[ Sh , bins, Sc , mu , sig ] = qf.acdf( data  );

assertEqual( size(Sh) , [201,5] )
assertEqual( size(Sc) , [201,1] )
assertTrue( all(all(abs( (Sc - ((0:200)'/200)) ) < 0.1)) )
assertTrue( all(all(abs(mu-(1:9)'/10*1) < 0.05)) )

data = { rand(2999,1), rand(10000,1) };
[ Sh , bins, Sc , mu , sig ] = qf.acdf( data );

assertEqual( size(Sh) , [201,2] )
assertEqual( size(Sc) , [201,1] )
assertTrue( all(all(abs( (Sc - ((0:200)'/200)) ) < 0.1)) )
assertTrue( all(all(abs(mu-(1:9)'/10*1) < 0.05)) )
