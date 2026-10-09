function testChannel_mat_save_load

for n = 1:4
    coeff = rand(2,3,4,n) + 1j*rand(2,3,4,n);
    if n == 3
        delay = rand(2,3,4,n);
    else
        delay = rand(4,n);
    end
    c(1,n) = qd_channel(coeff,delay);
    c(1,n).name = ['bs1_mt',num2str(n)];
    c(1,n).center_frequency = rand*1e9;
    if n < 2.5
        c(1,n).tx_position = rand(3,1);
        c(1,n).rx_position = rand(3,n);
    elseif n == 3
        c(1,n).tx_position = rand(3,n);
        c(1,n).rx_position = rand(3,n);
    end
    if n > 1.5
        par.test1 = rand;
        par.test2 = rand(2,2,3,3)+1j*rand(2,2,3,3);
        par.test3 = 'hello_world';
        c(1,n).par = par;
    end
end

mat_save(c,'mat_save_test.mat');

d = qd_channel.mat_load('mat_save_test.mat');

% Compare results
assertEqual( cat(2,c.name), cat(2,d.name) );
assertEqual( cat(2,c.version), cat(2,d.version) );
assertEqual( single(cat(2,c.center_frequency)), cat(2,d.center_frequency) );
assertEqual( single(cat(2,c.tx_position)), cat(2,d.tx_position) );
assertEqual( single(cat(2,c.rx_position)), cat(2,d.rx_position) );
assertEqual( cat(2,c.individual_delays), cat(2,d.individual_delays) );

for n = 1:numel(c)
    assertEqual( single(c(1,n).coeff), d(1,n).coeff );
    assertEqual( single(c(1,n).delay), d(1,n).delay );
    assertEqual( isempty(c(1,n).par), isempty(d(1,n).par) );
    if ~isempty(c(1,n).par)
        assertEqual( single(c(1,n).par.test1), d(1,n).par.test1 );
        assertEqual( single(c(1,n).par.test2), d(1,n).par.test2 );
        assertEqual( c(1,n).par.test3, d(1,n).par.test3 );
    end
end

[d,dims] = qd_channel.mat_load('mat_save_test.mat',0);
assertEqual( dims, [1,4,1,1] );

d = qd_channel.mat_load('mat_save_test.mat',[],[2,2]);

% Compare results
assertEqual( c(1,2).name, d(1,1).name );
assertEqual( c(1,2).name, d(1,2).name );
assertEqual( single(c(1,2).tx_position), d(1,1).tx_position );
assertEqual( single(c(1,2).tx_position), d(1,2).tx_position );
assertEqual( single(c(1,2).rx_position), d(1,1).rx_position );
assertEqual( single(c(1,2).rx_position), d(1,2).rx_position );
assertEqual( single(c(1,2).coeff), d(1,1).coeff );
assertEqual( single(c(1,2).coeff), d(1,2).coeff );
assertEqual( single(c(1,2).delay), d(1,1).delay );
assertEqual( single(c(1,2).delay), d(1,2).delay );
assertEqual( single(c(1,2).par.test2), d(1,1).par.test2 );
assertEqual( single(c(1,2).par.test2), d(1,2).par.test2 );

d = qd_channel.mat_load('mat_save_test.mat',[],[],[],[],'par');

% Compare results
for n = 1:numel(c)
    assertTrue( isempty( d(1,n).coeff ) );
    assertTrue( isempty( d(1,n).delay ) );
    assertEqual( isempty(c(1,n).par), isempty(d(1,n).par) );
    if ~isempty(c(1,n).par)
        assertEqual( single(c(1,n).par.test1), d(1,n).par.test1 );
        assertEqual( single(c(1,n).par.test2), d(1,n).par.test2 );
        assertEqual( c(1,n).par.test3, d(1,n).par.test3 );
    end
end

delete('mat_save_test.mat')

