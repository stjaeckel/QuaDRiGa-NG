function testChannel_HDF5

% Delete file
fn = 'test.hdf5';
if exist( fn,'file' )
    delete(fn);
end

c = qd_channel;

c.name = 'Buy_Bitcoin';
c.center_frequency = 21e6;
c.coeff = rand( 2,3,4,5 ) + 1j* rand(2,3,4,5);
c.delay = rand( 2,3,4,5 );
c.initial_position = 3;
c.tx_position = [0;0;1];
c.rx_position = rand(3,5);

c(1,2) = c(1,1).copy;
c(1,2).name = 'Not_ethereum';
c(1,2).center_frequency = 100e6;
c(1,2).coeff = rand( 2,3,4,6 ) + 1j* rand(2,3,4,6);
c(1,2).delay = rand( 2,3,4,6 );
c(1,2).initial_position = 4;
c(1,2).tx_position = [0;0;2];
c(1,2).rx_position = rand(3,6);

hdf5_write(c,fn);

% Load layout
[~,layout,has_data] = qd_channel.hdf5_read(fn,0);
assertEqual( layout, uint32([1024,64,1,1]) );
assertEqual( size(has_data), [1024,64] );

% Load channel
d = qd_channel.hdf5_read(fn);

% Compare
assertEqual( size(c), size(d) );
for n = 1:2
    assertEqual( c(1,n).name, d(1,n).name );
    assertEqual( double(single(c(1,n).coeff)), d(1,n).coeff );
    assertEqual( double(single(c(1,n).delay)), d(1,n).delay );
    assertEqual( double(single(c(1,n).tx_position)), d(1,n).tx_position );
    assertEqual( double(single(c(1,n).rx_position)), d(1,n).rx_position );
    assertEqual( int32(c(1,n).initial_position), d(1,n).initial_position );
end

% Load single channel
d = qd_channel.hdf5_read(fn,1,2);
assertEqual( 'Not_ethereum', d(1,1).name );

% Add another row of channels
c(1,3) = c(1,2).copy;
c(1,3).name = 'Bitcoin_Cash_Is_Trash';
hdf5_write(c,fn,[],2,2:4);

% Check layout
[~,~,has_data] = qd_channel.hdf5_read(fn,0);
assertEqual( sum(has_data(:)), 5 );

% Read all
d = qd_channel.hdf5_read(fn);
assertEqual( [2,4], size(d) );
assertEqual( 'Buy_Bitcoin', d(1,1).name );
assertEqual( 'Not_ethereum', d(1,2).name );
assertEqual( 'empty', d(1,3).name );
assertEqual( 'empty', d(1,4).name );
assertEqual( 'empty', d(2,1).name );
assertEqual( 'Buy_Bitcoin', d(2,2).name );
assertEqual( 'Not_ethereum', d(2,3).name );
assertEqual( 'Bitcoin_Cash_Is_Trash', d(2,4).name );

% Read second row with 2 snapshots in reverse order
d = qd_channel.hdf5_read(fn,2,[],1,1,[2,1]);
assertEqual( [1,3], size(d) );
assertEqual( 'Buy_Bitcoin', d(1,1).name );
assertEqual( 'Not_ethereum', d(1,2).name );
assertEqual( 'Bitcoin_Cash_Is_Trash', d(1,3).name );

assertTrue( d(1,1).no_snap == 2 );
assertTrue( d(1,2).no_snap == 2 );
assertTrue( d(1,3).no_snap == 2 );

assertEqual( double(single(c(1,1).coeff(:,:,:,[2,1]))), d(1,1).coeff );
assertEqual( double(single(c(1,2).coeff(:,:,:,[2,1]))), d(1,2).coeff );

assertEqual( double(single(c(1,1).delay(:,:,:,[2,1]))), d(1,1).delay );
assertEqual( double(single(c(1,2).delay(:,:,:,[2,1]))), d(1,2).delay );

assertEqual( double(single(c(1,1).rx_position(:,[2,1]))), d(1,1).rx_position );
assertEqual( double(single(c(1,2).rx_position(:,[2,1]))), d(1,2).rx_position );

% Add a par-struct
par = struct;
par.string = 'Buy Bitcoin!';
par.double = 21e6;
par.single = single(pi);
par.uint32 = uint32(21);
par.int32 = int32(-11001001);
par.uint64 = uint64( 21e6*100e6 );
par.int64 = -int64( 21e6*100e6 );
par.double_Col = [0:0.1:1]';
par.single_Col = -single([0:0.1:1]');
par.uint32_Col = uint32([14:18]');
par.int32_Col = -int32([14:18]');
par.uint64_Col = uint64( 21e6*100e6 + [0,1]' );
par.int64_Col = -int64( 21e6*100e6 + [0,1]' );
par.double_Row = [1:0.1:2];
par.single_Row = -single([1:0.1:2]);
par.uint32_Row = uint32([17:19]);
par.int32_Row = -int32([12:19]);
par.uint64_Row = uint64( 21e6*100e6 + [2,3] );
par.int64_Row = -int64( 21e6*100e6 + [3,4] );
par.double_Mat = rand(4);
par.single_Mat = -single(rand(5));
par.uint32_Mat = randi(10,3,'uint32');
par.int32_Mat = -randi(10,4,'int32');
par.uint64_Mat = uint64(randi(10,5,'uint32')) + uint64( 21e6*100e6) ;
par.int64_Mat = int64(randi(10,6,'int32'))  - int64( 21e6*100e6) ;
par.double_Cube = rand(4,3,2);
par.single_Cube = -single(rand(5,4,3));
par.uint32_Cube = randi(10,3,3,4,'uint32');
par.int32_Cube = -randi(10,4,5,6,'int32');
par.uint64_Cube = uint64(randi(10,5,6,7,'uint32')) + uint64( 21e6*100e6) ;
par.int64_Cube = int64(randi(10,6,7,8,'int32'))  - int64( 21e6*100e6) ;
c(1,1).par = par;

% Write to a new row
hdf5_write(c(1,1),fn,[],3,1);
d = qd_channel.hdf5_read(fn,3,1);

fieldsPar = fieldnames(c(1,1).par);
fieldsParR = fieldnames(d(1,1).par);
for n = 1:length(fieldsPar)
    field = fieldsPar{n};
    assertEqual( class(c(1,1).par.(field)), class(d(1,1).par.(field)));  % Same data type
    assertTrue(  isequal(c(1,1).par.(field), d(1,1).par.(field)) );      % Same data
end

% Reading from empty data
d = qd_channel.hdf5_read(fn,[3,4],[2,3]);
assertEqual( [2,2], size(d) );
assertEqual( 'empty', d(1,1).name );
assertEqual( 'empty', d(1,2).name );
assertEqual( 'empty', d(2,1).name );
assertEqual( 'empty', d(2,2).name );

% Read first column
d = qd_channel.hdf5_read(fn,[],1);
assertEqual( [2,1], size(d) );
assertEqual( 'Buy_Bitcoin', d(1,1).name );
assertEqual( 'Buy_Bitcoin', d(2,1).name );
assertTrue( isempty(d(1,1).par) );
assertFalse( isempty(d(2,1).par) );

% Check out-of-bound error
try
    d = qd_channel.hdf5_read(fn,1,1,2);
    error('moxunit:exceptionNotRaised', 'Expected an error!');
catch ME
    if (strcmp(ME.identifier, 'moxunit:exceptionNotRaised'))
        error('moxunit:exceptionNotRaised', 'Expected an error!');
    end
end

% Check snapshot out-of-bound error
try
    d = qd_channel.hdf5_read(fn,1,[],1,1,6);
    error('moxunit:exceptionNotRaised', 'Expected an error!');
catch ME
    if (strcmp(ME.identifier, 'moxunit:exceptionNotRaised'))
        error('moxunit:exceptionNotRaised', 'Expected an error!');
    end
end

% Trying to alter the layout should fail
try
    hdf5_save(c,fn,size(c));
    error('moxunit:exceptionNotRaised', 'Expected an error!');
catch ME
    if (strcmp(ME.identifier, 'moxunit:exceptionNotRaised'))
        error('moxunit:exceptionNotRaised', 'Expected an error!');
    end
end

% Check optional parameters
path_gain = rand(5,5);
path_length = rand(5,5);
path_polarization = rand(8,5,5);
path_angles = rand(5,4,5);
fbs_pos = rand(3,5,5);
lbs_pos = rand(3,5,5);
no_interact = [1 2 3 4 5 ; 5 4 3 2 1 ; 1 1 1 1 1 ; 1 2 3 4 5 ; 1 2 3 4 5]';
interact_coord = rand(3,15,5);
interact_coord(:,6:end,3) = 0;
rx_orientation = rand(3,5);
tx_orientation = rand(3,5);
tx_position = rand(3,5);

chan = struct( 'name', 'xxx', 'rx_position', c(1,1).rx_position, 'tx_position', tx_position, ...
    'path_gain', path_gain, 'path_length', path_length, 'path_polarization', path_polarization, ...
    'path_angles', path_angles, 'fbs_pos', fbs_pos, 'lbs_pos', lbs_pos, 'no_interact', uint32(no_interact), ...
    'interact_coord', interact_coord, 'rx_orientation', rx_orientation, 'tx_orientation', tx_orientation );
quadriga_lib.hdf5_write_channel(fn, chan, [], 5, 1, 1, 1 );

d = qd_channel.hdf5_read(fn,5);
assertEqual( [1,1], size(d) );

assertEqual( 'xxx', d.name );
assertEqual( double(single(c(1,1).rx_position)), d(1,1).rx_position );
assertEqual( double(single(tx_position)), d(1,1).tx_position );
assertEqual( double(single(path_gain)), d(1,1).par.path_gain );
assertEqual( double(single(path_length)), d(1,1).par.path_length );
assertEqual( double(single(path_polarization)), d(1,1).par.path_polarization );
assertEqual( double(single(path_angles)), d(1,1).par.path_angles );
assertEqual( double(single(fbs_pos)), d(1,1).par.fbs_pos );
assertEqual( double(single(lbs_pos)), d(1,1).par.lbs_pos );
assertEqual( uint32(no_interact), d(1,1).par.no_interact );
assertEqual( double(single(interact_coord)), d(1,1).par.interact_coord );
assertEqual( double(single(rx_orientation)), d(1,1).rx_orientation );
assertEqual( double(single(tx_orientation)), d(1,1).tx_orientation );

% Delete file
delete(fn);

% Set custom storage layout
hdf5_write(c(1,1),fn,[1,2,3,4],1,2,3,4);
[d,layout,has_data] = qd_channel.hdf5_read(fn);
assertEqual( 'Buy_Bitcoin', d(1,1).name );
assertEqual( layout, uint32([1,2,3,4]) );
assertTrue( has_data(end)==uint32(1) );
delete(fn);

end