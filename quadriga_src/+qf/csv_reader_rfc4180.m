function rows = csv_reader_rfc4180(filename, delimiter, quotechar)
% RFC 4180 CSV reader: handles quotes, commas, and newlines in quoted fields.
if nargin < 2 || isempty(delimiter),  delimiter  = ',';  end
if nargin < 3 || isempty(quotechar),  quotechar  = '"';  end

fid = fopen(filename,'r','n','UTF-8');
assert(fid > 0, 'Cannot open file: %s', filename);
data = fread(fid,'*char')';
fclose(fid);

rows = {};
row  = {};
field = '';
inQ = false;
i = 1; n = numel(data);

while i <= n
    ch = data(i);

    if inQ
        if ch == quotechar
            if i < n && data(i+1) == quotechar
                field = [field quotechar];  % escaped quote ("")
                i = i + 2;  continue
            else
                inQ = false;  i = i + 1;  continue
            end
        else
            field = [field ch];  i = i + 1;  continue
        end
    else
        if ch == quotechar
            inQ = true;  i = i + 1;  continue
        elseif ch == delimiter
            row{end+1} = field;  field = '';  i = i + 1;  continue
        elseif ch == 10 || ch == 13   % LF or CR
            row{end+1} = field;  field = '';
            if ch == 13 && i < n && data(i+1) == 10  % CRLF
                i = i + 2;
            else
                i = i + 1;
            end
            rows{end+1} = row;  row = {};  continue
        else
            field = [field ch];  i = i + 1;  continue
        end
    end
end
end
