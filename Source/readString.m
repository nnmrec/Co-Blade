function string = readString(fid)
% This function reads a single string from the file with identifier FID and
% then moves the cursor to the beginning of the next line.

string = fscanf(fid, '%s', 1);
% string = lower(string);         % convert to lowercase
fgetl(fid);                     

end % function readString

