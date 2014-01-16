function scalar = readScalar(fid)
% This function reads a single scalar from the file with identifier FID and
% then moves the cursor to the beginning of the next line.

scalar = fscanf(fid, '%f', 1);
fgetl(fid);                     

end % function readScalar

