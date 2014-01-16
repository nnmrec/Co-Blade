function scalar = readVector(fid,N)
% This function reads a vector (1D array of length N) from the file with identifier FID and
% then moves the cursor to the beginning of the next line.

scalar = fscanf(fid, '%f', [1 N]);
fgetl(fid);                     

end % function readVector

