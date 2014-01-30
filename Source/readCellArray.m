function cellArray = readCellArray(fid, aryFormat, aryRows)
% This function reads a cell array from the file with identifier FID and 
% then moves the cursor to the beginning of the next line.  The size of 
% CELLARRAY is 1 row by "the number of conversion specifiers defined in 
% ARYFORMAT" columns.  Each element of CELLARRAY is an array with size
% ARYROWS rows by 1 column.  

cellArray = textscan(fid, aryFormat, aryRows);
fgetl(fid); 

end % function readCellArray

