function NormAFcoords = readAirfoilCoordFiles(SIM, BLADE)

NormAFcoords(BLADE.NUM_SEC,1).x = [];    % preallocate structure array
NormAFcoords(BLADE.NUM_SEC,1).y = [];    % preallocate structure array
for i = 1:BLADE.NUM_SEC
    
    % open the airfoil file
    fid = fopen([SIM.airfoilDir filesep BLADE.afFile{i}], 'r');
    if fid == -1
        error(['ERROR: Could not locate and open file ' [SIM.airfoilDir filesep BLADE.afFile{i}]]);
    end
    
    % skip the header lines
    for j = 1:4
        fgetl(fid);
    end
    
    % read the table of x-y coordinates
    xy = fscanf(fid, '%f %f', [2, Inf]);
    
    % close the airfoil file
    fclose(fid);
    
    % transpose, so the xy data matches the orientation of the file
    xy = xy';
    
    NormAFcoords(i).x = xy(:,1);
    NormAFcoords(i).y = xy(:,2);
    
    % check for errors in the user inputs
    if length(xy) < 4
        error(['Error: In file ' BLADE.afFile{i} ' the coordinates should contain at least 4 points.']);
    end
    if xy(1,1) ~= 0 || xy(1,2) ~= 0
        error(['Error: In file ' BLADE.afFile{i} ' the leading edge coordinate is not located at (x,y) = (0,0).']);
    end
    [x_max ii] = max(xy(:,1)); 
    if x_max > 1
        error(['Error: In file ' BLADE.afFile{i} ' the maximum x-coordinate exceeds the chord boundary (it should be equal to 1).']);
    end
    x_upper = xy(1:ii,1);
    y_upper = xy(1:ii,2);
    x_lower = xy(ii:end,1);
    y_lower = xy(ii:end,2); 
    if y_upper(2) < y_lower(end)
        % note: this is not a perfect test for clockwise ordering, but I think it
        % will be sufficient given the other requirements for coordinate ordering
        error(['Error: In file ' BLADE.afFile{i} ' the coordinates are not labeled in clock-wise order.']);
    end
    if any( diff(x_upper) <= 0 )
        error(['Error: In file ' BLADE.afFile{i} ' the upper surface is not a single valued function.']);
    end
    if any( diff(x_lower) >= 0 )
        % note: we excluded the trailing edge in this check, to allow for finite thickness trailing edges
        error(['Error: In file ' BLADE.afFile{i} ' the lower surface is not a single valued function.']);
    end
      
end

end % function readAirfoilCoordFiles
