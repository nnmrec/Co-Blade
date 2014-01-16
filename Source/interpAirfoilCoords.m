function NormAFcoords = interpAirfoilCoords(OrigNormAFcoords, BLADE)

switch BLADE.INTERP_AF 
    case 'none'
        % we still need to interpolate the coordinates to measure the percent
        % thickness, but we will use the original airfoil coordinates in the end
        x_cUp = linspace(0, 1, 200 + 1)';
        
    case 'equal'
        x_cUp = linspace(0, 1, BLADE.N_AF/2 + 1)';
        
    case 'cosine'
        x_cUp = cosspace(0, 1, BLADE.N_AF/2 + 1, 'both');
end
x_cLo = flipud( x_cUp );

NormAFcoords(BLADE.NUM_SEC,1).x           = [];       
NormAFcoords(BLADE.NUM_SEC,1).y           = [];      
NormAFcoords(BLADE.NUM_SEC,1).maxPerThick = [];       
for i = 1:BLADE.NUM_SEC
        
    oldX           = OrigNormAFcoords(i).x;
    oldY           = OrigNormAFcoords(i).y;
    [unused i_TE]  = min(abs(oldX - 1)); % find the index of the trailing edge 
    oldXup         = oldX(1:i_TE);
    oldYup         = oldY(1:i_TE);
    oldXlo         = [oldX(i_TE:end); 0];
    oldYlo         = [oldY(i_TE:end); 0];

    newYup = interp1(oldXup, oldYup, x_cUp);
    newYlo = interp1(oldXlo, oldYlo, x_cLo);
    
    NormAFcoords(i).x           = [x_cUp; x_cLo(2:end-1)];
    NormAFcoords(i).y           = [newYup; newYlo(2:end-1)];  
    NormAFcoords(i).maxPerThick = max(newYup - flipud(newYlo));
end

switch BLADE.INTERP_AF
    case 'none'
        % overwrite the interpolated coordinates with the original coordinates
        for i = 1:BLADE.NUM_SEC
            NormAFcoords(i).x = OrigNormAFcoords(i).x;
            NormAFcoords(i).y = OrigNormAFcoords(i).y;
        end
    otherwise
        return
end

end % function interpAirfoilCoords

