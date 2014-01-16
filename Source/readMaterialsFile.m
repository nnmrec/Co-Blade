function MATS = readMaterialsFile(SIM, BLADE, OPT)

fid = fopen([SIM.materialDir filesep BLADE.MATS_FILE], 'r');
if fid == -1
    error(['ERROR: Could not locate and open file ' [SIM.materialDir filesep BLADE.MATS_FILE]]);
end
    
% skip the header lines
for j = 1:3
    fgetl(fid);
end

% read the table of materials data
matData = textscan(fid, '%f %f %f %f %f %f %f %f %f %f %f %s');

% close the materials file
fclose(fid);

%% check for errors on the materials data
if size(matData, 2) ~= 12
    error(['ERROR: Expected to read 12 columns in material data file, but the number of columns is ' num2str(size(matData, 2))]);
end

MATS.matID   = matData{:,1};  
MATS.E11     = matData{:,2};
MATS.E22     = matData{:,3};
MATS.G12     = matData{:,4};
MATS.nu12    = matData{:,5};
MATS.density = matData{:,6};
MATS.s11_fT  = matData{:,7};
MATS.s11_fC  = matData{:,8};
MATS.s22_yT  = matData{:,9};
MATS.s22_yC  = matData{:,10};
MATS.s12_y   = matData{:,11};
MATS.matName = matData{:,12};
MATS.numMats = size(MATS.matID, 1);
if OPT.OPTIMIZE && MATS.numMats ~= 8
    error(['ERROR: When OPTIMIZE = true, the number of rows for material data is expected to be equal to 8, but instead the number of rows is ' num2str(size(matData{:,1}, 1))]);
end
if numel( unique(MATS.matID) ) ~= numel(MATS.matID) || any(~mod(MATS.matID,1) == 0)
    error('ERROR: The values for matID must be unique positive integers.')
end
if any(MATS.E11 <=0)
    error('ERROR: Values for E11 must be positive.');
end
if any(MATS.E22 <=0)
    error('ERROR: Values for E22 must be positive.');
end
if any(MATS.G12 <=0)
    error('ERROR: Values for G12 must be positive.');
end
if any(MATS.density < 0)
    error('ERROR: Values for density must be positive.');
end
if any(MATS.s11_fT < 0)
    error('ERROR: Values for s11_fT are expected to be positive.');
end
if any(MATS.s11_fC > 0)
    error('ERROR: Values for s11_fC are expected to be negative.');
end
if any(MATS.s22_yT < 0)
    error('ERROR: Values for s22_yT are expected to be positive.');
end
if any(MATS.s22_yC > 0)
    error('ERROR: Values for s22_yC are expected to be negative.');
end
if any(MATS.s12_y < 0)
    error('ERROR: Values for s12_y are expected to be positive.');
end
if any( MATS.nu12 > sqrt(MATS.E11./MATS.E22) )
    error('ERROR: Material properies are not physically possible. The relation nu12 < sqrt(E11/E22) must be satisfied.');
end

end % function readMaterialsFile

