function fitVal = structFitness(xo, SIM, ANLS, OPT, ENV, BLADE, WEB, AF, MATS, Coord, z_oub, z_CP, nSegsTopBot)

    
%% make sure the xo is a column vector
xo = xo(:);

%% create the laminate data needed for the structural analysis
[WEB SECNODES LamData] = defineLaminateData(xo, SIM, OPT, BLADE, WEB, MATS, z_oub, z_CP, nSegsTopBot);
 
%% call the structural analysis
try
    [Panel StrProps AppLoads ResLoads Disp NormS ShearS Buckle MidPlane LaminaSS Modes] ...
     = structAnalysis(SIM, ANLS, ENV, BLADE, WEB, MATS, LamData, AF, SECNODES, Coord);
catch
    % something unspeakable happened, oh the horror
    fprintf(1, 'WARNING: a severe error occured within structAnalysis...continuing anyways.\n');
    fitVal = Inf;   % an error occured, return a placeholder value and hope the optimization can recover from it
    return
end

%% compute the fitness value
fitVal = fitnessFunction(SIM, ANLS, OPT, BLADE, WEB, StrProps, LaminaSS, Buckle, Disp, Modes);

%% write the laminate data as text laminate input files (if requested for debugging purposes)
if OPT.WRITE_STR                  
    writeInpFileLaminate(SIM, BLADE, WEB, SECNODES, LamData);
end

end % function structFitness

