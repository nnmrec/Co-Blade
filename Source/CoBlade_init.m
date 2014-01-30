function [SIM, ANLS, OPT, ENV, BLADE, WEB, OUT, MATS, AF, Coord] = CoBlade_init(inputFile)

%% set the Co-Blade version and directories
SIM.version = '2.0.0_rc1';

[SIM.pathStr SIM.case SIM.inpExt]      = fileparts(inputFile);
SIM.inputFile                          = [pwd filesep 'Inputs' filesep SIM.case SIM.inpExt];
[SIM, ANLS, OPT, ENV, BLADE, WEB, OUT] = readInput_CoBlade(SIM);

SIM.rootDir     = pwd;
SIM.sourceDir   = [pwd filesep 'Source'];
SIM.airfoilDir  = [SIM.rootDir filesep 'Inputs'  filesep 'Airfoil_Data'];
SIM.materialDir = [SIM.rootDir filesep 'Inputs'  filesep 'Material_Data'];
SIM.laminateDir = [SIM.rootDir filesep 'Inputs'  filesep 'Laminate_Data'];
SIM.optimDir    = [SIM.rootDir filesep 'Inputs'  filesep 'Optimization_Data'];
SIM.outputDir   = [SIM.rootDir filesep 'Outputs' filesep SIM.case];
SIM.logfile     = [SIM.outputDir filesep SIM.case '_Log.txt'];
SIM.matfile     = [SIM.outputDir filesep SIM.case '.mat'];

%% check for errors on user inputs within main input file
checkInpErrors(ANLS,OPT,ENV,BLADE,WEB)

%% Define the output directories
if exist([SIM.rootDir filesep 'Outputs'],'dir') == 0
    mkdir([SIM.rootDir filesep 'Outputs']);
end
if exist(SIM.outputDir,'dir') == 7
    rmdir(SIM.outputDir,'s');
    mkdir(SIM.outputDir);
else
    mkdir(SIM.outputDir);
end

% Start a log of all text and command line output. Read this file if
% something goes wrong and follow the error handling messages to debug.
diary(SIM.logfile);
fprintf(1, '\n\n Executing case %s with Co-Blade v%s.\n', SIM.case, SIM.version);

%% reset ALL the random number generators (rand, randn, randi)
rng(SIM.RAND_SEED);

% echo user parameters to output directory
copyfile([SIM.inputFile],[SIM.outputDir filesep SIM.case '_echo' SIM.inpExt]);

% Deal with the BModes executable
SIM.BmodesExe = [SIM.rootDir filesep 'Source' filesep 'BModes.exe'];
if exist(SIM.BmodesExe) == 0
    error(['ERROR: could not find BModes executable ' SIM.BmodesExe])
end
if ANLS.N_MODES > 0 % then we will use BModes
    copyfile(SIM.BmodesExe,SIM.outputDir);
end

%% Read the materials input file
MATS = readMaterialsFile(SIM,BLADE,OPT);

%% Read the normalized airfoil coordinates (*.prof files) and interpolate (if requested by user) 
AF.OrigNormAFcoords = readAirfoilCoordFiles(SIM, BLADE);
AF.NormAFcoords     = interpAirfoilCoords(AF.OrigNormAFcoords, BLADE);	% structure array of the interpolated chord normalized airfoil x-y coordinates, in the reference airfoil coordinate system

%% Define the coordinate systems and the tranformations matrices between them
Coord = defineCoordSystems(BLADE);

%% Determine the number of webs and cells at each station
existWeb = zeros(BLADE.NUM_SEC, WEB.NUM_WEBS);
if OPT.OPTIMIZE
    existWeb(OPT.INB_STN:OPT.OUB_STN,:) = 1;    
else
    for n = 1:WEB.NUM_WEBS
        existWeb(WEB.inbStn(n):WEB.oubStn(n),n) = 1;
    end
end
WEB.nWebs    = sum(existWeb,2);
BLADE.nCells = WEB.nWebs + 1;
BLADE.zSec   = BLADE.zFrac .* BLADE.BLD_LENGTH; % z-coordinate in blade coordinate system (m)

%% if pitch axis was not defined, define it now via optimization
if OPT.OPT_PITAXIS
    BLADE.pitAxis = definePitchAxis(BLADE, OPT);
end

end % function CoBlade_init

