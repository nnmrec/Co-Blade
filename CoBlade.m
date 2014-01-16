function CB = CoBlade(inputFile)
%Co-Blade - Analysis & Design of Composite Rotors.
%   _____               ______ _           _      
%  /  __ \              | ___ \ |         | |     
%  | /  \/ ___  ______  | |_/ / | __ _  __| | ___ 
%  | |    / _ \ |_____| | |__ \ |/ _` |/ _` |/ _ \
%  | \__/\ (_) |        | |_/ / | (_| | (_| |  __/
%   \____/\___/         \____/|_|\__,_|\__,_|\___|                                      
%
%   CoBlade(inputFile) executes the Co-Blade code for the input file 
%   named by the string inputFile (with the file extension included).
%
%   If using the compiled version of CoBlade (because you do not have
%   access to the required toolboxes) type the command:
%   CoBlade.exe inputFile
%   at the command prompt, where inputFile is the input file.

%% TO-DO:
% % add diary output
% % add Output folders to seperate data more cleanly? creat output folders for each case
% % change one sided color maps to use grey @ 0 (instead of yellow), red/yellow @ maximum
% % rest the RNG like in HARP_opt
% % actuall write the output files to the Output directory
% % add save for .mat file save (like in HARP_Opt)
% % add outputs to this function

% % vRC_1
% make sure TidalTurbine optimization algorithms are working still (for differen # of CP and different algorithms)
% normalize the fitness value by the initial fitness value
% clean-up all the GMREC Monte Carlo features and add input flags for creating additional plots
% make sure can couple to HARP_Opt as an external dependency

% % vRC_2
% add cross stiffnes terms missing from PreComp
% add Co-Blade to PreComp comparison tool?
% add particle generator for Material Point Method, add VTK particle writer (https://github.com/cfinch/Shocksolution_Examples/blob/master/Visualization/vtktools.py)

% tag_v2.0.0
% verify rc_2 is woooooorking 

% WISHLIST
% non-linear constraints to improve optimization effectiveness
% more physical failure criteria (e.g. Tsai-Wu and etc.)
% easier way to select between materials and regions in GUI feature
% assert statements to take care of all error handling (will reduce LOC significantly I think)


% DOCUMENTATION
% what it __is__, and what it __is not__
% is not a replacement for higher fidelity methods, such as FEM (ANSYS and VABS for example...)

%% clear all variables and close files/figures
% clear all;
close('all');
fclose('all');
format compact;
diary off

if ~isdeployed 
    % Everything in this directory is maintained by the Co-Blade developers
    addpath([pwd filesep 'Source'])
    % These are external dependencies of Co-Blade
    addpath([pwd filesep 'Source' filesep 'bipolar_colormap'])
    addpath([pwd filesep 'Source' filesep 'cbrewer'])
    addpath([pwd filesep 'Source' filesep 'consolidator'])
    addpath([pwd filesep 'Source' filesep 'dispCorrectionFactor'])
    addpath([pwd filesep 'Source' filesep 'export_fig'])
    addpath([pwd filesep 'Source' filesep 'format_ticks'])
    addpath([pwd filesep 'Source' filesep 'polygeom'])
    addpath([pwd filesep 'Source' filesep 'psopt'])
    addpath([pwd filesep 'Source' filesep 'randraw'])
    addpath([pwd filesep 'Source' filesep 'ws2struct'])
end

%% set debugging breakpoints
% dbstop in structOptimize.m at 524
% dbstop in structAnalysis.m at 339
% dbstop in CoBlade_init.m at 95
% dbstop in CoBlade.m at 59
% dbstop in pso.m at 198
% dbstop in plotBestPoint.m at 1
% dbstop in plotBladeShearStress at 69
% dbstop in plotBladeNormStress at 65
% dbstop in plotLaminaStressFC at 28
% dbstop in structOptimize at 250
% dbstop in writeInpFileNewMain at 4

%% run main program
[SIM, ANLS, OPT, ENV, BLADE, WEB, OUT, MATS, AF, Coord] = CoBlade_init(inputFile);


if ~OPT.OPTIMIZE && ANLS.VARY_MAT_PROPS && ANLS.SAMPLE_SIZE > 1
    % Co-Blade is run in stochastic mode
    
    MC = doMonteCarloStuff();

else
    % Co-Blade is run in deterministic mode
    
    % Define the laminate data, either by optimization, or read data from the laminate input files
    if OPT.OPTIMIZE
        WEB.inbStn = OPT.INB_STN .* ones(WEB.NUM_WEBS, 1);  % all the webs begin and end at blade stations INB_STN and OUB_STN
        WEB.oubStn = OPT.OUB_STN .* ones(WEB.NUM_WEBS, 1);
        % Determine the laminate data from optimization routine
        [WEB, SECNODES, LamData] = structOptimize(SIM, ANLS, OPT, ENV, BLADE, WEB, AF, MATS, Coord, OUT);
    else
        % Read the pre-defined laminate data from the laminate input files
        [SECNODES, LamData] = readLaminateData(SIM, BLADE, WEB, MATS);
    end

    % execute the structural analysis
    [Panel, StrProps, AppLoads, ResLoads, Disp, NormS, ShearS, Buckle, MidPlane, LaminaSS, Modes] ...
    = structAnalysis(SIM, ANLS, ENV, BLADE, WEB, MATS, LamData, AF, SECNODES, Coord);
end

                                    

%% write output files
if OPT.OPTIMIZE || OPT.OPT_PITAXIS
    % create a copy of the main input file, but update any parameters that were changed/created by the optimization routine                 
    if OPT.OPTIMIZE
        BLADE.strFile = cell(BLADE.NUM_SEC, 1); % overwrite this variable with the new names
        for n = 1:BLADE.NUM_SEC
            BLADE.strFile{n} = [SIM.case '_OPT_' num2str(n) '.lam'];
        end 
    end
    writeInpFileNewMain(SIM, OPT, BLADE, WEB);
end

if OUT.PROPS_FILE
    writeOupFileProps(SIM, BLADE, StrProps, OUT);
end   

if OUT.LOAD_DSP_FILE
    writeOupFileLoads(SIM, BLADE, AppLoads, ResLoads, Disp, OUT);
end

if OUT.PANEL_FILE
    writeOupFilePanel(SIM, BLADE, WEB, SECNODES, Panel, NormS, ShearS, Buckle, OUT);
end

if OUT.LAMINA_FILE
    writeOupFileLamina(SIM, BLADE, WEB, SECNODES, Panel, Buckle, LaminaSS, OUT);
end

%% create output plots
if OUT.DATA_GUI
    startDataGUI(SIM, OPT, BLADE, WEB, AF, Panel, LaminaSS);
end

if OUT.PLOT_F_BLD || OUT.PLOT_DISP_BLD
    plotLoadsDisp(SIM, BLADE, AF, Coord, AppLoads, Disp, StrProps, OUT)
end

if OUT.PLOT_YMOD
    plotBladeYModulus(SIM, BLADE, WEB, OUT, Panel)
end

if OUT.PLOT_GMOD
    plotBladeGModulus(SIM, BLADE, WEB, OUT, Panel)
end

if OUT.PLOT_MASS_DEN || OUT.PLOT_PRIN_ANG || OUT.PLOT_AT_STFF || ...
   OUT.PLOT_BSTFF    || OUT.PLOT_INER     || OUT.PLOT_CENTERS
    plotStructProps(SIM, BLADE, StrProps, OUT);
end

if OUT.PLOT_NORMS
    plotBladeNormStress(SIM, BLADE, WEB, OUT, Panel, NormS)
end

if OUT.PLOT_SHEARS
    plotBladeShearStress(SIM, BLADE, WEB, OUT, Panel, ShearS)
end

if OUT.PLOT_BCRIT
    plotBladeBuckleCrit(SIM, BLADE, WEB, OUT, Panel, Buckle)
end

if OUT.PLOT_E11 || OUT.PLOT_E22 || OUT.PLOT_E12
    plotLaminaStrain(SIM, BLADE, WEB, MATS, Panel, LaminaSS, OUT)
end

if OUT.PLOT_S11 || OUT.PLOT_S22 || OUT.PLOT_S12
    plotLaminaStress(SIM, BLADE, WEB, MATS, Panel, LaminaSS, OUT)
end

if OUT.PLOT_S11_FC || OUT.PLOT_S22_FC || OUT.PLOT_S12_FC
    plotLaminaStressFC(SIM, BLADE, WEB, MATS, Panel, LaminaSS, OUT)
end

if ( OUT.PLOT_MODE_S || OUT.PLOT_MODE_D ) && ANLS.N_MODES >= 1
    plotModes(SIM, ANLS, Modes, OUT)
end

if OUT.PLOT_APPLOADS || OUT.PLOT_RESLOADS
    plotLoads(SIM, BLADE, AppLoads, ResLoads, OUT)
end

if OUT.PLOT_DEFLECT
    plotDeflect(SIM, BLADE, Disp, OUT)
end

%% Final saving of data and clean-up
CB = ws2struct();   % create a single structure for all the outputs (lots of them!)
save(SIM.matfile);  % save the workspace to a file
fprintf(1, 'All data in Co-Blade workspace saved to MAT-file: \n %s \n\n', SIM.matfile);
fprintf(1, ' Co-Blade v%s terminated normally for case %s. \r\n', SIM.version, SIM.case);
diary OFF
fclose all;
