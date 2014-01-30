function [SIM, ANLS, OPT, ENV, BLADE, WEB, OUT] = readInput_CoBlade(SIM)

%% read the main input file
fid = fopen(SIM.inputFile, 'r');
if fid == -1
    error(['ERROR: Could not locate and open file ' SIM.inputFile]);
end

% Skip through the header of the file
for i = 1:3
    fgetl(fid);
end

% Read the "Develpment Settings" section
fgetl(fid);
SIM.DEBUG_LVL = readScalar(fid);
SIM.RAND_SEED = readScalar(fid);

% Read the "Analysis Options" section
fgetl(fid);
ANLS.SELF_WEIGHT    = readLogical(fid);
ANLS.BUOYANCY       = readLogical(fid);     
ANLS.CENTRIF        = readLogical(fid);    
ANLS.DISP_CF        = readLogical(fid);    
ANLS.N_MODES        = readScalar(fid);
ANLS.N_ELEMS        = readScalar(fid);
ANLS.VARY_MAT_PROPS = readLogical(fid);
ANLS.SAMPLE_SIZE    = readScalar(fid);
ANLS.MAT_COV        = readScalar(fid);

% Read the "Optimization" section
fgetl(fid);
OPT.OPTIMIZE     = readLogical(fid);
OPT.OPT_METHOD   = readString(fid);
OPT.OPT_PITAXIS  = readLogical(fid);
OPT.PITAXIS_VAL  = readScalar(fid);
OPT.INB_STN      = readScalar(fid);
OPT.TRAN_STN     = readScalar(fid);
OPT.OUB_STN      = readScalar(fid);
OPT.NUM_CP       = readScalar(fid);
OPT.READ_INITX   = readLogical(fid);
OPT.INITX_FILE   = readString(fid);
OPT.WRITE_STR    = readLogical(fid);
OPT.WRITE_F_ALL  = readLogical(fid);
OPT.WRITE_X_ALL  = readLogical(fid);
OPT.WRITE_X_ITER = readLogical(fid);

% Read the "Constraints" section
fgetl(fid);
OPT.MAX_TIP_D    = readScalar(fid);
OPT.MIN_FREQ_SEP = readScalar(fid);

% Read the "Environmental Data" section
fgetl(fid);
ENV.FLUID_DEN  = readScalar(fid);
ENV.GRAV       = readScalar(fid);

% Read the "Blade Data" section
fgetl(fid);
BLADE.NUM_SEC     = readScalar(fid);
BLADE.BLD_LENGTH  = readScalar(fid);
BLADE.HUB_RAD     = readScalar(fid);
BLADE.SHAFT_TILT  = readScalar(fid);
BLADE.PRE_CONE    = readScalar(fid);
BLADE.AZIM        = readScalar(fid);
BLADE.BLD_PITCH   = readScalar(fid);
BLADE.ROT_SPD     = readScalar(fid);
BLADE.INTERP_AF   = readString(fid);
BLADE.N_AF        = readScalar(fid);
BLADE.MATS_FILE   = readString(fid);
BLADE.FILLER_DENS = readScalar(fid);

% Read the array of "Blade Data"
fgetl(fid);
fgetl(fid);
bladeDataArray = readCellArray(fid, '%f %f %f %f %f %f %f %q %q', BLADE.NUM_SEC);
BLADE.zFrac     = bladeDataArray{1}; 
BLADE.aeroTwst  = bladeDataArray{2}; 
BLADE.chord     = bladeDataArray{3}; 
BLADE.pitAxis   = bladeDataArray{4};
BLADE.px_a      = bladeDataArray{5}; 
BLADE.py_a      = bladeDataArray{6}; 
BLADE.qz_a      = bladeDataArray{7}; 
BLADE.afFile    = bladeDataArray{8}; 
BLADE.strFile   = bladeDataArray{9}; 

% Read the "Shear Web Data" section
fgetl(fid);
WEB.NUM_WEBS  = readScalar(fid);
WEB.WEB_NODES = readScalar(fid);

if OPT.OPTIMIZE
    % skip the array of "Shear Web Data", it is unused
    % move to the "Output Options" section
    line = [];
    while ~strncmpi(line,'-----  Output Options', 21)
        line = fgetl(fid);
    end

else
    % Read the array of "Shear Web Data"
    fgetl(fid);
    fgetl(fid);
    if WEB.NUM_WEBS >= 1
        webArray     = readCellArray(fid, '%f %f %f %f %f', WEB.NUM_WEBS);
        WEB.webNum   = webArray{1}; 
        WEB.inbStn   = webArray{2};
        WEB.oubStn   = webArray{3};
        WEB.inbChLoc = webArray{4}; 
        WEB.oubChLoc = webArray{5};
    else
        WEB.webNum   = NaN;
        WEB.inbStn   = NaN;
        WEB.oubStn   = NaN;
        WEB.inbChLoc = NaN;
        WEB.oubChLoc = NaN;
    end
    fgetl(fid);

end

% Read the "Output Options" section
OUT.TAB_DEL       = readLogical(fid);
OUT.PROPS_FILE    = readLogical(fid);
OUT.LOAD_DSP_FILE = readLogical(fid);
OUT.PANEL_FILE    = readLogical(fid);
OUT.LAMINA_FILE   = readLogical(fid);
OUT.DATA_GUI      = readLogical(fid);
OUT.SAVE_PLOTS    = readLogical(fid);
OUT.SAVE_FIG_FMT  = readString(fid);
OUT.PLOT_OPT_ITER = readLogical(fid);
OUT.PLOT_F_BLD    = readLogical(fid);
OUT.PLOT_DISP_BLD = readLogical(fid);
OUT.PLOT_GBL_SYS  = readLogical(fid);
OUT.PLOT_YMOD     = readLogical(fid);
OUT.PLOT_GMOD     = readLogical(fid);
OUT.PLOT_MASS_DEN = readLogical(fid);
OUT.PLOT_PRIN_ANG = readLogical(fid);
OUT.PLOT_AT_STFF  = readLogical(fid);
OUT.PLOT_BSTFF    = readLogical(fid);
OUT.PLOT_INER     = readLogical(fid);
OUT.PLOT_CENTERS  = readLogical(fid);                   
OUT.PLOT_NORMS    = readLogical(fid);
OUT.PLOT_SHEARS   = readLogical(fid);
OUT.PLOT_BCRIT    = readLogical(fid);
OUT.PLOT_E11      = readLogical(fid);
OUT.PLOT_E22      = readLogical(fid);
OUT.PLOT_E12      = readLogical(fid);
OUT.PLOT_S11      = readLogical(fid);
OUT.PLOT_S22      = readLogical(fid);
OUT.PLOT_S12      = readLogical(fid);
OUT.PLOT_S11_FC   = readLogical(fid);
OUT.PLOT_S22_FC   = readLogical(fid);
OUT.PLOT_S12_FC   = readLogical(fid);
OUT.PLOT_MODE_D   = readLogical(fid);
OUT.PLOT_MODE_S   = readLogical(fid);
OUT.PLOT_APPLOADS = readLogical(fid);
OUT.PLOT_RESLOADS = readLogical(fid);
OUT.PLOT_DEFLECT  = readLogical(fid);

fclose(fid); % close main input file

