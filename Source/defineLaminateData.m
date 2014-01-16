function [WEB SECNODES LamData] = defineLaminateData(xo, SIM, OPT, BLADE, WEB, MATS, z_oub, z_CP, nSegsTopBot)

%% re-assign some structure variable names (for convenience)
INB_STN  = OPT.INB_STN;
TRAN_STN = OPT.TRAN_STN;
OUB_STN  = OPT.OUB_STN;
NUM_SEC  = BLADE.NUM_SEC;
zSec     = BLADE.zSec;
chord    = BLADE.chord;
pitAxis  = BLADE.pitAxis;
nWebs    = WEB.nWebs;

%% assign the elements of the design vector into meaningful variable names
[xCapSt_inb, ...        
 xCapEnd_inb, ...   
 xCapSt_oub, ...
 xCapEnd_oub, ...
 inbChLoc, ...
 oubChLoc, ...
 t_blade_root, ...
 t_blade_skin, ...
 t_cap_uni, ...
 t_cap_core, ...
 t_lep_core, ...
 t_tep_core, ...
 t_web_skin, ...
 t_web_core] = assignDesignVars(xo, OPT, BLADE, WEB, z_oub, z_CP);

if OPT.WRITE_X_ALL
    % the user has requested to write the design variables to a text file
    fmt = [ repmat('%18.16f  ', 1, numel(xo)), '\r\n' ];
    fid = fopen([SIM.outputDir filesep SIM.case '_allX.out'], 'a');
    fprintf(fid, fmt, xo);       
    fclose(fid);
end
                            
%%           
xCapSt        = planformLine(zSec, chord, pitAxis, xCapSt_inb,  xCapSt_oub,  INB_STN, OUB_STN);
xCapEnd       = planformLine(zSec, chord, pitAxis, xCapEnd_inb, xCapEnd_oub, INB_STN, OUB_STN);
xsec_node_st  = xCapSt  ./ chord + pitAxis; 
xsec_node_end = xCapEnd ./ chord + pitAxis; 

xNodeTop           = cell(NUM_SEC, 1);
xNodeBot           = cell(NUM_SEC, 1);
LamData(NUM_SEC,1) = struct('Top',[],'Bot',[],'Web',[]);
for i = 1:NUM_SEC
    % on the top and bottom, assign the number of segments and chordwise
    % locations of the segment boundaries
    if i <= INB_STN || i > OUB_STN
        % only 1 sector
        xNodeTop{i} = [0; 1];
        xNodeBot{i} = [0; 1];
    else
        % 3 sectors
        xNodeTop{i} = [0; xsec_node_st(i); xsec_node_end(i); 1];
        xNodeBot{i} = [0; xsec_node_st(i); xsec_node_end(i); 1];
    end

    % top surface laminate data
    LamData(i).Top.nLam    = zeros(nSegsTopBot(i), 1); % number of laminas in the laminate
    LamData(i).Top.lamNum  =  cell(nSegsTopBot(i), 1); % lamina number
    LamData(i).Top.nPlies  =  cell(nSegsTopBot(i), 1); % number of plies
    LamData(i).Top.tPly    =  cell(nSegsTopBot(i), 1); % thickness of each ply
    LamData(i).Top.fibAng  =  cell(nSegsTopBot(i), 1); % fiber angle
    LamData(i).Top.matID   =  cell(nSegsTopBot(i), 1); % material identification number
    LamData(i).Top.matName =  cell(nSegsTopBot(i), 1); % material name
    LamData(i).Top.t       = zeros(nSegsTopBot(i), 1); % total thickness of the laminate
    LamData(i).Top.E_eff   = zeros(nSegsTopBot(i), 1); % effective Young's modulus of the laminate
    LamData(i).Top.G_eff   = zeros(nSegsTopBot(i), 1); % effective shear modulus of the laminate
    LamData(i).Top.nu_eff  = zeros(nSegsTopBot(i), 1); % effective Poisson ratio of the laminate
    LamData(i).Top.massByV = zeros(nSegsTopBot(i), 1); % effective density of the laminate
    LamData(i).Top.A       =  cell(nSegsTopBot(i), 1); % extensional stiffness matrix of the laminate
    LamData(i).Top.B       =  cell(nSegsTopBot(i), 1); % coupling matrix of the laminate
    LamData(i).Top.D       =  cell(nSegsTopBot(i), 1); % bending stiffness matrix of the laminate
    LamData(i).Top.ABD     =  cell(nSegsTopBot(i), 1); % ABD matrix of the laminate
    LamData(i).Top.Qbar    =  cell(nSegsTopBot(i), 1); % transformed reduced stiffness matrix for each lamina in the laminate (a multi-dimensional array)
    LamData(i).Top.T       =  cell(nSegsTopBot(i), 1); % transformation matrix, which transforms strains/stresses in the x-y axes to the principal axes
    LamData(i).Top.z_i     =  cell(nSegsTopBot(i), 1); % ply interface locations (z-coordinate in plate coordinate system)   
    for n = 1:nSegsTopBot(i)     
        if     i <= INB_STN
            % blade root build-up
            tPly  = [t_blade_skin(i);
                     t_blade_root(i) * 2;
                     t_blade_skin(i)];
            matID = [2; 1; 2];
            
        elseif i > INB_STN && i < TRAN_STN
            % transition from root build-up
            if     n == 1  % leading edge panel
                tPly  = [t_blade_skin(i);
                         t_blade_root(i);
                         t_lep_core(i) * 2;
                         t_blade_root(i);
                         t_blade_skin(i)];
                matID = [2; 1; 5; 1; 2];
                
            elseif n == 2  % spar cap panel
                tPly  = [t_blade_skin(i);
                         t_blade_root(i);
                         t_cap_uni(i);
                         t_cap_core(i) * 2;
                         t_cap_uni(i);
                         t_blade_root(i);
                         t_blade_skin(i)];
                matID = [2; 1; 3; 4; 3; 1; 2];
                 
            elseif n == 3  % trailing edge panel
                tPly  = [t_blade_skin(i);
                         t_blade_root(i);
                         t_tep_core(i) * 2;
                         t_blade_root(i);
                         t_blade_skin(i)];
                matID = [2; 1; 6; 1; 2];
            end
            
        elseif i >= TRAN_STN && i <= OUB_STN  
            % root material no longer exists
            if     n == 1  % leading edge panel
                tPly  = [t_blade_skin(i);
                         t_lep_core(i) * 2;
                         t_blade_skin(i)];
                matID = [2; 5; 2];
                     
            elseif n == 2  % spar cap panel
                tPly  = [t_blade_skin(i);
                         t_cap_uni(i);
                         t_cap_core(i) * 2;
                         t_cap_uni(i);
                         t_blade_skin(i)];
                matID = [2; 3; 4; 3; 2];
                     
            elseif n == 3  % trailing edge panel
                tPly  = [t_blade_skin(i);
                         t_tep_core(i) * 2;
                         t_blade_skin(i)];
                matID = [2; 6; 2];
                     
            end
            
        elseif i > OUB_STN
            % blade tip
            tPly  = t_blade_skin(i) * 2;
            matID = 2;
            
        end % if i <= INB_STN 
        
        lamNum  = (1:numel(tPly))';
        nPlies  =  ones(numel(tPly), 1);
        fibAng  = zeros(numel(tPly), 1);
        E11     = MATS.E11(matID);
        E22     = MATS.E22(matID);
        G12     = MATS.G12(matID);
        nu12    = MATS.nu12(matID);
        density = MATS.density(matID);
        matName = MATS.matName(matID);

        % determine laminate properties by classical lamination theory (CLT)
        CLT = cltAnalysis(lamNum,nPlies,tPly,fibAng,E11,E22,G12,nu12,density);
        LamData(i).Top.nLam(n)    = length(lamNum);
        LamData(i).Top.lamNum{n}  = lamNum;
        LamData(i).Top.nPlies{n}  = nPlies;
        LamData(i).Top.tPly{n}    = tPly;
        LamData(i).Top.fibAng{n}  = fibAng;
        LamData(i).Top.matID{n}   = matID;
        LamData(i).Top.matName{n} = matName;
        LamData(i).Top.t(n)       = CLT.t;
        LamData(i).Top.E_eff(n)   = CLT.E_eff;
        LamData(i).Top.G_eff(n)   = CLT.G_eff;
        LamData(i).Top.nu_eff(n)  = CLT.nu_eff;
        LamData(i).Top.massByV(n) = CLT.massByV;
        LamData(i).Top.A{n}       = CLT.A;
        LamData(i).Top.B{n}       = CLT.B;
        LamData(i).Top.D{n}       = CLT.D;
        LamData(i).Top.ABD{n}     = CLT.ABD;
        LamData(i).Top.Qbar{n}    = CLT.Qbar;
        LamData(i).Top.T{n}       = CLT.T;
        LamData(i).Top.z_i{n}     = CLT.z_i;
    end % for n = 1:nSegsTopBot(i)
    
    % bottom surface laminate data is same as top in optimization mode
    LamData(i).Bot = LamData(i).Top;
    
    % web laminate data
    if nWebs(i) >= 1

        LamData(i).Web.nLam    = zeros(nWebs(i), 1); 
        LamData(i).Web.lamNum  =  cell(nWebs(i), 1);  
        LamData(i).Web.nPlies  =  cell(nWebs(i), 1);  
        LamData(i).Web.tPly    =  cell(nWebs(i), 1);  
        LamData(i).Web.fibAng  =  cell(nWebs(i), 1); 
        LamData(i).Web.matID   =  cell(nWebs(i), 1);  
        LamData(i).Web.matName =  cell(nWebs(i), 1); 
        LamData(i).Web.t       = zeros(nWebs(i), 1); 
        LamData(i).Web.E_eff   = zeros(nWebs(i), 1); 
        LamData(i).Web.G_eff   = zeros(nWebs(i), 1); 
        LamData(i).Web.nu_eff  = zeros(nWebs(i), 1); 
        LamData(i).Web.massByV = zeros(nWebs(i), 1); 
        LamData(i).Web.A       =  cell(nWebs(i), 1);  
        LamData(i).Web.B       =  cell(nWebs(i), 1);
        LamData(i).Web.D       =  cell(nWebs(i), 1);
        LamData(i).Web.ABD     =  cell(nWebs(i), 1);
        LamData(i).Web.Qbar    =  cell(nWebs(i), 1);
        LamData(i).Web.T       =  cell(nWebs(i), 1);
        LamData(i).Web.z_i     =  cell(nWebs(i), 1);
        for n = 1:nWebs(i)            
            tPly  = [t_web_skin(i);
                     t_web_core(i) * 2;
                     t_web_skin(i)];
            matID = [7; 8; 7];
      
            lamNum  = (1:numel(tPly))';
            nPlies  =  ones(numel(tPly), 1);
            fibAng  = zeros(numel(tPly), 1);
            E11     = MATS.E11(matID);
            E22     = MATS.E22(matID);
            G12     = MATS.G12(matID);
            nu12    = MATS.nu12(matID);
            density = MATS.density(matID);
            matName = MATS.matName(matID);
                         
            % determine laminate properties by classical lamination theory (CLT)
            CLT = cltAnalysis(lamNum,nPlies,tPly,fibAng,E11,E22,G12,nu12,density);
            LamData(i).Web.nLam(n)    = length(lamNum);
            LamData(i).Web.lamNum{n}  = lamNum;
            LamData(i).Web.nPlies{n}  = nPlies;
            LamData(i).Web.tPly{n}    = tPly;
            LamData(i).Web.fibAng{n}  = fibAng;
            LamData(i).Web.matID{n}   = matID;
            LamData(i).Web.matName{n} = matName;
            LamData(i).Web.t(n)       = CLT.t;
            LamData(i).Web.E_eff(n)   = CLT.E_eff;
            LamData(i).Web.G_eff(n)   = CLT.G_eff;
            LamData(i).Web.nu_eff(n)  = CLT.nu_eff;
            LamData(i).Web.massByV(n) = CLT.massByV;
            LamData(i).Web.A{n}       = CLT.A;
            LamData(i).Web.B{n}       = CLT.B;
            LamData(i).Web.D{n}       = CLT.D;
            LamData(i).Web.ABD{n}     = CLT.ABD;
            LamData(i).Web.Qbar{n}    = CLT.Qbar;
            LamData(i).Web.T{n}       = CLT.T;
            LamData(i).Web.z_i{n}     = CLT.z_i;
                                                                      
        end
        
    else
        % create placeholder values
        LamData(i).Web.nLam    = 0;
        LamData(i).Web.lamNum  = [];
        LamData(i).Web.nPlies  = [];
        LamData(i).Web.tPly    = [];
        LamData(i).Web.fibAng  = [];
        LamData(i).Web.matID   = [];
        LamData(i).Web.matName = [];
        LamData(i).Web.t       = [];
        LamData(i).Web.E_eff   = [];
        LamData(i).Web.G_eff   = [];
        LamData(i).Web.nu_eff  = [];
        LamData(i).Web.massByV = [];
        LamData(i).Web.A       = [];
        LamData(i).Web.B       = [];
        LamData(i).Web.D       = [];
        LamData(i).Web.ABD     = [];
        LamData(i).Web.Qbar    = [];
        LamData(i).Web.T       = [];
        LamData(i).Web.z_i     = [];
    end
     
end

%% Collect output 
WEB.inbChLoc = inbChLoc;
WEB.oubChLoc = oubChLoc;

%% Determine the normalized chord locations of the webs and perform error checking
xWebNode = defineWebNodes(BLADE, WEB);
% check to see if the webs extend beyond the leading or trailing edge of the blade
if any(any(xWebNode < 0))
    error('ERROR: some of the shear webs extend beyond the leading edge. x/c < 0');
end
if any(any(xWebNode > 1))
    error('ERROR: some of the shear webs extend beyond the trailing edge. x/c > 1');
end

webLocs   = cell(NUM_SEC, 1);	% normalized chord locations of the shear web nodes
embNdsTop = cell(NUM_SEC, 1);	% normalized chord locations of the top panel endpoints, including the chord locations of the webs
embNdsBot = cell(NUM_SEC, 1);	% normalized chord locations of the bottom panel endpoints, including the chord locations of the webs
for i = 1:NUM_SEC
    % web nodes
    wl           = xWebNode(i,:);
    webLocs{i}   = wl(~isnan(wl))'; % an empty array if no webs exist at this station
    % top of the airfoil
    embNdsTop{i} = consolidator([xNodeTop{i}; webLocs{i}], [] ,[], 1e-4);  	% a sorted array in ascending order
    % bottom of the airfoil
    embNdsBot{i} = consolidator([xNodeBot{i}; webLocs{i}], [] ,[], 1e-4);	% a sorted array in ascending order  
end

%% Collect output 
SECNODES.nSegsTop  = nSegsTopBot;
SECNODES.nSegsBot  = nSegsTopBot;
SECNODES.xNodeTop  = xNodeTop;
SECNODES.xNodeBot  = xNodeBot;
SECNODES.webLocs   = webLocs;
SECNODES.embNdsTop = embNdsTop;
SECNODES.embNdsBot = embNdsBot;

end % function defineLaminateData
