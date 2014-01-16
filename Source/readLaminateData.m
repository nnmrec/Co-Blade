function [SECNODES LamData] = readLaminateData(SIM, BLADE, WEB, MATS)

%% re-assign some structure variable names (for convenience)
NUM_SEC  = BLADE.NUM_SEC;
nWebs    = WEB.nWebs;
NUM_WEBS = WEB.NUM_WEBS;

%%
nSegsTop           = zeros(NUM_SEC, 1);
nSegsBot           = zeros(NUM_SEC, 1); 
xNodeTop           =  cell(NUM_SEC, 1);
xNodeBot           =  cell(NUM_SEC, 1);
LamData(NUM_SEC,1) = struct('Top',[],'Bot',[],'Web',[]);
for i = 1:NUM_SEC
    
    % open the laminate file
    fid = fopen([SIM.laminateDir filesep BLADE.strFile{i}], 'r');
    if fid == -1
        error(['ERROR: Could not locate and open file ' [SIM.laminateDir filesep BLADE.strFile{i}]]);
    end
    
    % move the cursor to the next line that starts with character '*'
    line = [];
    while ~strncmpi(line,'*',1)
        line = fgetl(fid);
    end
    
    % read the "Top Surface" section
    nSegsTop(i) = readScalar(fid);
    fgetl(fid);
    fgetl(fid);
    xNodeTop{i} = readVector(fid,nSegsTop(i)+1)';

    LamData(i).Top.nLam    = zeros(nSegsTop(i), 1); % number of laminas in the laminate
    LamData(i).Top.lamNum  =  cell(nSegsTop(i), 1); % lamina number
    LamData(i).Top.nPlies  =  cell(nSegsTop(i), 1); % number of plies
    LamData(i).Top.tPly    =  cell(nSegsTop(i), 1); % thickness of each ply
    LamData(i).Top.fibAng  =  cell(nSegsTop(i), 1); % fiber angle
    LamData(i).Top.matID   =  cell(nSegsTop(i), 1); % material identification number
    LamData(i).Top.matName =  cell(nSegsTop(i), 1); % material name
    LamData(i).Top.t       = zeros(nSegsTop(i), 1); % total thickness of the laminate
    LamData(i).Top.E_eff   = zeros(nSegsTop(i), 1); % effective Young's modulus of the laminate
    LamData(i).Top.G_eff   = zeros(nSegsTop(i), 1); % effective shear modulus of the laminate
    LamData(i).Top.nu_eff  = zeros(nSegsTop(i), 1); % effective Poisson ratio of the laminate
    LamData(i).Top.massByV = zeros(nSegsTop(i), 1); % effective density of the laminate
    LamData(i).Top.A       =  cell(nSegsTop(i), 1); % extensional stiffness matrix of the laminate
    LamData(i).Top.B       =  cell(nSegsTop(i), 1); % coupling matrix of the laminate
    LamData(i).Top.D       =  cell(nSegsTop(i), 1); % bending stiffness matrix of the laminate
    LamData(i).Top.ABD     =  cell(nSegsTop(i), 1); % ABD matrix of the laminate
    LamData(i).Top.Qbar    =  cell(nSegsTop(i), 1); % transformed reduced stiffness matrix for each lamina in the laminate (a multi-dimensional array)
    LamData(i).Top.T       =  cell(nSegsTop(i), 1); % transformation matrix, which transforms strains/stresses in the x-y axes to the principal axes
    LamData(i).Top.z_i     =  cell(nSegsTop(i), 1); % ply interface locations (z-coordinate in plate coordinate system)
    for n = 1:nSegsTop(i)
        fgetl(fid);
        fgetl(fid);
        header   = readVector(fid,2);
        nLaminas = header(2);
        fgetl(fid);
        fgetl(fid);
        fgetl(fid);
        fgetl(fid);
        lamData  = readCellArray(fid, '%f %f %f %f %f %q', nLaminas);
        
        lamNum  = lamData{1};
        nPlies  = lamData{2};
        tPly    = lamData{3};
        fibAng  = lamData{4};
        matID   = lamData{5};
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
    end
    
    % move the cursor to the next line that starts with character '*'
    line = [];
    while ~strncmpi(line,'*',1)
        line = fgetl(fid);
    end
    % read the "Bottom Surface" section
    nSegsBot(i) = readScalar(fid);
    fgetl(fid);
    fgetl(fid);
    xNodeBot{i} = readVector(fid,nSegsBot(i)+1)';
    
    LamData(i).Bot.nLam    = zeros(nSegsBot(i), 1); 
    LamData(i).Bot.lamNum  =  cell(nSegsBot(i), 1);  
    LamData(i).Bot.nPlies  =  cell(nSegsBot(i), 1);  
    LamData(i).Bot.tPly    =  cell(nSegsBot(i), 1);  
    LamData(i).Bot.fibAng  =  cell(nSegsBot(i), 1); 
    LamData(i).Bot.matID   =  cell(nSegsBot(i), 1);  
    LamData(i).Bot.matName =  cell(nSegsBot(i), 1); 
    LamData(i).Bot.t       = zeros(nSegsBot(i), 1); 
    LamData(i).Bot.E_eff   = zeros(nSegsBot(i), 1); 
    LamData(i).Bot.G_eff   = zeros(nSegsBot(i), 1); 
    LamData(i).Bot.nu_eff  = zeros(nSegsBot(i), 1); 
    LamData(i).Bot.massByV = zeros(nSegsBot(i), 1); 
    LamData(i).Bot.A       =  cell(nSegsBot(i), 1);  
    LamData(i).Bot.B       =  cell(nSegsBot(i), 1);
    LamData(i).Bot.D       =  cell(nSegsBot(i), 1);
    LamData(i).Bot.ABD     =  cell(nSegsBot(i), 1);
    LamData(i).Bot.Qbar    =  cell(nSegsBot(i), 1);
    LamData(i).Bot.T       =  cell(nSegsBot(i), 1);
    LamData(i).Bot.z_i     =  cell(nSegsBot(i), 1);
    for n = 1:nSegsBot(i)
        fgetl(fid);
        fgetl(fid);
        header   = readVector(fid,2);
        nLaminas = header(2);
        fgetl(fid);
        fgetl(fid);
        fgetl(fid);
        fgetl(fid);
        lamData  = readCellArray(fid, '%f %f %f %f %f %q', nLaminas);
        
        lamNum  = lamData{1};
        nPlies  = lamData{2};
        tPly    = lamData{3};
        fibAng  = lamData{4};
        matID   = lamData{5};
        E11     = MATS.E11(matID);
        E22     = MATS.E22(matID);
        G12     = MATS.G12(matID);
        nu12    = MATS.nu12(matID);
        density = MATS.density(matID);
        matName = MATS.matName(matID);
        
        % determine laminate properties by classical lamination theory (CLT)
        CLT = cltAnalysis(lamNum,nPlies,tPly,fibAng,E11,E22,G12,nu12,density);
        LamData(i).Bot.nLam(n)    = length(lamNum);
        LamData(i).Bot.lamNum{n}  = lamNum;
        LamData(i).Bot.nPlies{n}  = nPlies;
        LamData(i).Bot.tPly{n}    = tPly;
        LamData(i).Bot.fibAng{n}  = fibAng;
        LamData(i).Bot.matID{n}   = matID;
        LamData(i).Bot.matName{n} = matName;
        LamData(i).Bot.t(n)       = CLT.t;
        LamData(i).Bot.E_eff(n)   = CLT.E_eff;
        LamData(i).Bot.G_eff(n)   = CLT.G_eff;
        LamData(i).Bot.nu_eff(n)  = CLT.nu_eff;
        LamData(i).Bot.massByV(n) = CLT.massByV;
        LamData(i).Bot.A{n}       = CLT.A;
        LamData(i).Bot.B{n}       = CLT.B;
        LamData(i).Bot.D{n}       = CLT.D;
        LamData(i).Bot.ABD{n}     = CLT.ABD;
        LamData(i).Bot.Qbar{n}    = CLT.Qbar;
        LamData(i).Bot.T{n}       = CLT.T;
        LamData(i).Bot.z_i{n}     = CLT.z_i;
    end 
    
    if nWebs(i) >= 1
        % read the "Web" section
        
        webActive = false(NUM_WEBS, 1); % an array of logicals
        lamNum    =  cell(NUM_WEBS, 1);
        matID     =  cell(NUM_WEBS, 1);
        nPlies    =  cell(NUM_WEBS, 1);
        tPly      =  cell(NUM_WEBS, 1);
        fibAng    =  cell(NUM_WEBS, 1);
        E11       =  cell(NUM_WEBS, 1);
        E22       =  cell(NUM_WEBS, 1);
        G12       =  cell(NUM_WEBS, 1);
        nu12      =  cell(NUM_WEBS, 1);
        density   =  cell(NUM_WEBS, 1);
        matName   =  cell(NUM_WEBS, 1);
        for n = 1:NUM_WEBS
            % move the cursor to the next line that starts with 'web_num'
            line = [];
            while ~strncmpi(line,'web_num',7)
                line = fgetl(fid);
                if line == -1
                    error(['ERROR: Error reading the web data in file ' [SIM.laminateDir filesep BLADE.strFile{i}]]);
                end
            end
        
            header   = readVector(fid,2);
            nLaminas = header(2);
            
            if nLaminas < 1
                webActive(n) = false;
                continue    
            else
                webActive(n) = true;
                
                fgetl(fid);
                fgetl(fid);
                fgetl(fid);
                fgetl(fid);

                lamData = readCellArray(fid, '%f %f %f %f %f %q', nLaminas);
                lamNum{n}  = lamData{1};
                nPlies{n}  = lamData{2};
                tPly{n}    = lamData{3};
                fibAng{n}  = lamData{4};
                matID{n}   = lamData{5};               
                E11{n}     = MATS.E11(matID{n});
                E22{n}     = MATS.E22(matID{n});
                G12{n}     = MATS.G12(matID{n});
                nu12{n}    = MATS.nu12(matID{n});
                density{n} = MATS.density(matID{n});
                matName{n} = MATS.matName(matID{n});
            end           
        end
              
        % only keep the data for the active webs at this station
        lamNum    = lamNum(webActive);
        matID     = matID(webActive);
        nPlies    = nPlies(webActive);
        tPly      = tPly(webActive);
        fibAng    = fibAng(webActive);
        E11       = E11(webActive);
        E22       = E22(webActive);
        G12       = G12(webActive);
        nu12      = nu12(webActive);
        density   = density(webActive);
        matName   = matName(webActive);
        
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
            % determine laminate properties by classical lamination theory (CLT)
            CLT = cltAnalysis(lamNum{n},nPlies{n},tPly{n},fibAng{n},E11{n},E22{n},G12{n},nu12{n},density{n});
            LamData(i).Web.nLam(n)    = length(lamNum{n});
            LamData(i).Web.lamNum{n}  = lamNum{n};
            LamData(i).Web.nPlies{n}  = nPlies{n};
            LamData(i).Web.tPly{n}    = tPly{n};
            LamData(i).Web.fibAng{n}  = fibAng{n};
            LamData(i).Web.matID{n}   = matID{n};
            LamData(i).Web.matName{n} = matName{n};
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
    
    fclose(fid);
end

%% Determine the normalized chord locations of the webs
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
SECNODES.nSegsTop  = nSegsTop;
SECNODES.nSegsBot  = nSegsBot;
SECNODES.xNodeTop  = xNodeTop;
SECNODES.xNodeBot  = xNodeBot;
SECNODES.webLocs   = webLocs;
SECNODES.embNdsTop = embNdsTop;
SECNODES.embNdsBot = embNdsBot;

end % function readLaminateData

