function [Panel StrProps AppLoads ResLoads Disp NormS ShearS Buckle MidPlane LaminaSS Modes] ...
         = structAnalysis(SIM, ANLS, ENV, BLADE, WEB, MATS, LamData, AF, SECNODES, Coord)   
       
%% re-assign some structure variable names (for convenience) NOTE: the FEX has some more automated ways to do this
SELF_WEIGHT  = ANLS.SELF_WEIGHT;
BUOYANCY     = ANLS.BUOYANCY;
CENTRIF      = ANLS.CENTRIF;
DISP_CF      = ANLS.DISP_CF;
N_MODES      = ANLS.N_MODES;
FLUID_DEN    = ENV.FLUID_DEN;
GRAV         = ENV.GRAV;
NUM_SEC      = BLADE.NUM_SEC;
HUB_RAD      = BLADE.HUB_RAD;
PRE_CONE     = BLADE.PRE_CONE;
ROT_SPD      = BLADE.ROT_SPD;
zFrac        = BLADE.zFrac;
zSec         = BLADE.zSec;
aeroTwst     = BLADE.aeroTwst;
chord        = BLADE.chord;
pitAxis      = BLADE.pitAxis;
px_a         = BLADE.px_a;
py_a         = BLADE.py_a;
qz_a         = BLADE.qz_a;
nCells       = BLADE.nCells;
WEB_NODES    = WEB.WEB_NODES;
nWebs        = WEB.nWebs;
NormAFcoords = AF.NormAFcoords;
nSegsTop     = SECNODES.nSegsTop;
nSegsBot     = SECNODES.nSegsBot;
xNodeTop     = SECNODES.xNodeTop;
xNodeBot     = SECNODES.xNodeBot;
webLocs      = SECNODES.webLocs;
embNdsTop    = SECNODES.embNdsTop;
embNdsBot    = SECNODES.embNdsBot;
Tran         = Coord.Tran;

%%
d2r = pi/180; % convert degrees to radians
                      
%% Create data structures for the panel geometry and properties
Panel(NUM_SEC,1) = struct('Top',[],'Bot',[],'Web',[]);  % structure array storing various data related to the panel geometry and properties   
xw_top           =  cell(NUM_SEC, 1);                   % x-coordinates of the middle of the shear webs on the top surface
yw_top           =  cell(NUM_SEC, 1);                   % y-coordinates of the middle of the shear webs on the top surface
xw_bot           =  cell(NUM_SEC, 1);                   % x-coordinates of the middle of the shear webs on the bottom surface
yw_bot           =  cell(NUM_SEC, 1);                   % y-coordinates of the middle of the shear webs on the bottom surface
A_disp           = zeros(NUM_SEC, 1);                   % area displaced by the dimensional airfoil shape
x_mc             = zeros(NUM_SEC, 1);                   % x-coordinate of the mid-chord point w.r.t the x-y blade reference axes
y_mc             = zeros(NUM_SEC, 1);                   % y-coordinate of the mid-chord point w.r.t the x-y blade reference axes
for i = 1:NUM_SEC
    % airfoil node coordinates 
    x             = NormAFcoords(i).x;  % normalized x-coordinate of airfoil nodes
    y             = NormAFcoords(i).y;  % normalized y-coordinate of airfoil nodes
    [unused i_TE] = min(abs(x - 1));    % index of the trailing edge airfoil node
    
    % top of airfoil
    x_top         = consolidator([x(1:i_TE); embNdsTop{i}], [] ,[], 1e-4);   % a sorted array in ascending order
    y_top         = interp1(x(1:i_TE), y(1:i_TE), x_top);
    i_sNodesTop   = findClosest(x_top, xNodeTop{i});	% indices of the sector nodes on the airfoil top 
    i_wNodesTop   = findClosest(x_top, webLocs{i});   	% indices of the web nodes on the airfoil top 
    
    % bottom of the airfoil
    x_bot        = consolidator([flipud([x(i_TE:end); 0]); embNdsBot{i}], [] ,[], 1e-4);	% a sorted array in ascending order
    y_bot        = interp1([x(i_TE:end); 0], [y(i_TE:end); 0], x_bot);
    i_sNodesBot  = findClosest(x_bot, xNodeBot{i});     	% indices of the sector nodes on the airfoil bottom 
    i_wNodesBot  = findClosest(x_bot, webLocs{i});          % indices of the web nodes on the airfoil bottom 
        
    % Now scale, translate, and rotate the all coordinates
    [x_top y_top] = scaleAirfoilCoords(x_top, y_top, pitAxis(i), chord(i), aeroTwst(i));
    [x_bot y_bot] = scaleAirfoilCoords(x_bot, y_bot, pitAxis(i), chord(i), aeroTwst(i));
    
    % Compute the mid-chord point
    x_mc(i) = (x_top(1) + x_top(end)) / 2;
    y_mc(i) = (y_top(1) + y_top(end)) / 2;
    
    % endpoints of the shear webs
    xw_top{i} = x_top(i_wNodesTop);
    yw_top{i} = y_top(i_wNodesTop);
    xw_bot{i} = x_bot(i_wNodesBot);
    yw_bot{i} = y_bot(i_wNodesBot);
        
    % Compute the area displaced by the dimensional airfoil shape  
    A_disp(i) = polyarea([x_top; flipud(x_bot)], [y_top; flipud(y_bot)]);
    
    % Define the panel data structures 
    Panel(i).Top  = definePanel('top',    x_top, y_top, i_sNodesTop, i_wNodesTop, zSec(i), LamData(i).Top);
    Panel(i).Bot  = definePanel('bottom', x_bot, y_bot, i_sNodesBot, i_wNodesBot, zSec(i), LamData(i).Bot);   
    [x_web y_web] = webCoords(Panel(i), xw_top{i}, yw_top{i}, xw_bot{i}, yw_bot{i}, nWebs(i), WEB_NODES); % computes the x-y coordinates of the web mid-lines                     
    Panel(i).Web  = definePanel('web',    x_web, y_web,          [],    nWebs(i), zSec(i), LamData(i).Web);
                        
    % plot to verify location of nodes and panels (for debugging only)
%     figure
%     hold on
%     plot(x_top, y_top, '.-k')
%     plot(x_top(i_sNodesTop), y_top(i_sNodesTop), 'sr')
%     plot(x_top(i_wNodesTop), y_top(i_wNodesTop), 'or')
%     plot(x_bot, y_bot, '.-r')
%     plot(x_bot(i_sNodesBot), y_bot(i_sNodesBot), 'xb')
%     plot(x_bot(i_wNodesBot), y_bot(i_wNodesBot), 'ob')
%     legend('top: af nodes','top: sect nodes','top: web nodes', ...
%            'bot: af nodes','bot: sect nodes','bot: web nodes')
%     for n = 1:Panel(i).Top.nPanels
%         plot3(Panel(i).Top.x{n},  Panel(i).Top.y{n},  Panel(i).Top.z{n},  '.-r')
%         plot3(Panel(i).Top.xs{n}, Panel(i).Top.ys{n}, Panel(i).Top.zs{n}, '.:g')
%     end
%     for n = 1:Panel(i).Bot.nPanels
%         plot3(Panel(i).Bot.x{n},  Panel(i).Bot.y{n},  Panel(i).Bot.z{n},  '.-b')
%         plot3(Panel(i).Bot.xs{n}, Panel(i).Bot.ys{n}, Panel(i).Bot.zs{n}, '.:g')
%     end
%     for n = 1:Panel(i).Web.nPanels
%         plot3(Panel(i).Web.x{n},  Panel(i).Web.y{n},  Panel(i).Web.z{n},  '.-k')
%         plot3(Panel(i).Web.xs{n}, Panel(i).Web.ys{n}, Panel(i).Web.zs{n}, '.:g')
%     end
%     axis equal
    
end

%% Compute structural properties w.r.t. the reference x-y axes
E_ref      = zeros(NUM_SEC,1);    % reference Young's modulus
A_crs      = zeros(NUM_SEC,1);    % geometric cross sectional area
A_e        = zeros(NUM_SEC,1);    % Young's modulus weighted effective cross sectional area
axial_stff = zeros(NUM_SEC,1);    % Young's modulus weighted effective axial stiffness
Ix         = zeros(NUM_SEC,1);    % second area moment of inertia w.r.t. the reference x-axis of the blade coord. sys.
Iy         = zeros(NUM_SEC,1);    % second area moment of inertia w.r.t. the reference y-axis of the blade coord. sys.
Ixy        = zeros(NUM_SEC,1);    % product area moment of inertia w.r.t. the reference x-y axes of the blade coord. sys.
Ix_e       = zeros(NUM_SEC,1);    % Young's modulus weighted effective second area moment of inertia w.r.t. the reference x-axis of the blade coord. sys.
Iy_e       = zeros(NUM_SEC,1);    % Young's modulus weighted effective second area moment of inertia w.r.t. the reference y-axis of the blade coord. sys.
Ixy_e      = zeros(NUM_SEC,1);    % Young's modulus weighted effective product area moment of inertia w.r.t. the reference x-y axes of the blade coord. sys.
mIx        = zeros(NUM_SEC,1);    % mass moment of inertia per unit length (span) w.r.t. the reference x-axis of the blade coord. sys.
mIy        = zeros(NUM_SEC,1);    % mass moment of inertia per unit length (span) w.r.t. the reference y-axis of the blade coord. sys.
mIxy       = zeros(NUM_SEC,1);    % product mass moment of inertia per unit length (span) w.r.t. the reference x-y axes of the blade coord. sys.
mass_den   = zeros(NUM_SEC,1);    % mass per unit length (span)
x_cm       = zeros(NUM_SEC,1);    % x-coordinate center of mass w.r.t. the reference blade coord. sys.
y_cm       = zeros(NUM_SEC,1);    % y-coordinate center of mass w.r.t. the reference blade coord. sys.
x_tc       = zeros(NUM_SEC,1);    % x-coordinate of Young's modulus weighted effective centroid (tension center) w.r.t. the reference blade coord. sys.
y_tc       = zeros(NUM_SEC,1);    % y-coordinate of Young's modulus weighted effective centroid (tension center) w.r.t. the reference blade coord. sys.
for i = 1:NUM_SEC
    % create vectors of the panel properties (we will perform summation of these vectors to compute effective cross section properties)
    crsArea = [Panel(i).Top.crsArea; Panel(i).Bot.crsArea; Panel(i).Web.crsArea];
    E_eff   = [Panel(i).Top.E_eff;   Panel(i).Bot.E_eff;   Panel(i).Web.E_eff];
    I_x     = [Panel(i).Top.Ix;      Panel(i).Bot.Ix;      Panel(i).Web.Ix];
    I_y     = [Panel(i).Top.Iy;      Panel(i).Bot.Iy;      Panel(i).Web.Iy];
    I_xy    = [Panel(i).Top.Ixy;     Panel(i).Bot.Ixy;     Panel(i).Web.Ixy];
    massByL = [Panel(i).Top.massByL; Panel(i).Bot.massByL; Panel(i).Web.massByL];
    massByV = [Panel(i).Top.massByV; Panel(i).Bot.massByV; Panel(i).Web.massByV];
    x_cen   = [Panel(i).Top.x_cen;   Panel(i).Bot.x_cen;   Panel(i).Web.x_cen];
    y_cen   = [Panel(i).Top.y_cen;   Panel(i).Bot.y_cen;   Panel(i).Web.y_cen];
    
    % reference effective Young's modulus
    E_ref(i) = max(E_eff);    % can choose any value that exists in the current cross section, just need to be consistent
    
    % cross sectional area
    A_crs(i) = sum(crsArea);
        
    % Young's modulus weighted effective cross sectional area
    A_e(i) = sum( E_eff.*crsArea ) / E_ref(i);
    
    % Young's modulus weighted effective axial stiffness
    axial_stff(i) = E_ref(i)*A_e(i);
    
    % second area moments of inertia w.r.t the x-y reference axes
    Ix(i)  = sum(I_x);
    Iy(i)  = sum(I_y);
    Ixy(i) = sum(I_xy);
     
    % Young's modulus weighted effective second area moments of inertia w.r.t. the x-y reference axes
    Ix_e(i)  = sum( E_eff.*I_x)  / E_ref(i);
    Iy_e(i)  = sum( E_eff.*I_y)  / E_ref(i);
    Ixy_e(i) = sum( E_eff.*I_xy) / E_ref(i);
    
    % mass moments of inertia per unit length (span) w.r.t the x-y reference axes
    mIx(i)  = sum(massByV.*I_x);
    mIy(i)  = sum(massByV.*I_y);
    mIxy(i) = sum(massByV.*I_xy);
    
    % mass per unit blade length
    mass_den(i) = sum(massByL);
    
    % center of mass
    x_cm(i) = sum(massByL .* x_cen) / mass_den(i);
    y_cm(i) = sum(massByL .* y_cen) / mass_den(i);
    
    % Young's modulus weighted effective centroid (tension center) 
    x_tc(i) = sum( E_eff.*crsArea.*x_cen ) / (E_ref(i)*A_e(i));
    y_tc(i) = sum( E_eff.*crsArea.*y_cen ) / (E_ref(i)*A_e(i));       
end

% mass moments of inertia per unit length (span) w.r.t. the inertial (center of mass) axes of the blade coord. sys.
mIx_cm  = mIx  - mass_den.*(y_cm.^2);   
mIy_cm  = mIy  - mass_den.*(x_cm.^2);   
mIxy_cm = mIxy - mass_den.*x_cm.*y_cm;  
  
% mass moments of inertia per unit length (span) w.r.t. the centroidal (tension center) axes of the blade coord. sys.
mIx_tc  = mIx_cm  + mass_den.*((y_tc - y_cm).^2);           
mIy_tc  = mIy_cm  + mass_den.*((x_tc - x_cm).^2);           
mIxy_tc = mIxy_cm + mass_den.*(x_tc - x_cm).*(y_tc - y_cm);		

% mass moments of inertia per unit length (span) w.r.t. the principal inertial (center of mass) axes
% flapIner_cm : flapwise mass moment of inertia per unit length (span) w.r.t. the principal inertial (mass center) x-y axes of the blade coord. sys.
% edgeIner_cm : edgewise mass moment of inertia per unit length (span) w.r.t. the principal inertial (mass center) x-y axes of the blade coord. sys.
% iner_tw     : inertial twist, the angle between the reference x-axis and the inertial 1st principal (flapwise) axis
[flapIner_cm edgeIner_cm iner_tw] = principalValues(mIx_cm, mIy_cm, mIxy_cm, aeroTwst);
    
% mass moments of inertia per unit length (span) w.r.t. the principal centroidal (tension center) axes
% flapIner_tc = zeros(NUM_SEC,1);    % flapwise mass moment of inertia per unit length (span) w.r.t. the principal centroidal (tension center) x-y axes of the blade coord. sys.
% edgeIner_tc = zeros(NUM_SEC,1);    % edgewise mass moment of inertia per unit length (span) w.r.t. the principal centroidal (tension center) x-y axes of the blade coord. sys.
[flapIner_tc edgeIner_tc] = principalValues(mIx_tc, mIy_tc, mIxy_tc, aeroTwst);

% Young's modulus weighted effective second area moments of inertia w.r.t. the inertial (center of mass) axes of the blade coord. sys.
Ix_e_cm  = Ix_e  - A_e.*(y_cm.^2);  
Iy_e_cm  = Iy_e  - A_e.*(x_cm.^2);  
Ixy_e_cm = Ixy_e - A_e.*x_cm.*y_cm; 

% Young's modulus weighted effective second area moments of inertia w.r.t. the centroidal (tension center) axes of the blade coord. sys.
% Ix_e_tc  = Ix_e_cm  + A_e.*((y_tc - y_cm).^2);            % this was a mistake!    
% Iy_e_tc  = Iy_e_cm  + A_e.*((x_tc - x_cm).^2);            % this was a mistake!      
% Ixy_e_tc = Ixy_e_cm + A_e.*(x_tc - x_cm).*(y_tc - y_cm);  % this was a mistake!
Ix_e_tc  = Ix_e  - A_e.*(y_tc.^2);              
Iy_e_tc  = Iy_e  - A_e.*(x_tc.^2);              
Ixy_e_tc = Ixy_e - A_e.*x_tc.*y_tc;

% Young's modulus weighted effective second area moments of inertia w.r.t. the mid-chord axes of the blade coord. sys. (this is used for the displacement correction factors)
% Ix_e_mc  = Ix_e_cm  + A_e.*((y_mc - y_cm).^2);            % this was a mistake!         
% Iy_e_mc  = Iy_e_cm  + A_e.*((x_mc - x_cm).^2);            % this was a mistake!          
% Ixy_e_mc = Ixy_e_cm + A_e.*(x_mc - x_cm).*(y_mc - y_cm);  % this was a mistake!    
Ix_e_mc  = Ix_e_tc  + A_e.*((y_tc - y_mc).^2);              
Iy_e_mc  = Iy_e_tc  + A_e.*((x_tc - x_mc).^2);              
Ixy_e_mc = Ixy_e_tc + A_e.*(x_tc - x_mc).*(y_tc - y_mc);

% bending stiffnesses w.r.t. the reference x-y axes of the blade coord. sys.
EIx  = E_ref.*Ix_e;     
EIy  = E_ref.*Iy_e;    
EIxy = E_ref.*Ixy_e;

% bending stiffnesses w.r.t. the inertial (center of mass) axes of the blade coord. sys.
EIx_cm  = E_ref.*Ix_e_cm;     
EIy_cm  = E_ref.*Iy_e_cm;     
EIxy_cm = E_ref.*Ixy_e_cm;   

% bending stiffnesses w.r.t. the centroidal (tension center) axes of the blade coord. sys.
EIx_tc  = E_ref.*Ix_e_tc;    
EIy_tc  = E_ref.*Iy_e_tc;    
EIxy_tc = E_ref.*Ixy_e_tc;   

% principal Young's modulus weighted effective second area moments of inertia w.r.t. the mid-chord axes (this is used for the displacement correction factors)
% flapI_mc : flapwise bending stiffness w.r.t. the principal inertial (center of mass) x-y axes of the blade coord. sys.
% edgeI_mc : edgewise bending stiffness w.r.t. the principal inertial (center of mass) x-y axes of the blade coord. sys.
[flapI_mc edgeI_mc] = principalValues(Ix_e_mc, Iy_e_mc, Ixy_e_mc, aeroTwst);

% bending stiffnesses w.r.t. the principal inertial (center of mass) axes
% flapEI_cm : flapwise bending stiffness w.r.t. the principal inertial (center of mass) x-y axes of the blade coord. sys.
% edgeEI_cm : edgewise bending stiffness w.r.t. the principal inertial (center of mass) x-y axes of the blade coord. sys.
[flapEI_cm edgeEI_cm] = principalValues(EIx_cm, EIy_cm, EIxy_cm, aeroTwst);

% bending stiffnesses w.r.t. the principal centroidal (tension center) axes
% flapEI_tc : flapwise bending stiffness w.r.t. the principal centroidal (tension center) x-y axes of the blade coord. sys.
% edgeEI_tc : edgewise bending stiffness w.r.t. the principal centroidal (tension center) x-y axes of the blade coord. sys.
% cent_tw   : centroidal twist, the angle between the reference x-axis and the centroidal 1st principal (flapwise) axis
[flapEI_tc edgeEI_tc cent_tw] = principalValues(EIx_tc, EIy_tc, EIxy_tc, aeroTwst);

%% Pre-process some of the panel data prior to calculating the elastic properties, which depend on shear flow
[Top Bot Web Cell] = preProcessShearFlowData(NUM_SEC, x_tc, y_tc, Panel, xw_top, yw_top, xw_bot, yw_bot, nCells);                                        

%% Compute structural properties depending on shear flow                                     
x_sc     = zeros(NUM_SEC,1);    % x-coordinate of shear center w.r.t. the reference blade coord. sys.
y_sc     = zeros(NUM_SEC,1);    % y-coordinate of shear center w.r.t. the reference blade coord. sys.
tor_stff = zeros(NUM_SEC,1);    % torsional stiffness
for i = 1:NUM_SEC
    % shear center
    [x_sc(i) y_sc(i)] = shearCenter(Panel(i), ...
                                    Top(i), ...
                                    Bot(i), ...
                                    Web(i), ...
                                    Cell(i), ...
                                    nSegsTop(i), ...
                                    nSegsBot(i), ...
                                    nCells(i), ...
                                    x_cm(i), ...
                                    y_cm(i), ...
                                    x_tc(i), ...
                                    y_tc(i), ...
                                    axial_stff(i), ...
                                    EIx_tc(i), ...
                                    EIy_tc(i), ...
                                    EIxy_tc(i));

    % torsional stiffness     
    [unused tor_stff(i)] = stressShear(Top(i), ...
                                       Bot(i), ...
                                       Web(i), ...
                                       Cell(i), ...
                                       Panel(i).Top.nPanels, ...
                                       Panel(i).Bot.nPanels, ...
                                       nCells(i), ...
                                       x_cm(i), ...
                                       y_cm(i), ...
                                       x_tc(i), ...
                                       y_tc(i), ...
                                       axial_stff(i), ...
                                       EIx_tc(i), ...
                                       EIy_tc(i), ...
                                       EIxy_tc(i), ...
                                       0, ...
                                       0, ...
                                       1, ...
                                       0, ...
                                       0);      
end

% mass moments of inertia per unit length (span) w.r.t. the elastic (shear center) axes of the blade coord. sys.
mIx_sc  = mIx_cm  + mass_den.*((y_sc - y_cm).^2);          
mIy_sc  = mIy_cm  + mass_den.*((x_sc - x_cm).^2);           
mIxy_sc = mIxy_cm + mass_den.*(x_sc - x_cm).*(y_sc - y_cm);	

% mass moments of inertia per unit length (span) w.r.t. the principal elastic (shear center) axes
% flapIner_sc : flapwise mass moment of inertia per unit length (span) w.r.t. the principal elastic (shear center) x-y axes of the blade coord. sys.
% edgeIner_sc : edgewise mass moment of inertia per unit length (span) w.r.t. the principal elastic (shear center) x-y axes of the blade coord. sys.
[flapIner_sc edgeIner_sc] = principalValues(mIx_sc, mIy_sc, mIxy_sc, aeroTwst);

% Young's modulus weighted effective second area moments of inertia w.r.t. the elastic (shear center) axes of the blade coord. sys.
% Ix_e_sc  = Ix_e_cm  + A_e.*((y_sc - y_cm).^2);            % this was a mistake!           
% Iy_e_sc  = Iy_e_cm  + A_e.*((x_sc - x_cm).^2);            % this was a mistake!              
% Ixy_e_sc = Ixy_e_cm + A_e.*(x_sc - x_cm).*(y_sc - y_cm);  % this was a mistake!    
% Ix_e_sc  = Ix_e  - A_e.*(y_sc.^2);                        % this was a mistake!
% Iy_e_sc  = Iy_e  - A_e.*(x_sc.^2);                        % this was a mistake!
% Ixy_e_sc = Ixy_e - A_e.*x_sc.*y_sc;                       % this was a mistake!
Ix_e_sc  = Ix_e_tc  + A_e.*((y_tc - y_sc).^2);              
Iy_e_sc  = Iy_e_tc  + A_e.*((x_tc - x_sc).^2);              
Ixy_e_sc = Ixy_e_tc + A_e.*(x_tc - x_sc).*(y_tc - y_sc);

% bending stiffnesses w.r.t. the elastic (shear center) axes of the blade coord. sys.
EIx_sc  = E_ref.*Ix_e_sc;     
EIy_sc  = E_ref.*Iy_e_sc;     
EIxy_sc = E_ref.*Ixy_e_sc;    

% bending stiffnesses w.r.t. the principal elastic (shear center) axes
% flapEI_sc : flapwise bending stiffness w.r.t. the principal elastic (shear center) x-y axes of the blade coord. sys.
% edgeEI_sc : edgewise bending stiffness w.r.t. the principal elastic (shear center) x-y axes of the blade coord. sys.
% elas_tw   : elastic twist, the angle between the reference x-axis and the elastic 1st principal (flapwise) axis
[flapEI_sc edgeEI_sc elas_tw] = principalValues(EIx_sc, EIy_sc, EIxy_sc, aeroTwst);

%% Chordwise offsets for the center of mass, tension center, and shear center
thetaC   = atan2(y_cm, x_cm);
thetaT   = atan2(y_tc, x_tc);
thetaS   = atan2(y_sc, x_sc);
rC       = sqrt( x_cm.^2 + y_cm.^2 );
rT       = sqrt( x_tc.^2 + y_tc.^2 );
rS       = sqrt( x_sc.^2 + y_sc.^2 );
cm_offst = rC.*cos(thetaC - aeroTwst.*d2r) ./ chord; % center of mass chordwise offset from the pitch axis
tc_offst = rT.*cos(thetaT - aeroTwst.*d2r) ./ chord; % tension center chordwise offset from the pitch axis
sc_offst = rS.*cos(thetaS - aeroTwst.*d2r) ./ chord; % shear center chordwise offset from the pitch axis

%% compute the total blade mass
bladeMass = trapzf(zSec, mass_den);

%% For convenience, store some of the strucutral properties in a strucuture
StrProps.mass_den    = mass_den;
StrProps.iner_tw     = iner_tw;
StrProps.cent_tw     = cent_tw;
StrProps.elas_tw     = elas_tw;	
StrProps.axial_stff  = axial_stff;	
StrProps.tor_stff    = tor_stff;	
StrProps.EIx         = EIx;	
StrProps.EIy         = EIy;	
StrProps.EIxy        = EIxy;	
StrProps.flapEI_cm   = flapEI_cm;	
StrProps.edgeEI_cm   = edgeEI_cm;	
StrProps.flapEI_tc   = flapEI_tc;	
StrProps.edgeEI_tc   = edgeEI_tc;	
StrProps.flapEI_sc   = flapEI_sc;	
StrProps.edgeEI_sc   = edgeEI_sc;	
StrProps.mIx         = mIx;	
StrProps.mIy         = mIy;	
StrProps.mIxy        = mIxy;	
StrProps.flapIner_cm = flapIner_cm;	
StrProps.edgeIner_cm = edgeIner_cm;	
StrProps.flapIner_tc = flapIner_tc;	
StrProps.edgeIner_tc = edgeIner_tc;
StrProps.flapIner_sc = flapIner_sc;
StrProps.edgeIner_sc = edgeIner_sc;
StrProps.cm_offst    = cm_offst;	
StrProps.tc_offst    = tc_offst;	
StrProps.sc_offst    = sc_offst;	
StrProps.x_cm        = x_cm;	
StrProps.y_cm        = y_cm;	
StrProps.x_tc        = x_tc;	
StrProps.y_tc        = y_tc;	
StrProps.x_sc        = x_sc;
StrProps.y_sc        = y_sc;
StrProps.bladeMass   = bladeMass;

%% Compute applied forces & moments, and the resultant forces & moments
% applied loads    
AppLoads = loadsApplied(zSec, ...
                        NUM_SEC, ...
                        HUB_RAD, ...
                        PRE_CONE, ...
                        ROT_SPD, ...
                        FLUID_DEN, ...
                        GRAV, ...
                        mass_den, ...
                        A_disp, ...
                        px_a, ...
                        py_a, ...
                        qz_a, ...
                        Tran, ...
                        SELF_WEIGHT, ...
                        BUOYANCY, ...
                        CENTRIF);
% resultant loads
ResLoads = loadsResultant(zSec, ...
                          NUM_SEC, ...
                          AppLoads,...
                          x_cm, ...
                          y_cm, ...
                          x_tc, ...
                          y_tc, ...
                          x_sc, ...
                          y_sc);

%% Compute the blade displacements
Disp = bladeDisplacement(zSec, ...
                         zFrac, ...
                         NUM_SEC, ...
                         ResLoads, ...
                         flapI_mc, ...
                         edgeI_mc, ...
                         EIx_tc, ...
                         EIy_tc, ...
                         EIxy_tc, ...
                         axial_stff, ...
                         tor_stff, ...
                         Tran, ...
                         DISP_CF);
   
%% Compute the strains and stresses
NormS(NUM_SEC, 1)    = struct('Top',[],'Bot',[],'Web',[]);  % structure array storing data related to the normal stresses in the panels
ShearS(NUM_SEC, 1)   = struct('Top',[],'Bot',[],'Web',[]);  % structure array storing data related to the shear stresses in the panels
Buckle(NUM_SEC, 1)   = struct('Top',[],'Bot',[],'Web',[]);  % structure array storing data related to the buckling criteria of the panels
MidPlane(NUM_SEC, 1) = struct('Top',[],'Bot',[],'Web',[]);  % structure array storing data related to the midplane strains and curvatures of panel laminates
LaminaSS(NUM_SEC, 1) = struct('Top',[],'Bot',[],'Web',[]);  % structure array storing data related to the lamina strains and stresses of the panel laminates
for i = 1:NUM_SEC
    % normal stresses
    NormS(i).Top = stressNormal(Panel(i).Top, x_tc(i), y_tc(i), axial_stff(i), EIx_tc(i), EIy_tc(i), EIxy_tc(i), ResLoads.Vz(i), ResLoads.Mx(i), ResLoads.My(i));
    NormS(i).Bot = stressNormal(Panel(i).Bot, x_tc(i), y_tc(i), axial_stff(i), EIx_tc(i), EIy_tc(i), EIxy_tc(i), ResLoads.Vz(i), ResLoads.Mx(i), ResLoads.My(i));
    NormS(i).Web = stressNormal(Panel(i).Web, x_tc(i), y_tc(i), axial_stff(i), EIx_tc(i), EIy_tc(i), EIxy_tc(i), ResLoads.Vz(i), ResLoads.Mx(i), ResLoads.My(i));
                            
    % shear stresses
    ShearS(i) = stressShear(Top(i), ...
                            Bot(i), ...
                            Web(i), ...
                            Cell(i), ...
                            Panel(i).Top.nPanels, ...
                            Panel(i).Bot.nPanels, ...
                            nCells(i), ...
                            x_cm(i), ...
                            y_cm(i), ...
                            x_tc(i), ...
                            y_tc(i), ...
                            axial_stff(i), ...
                            EIx_tc(i), ...
                            EIy_tc(i), ...
                            EIxy_tc(i), ...
                            ResLoads.Vx(i), ...
                            ResLoads.Vy(i), ...
                            ResLoads.Mz(i), ...
                            AppLoads.pz_w(i), ...
                            AppLoads.pz_c(i));
    
    % Compute panel buckling criterias
    Buckle(i).Top = stressBuckle(Panel(i).Top, NormS(i).Top, ShearS(i).Top, 'top');
    Buckle(i).Bot = stressBuckle(Panel(i).Bot, NormS(i).Bot, ShearS(i).Bot, 'bottom');
    Buckle(i).Web = stressBuckle(Panel(i).Web, NormS(i).Web, ShearS(i).Web, 'web');
    
	% plate midplane strains and curvatures
    MidPlane(i).Top = strainPlate(Panel(i).Top, NormS(i).Top, ShearS(i).Top);
    MidPlane(i).Bot = strainPlate(Panel(i).Bot, NormS(i).Bot, ShearS(i).Bot);
    MidPlane(i).Web = strainPlate(Panel(i).Web, NormS(i).Web, ShearS(i).Web);

    % lamina principal stresses and principal strains
    LaminaSS(i).Top = strainStressPly(Panel(i).Top, MidPlane(i).Top, MATS);
    LaminaSS(i).Bot = strainStressPly(Panel(i).Bot, MidPlane(i).Bot, MATS);
    LaminaSS(i).Web = strainStressPly(Panel(i).Web, MidPlane(i).Web, MATS);

end              

% check for errors on the inputs first though, this is an indicator if
% something failed during the structural analysis  NOTE: it would be much
% better to use ASSERT statements so we dont have to carry around error
% checking code everywhere...
data = [mass_den, ...
        iner_tw, ... 
        cent_tw, ... 
        elas_tw, ... 
        axial_stff, ...
        tor_stff, ...  
        EIx, ...  
        EIy, ...      
        EIxy, ...      
        flapEI_cm, ...
        edgeEI_cm, ...
        flapEI_tc, ...
        edgeEI_tc, ...
        flapEI_sc, ...
        edgeEI_sc, ...
        mIx, ...
        mIy, ...      
        mIxy, ...        
        flapIner_cm, ... 
        edgeIner_cm, ...
        flapIner_tc, ...
        edgeIner_tc, ...
        flapIner_sc, ...
        edgeIner_sc, ...
        cm_offst, ...    
        tc_offst, ...   
        sc_offst, ...   
        x_cm, ...        
        y_cm, ...        
        x_tc, ...        
        y_tc, ...        
        x_sc, ...        
        y_sc];   
if any(any(isnan(data))) || any(any(isinf(data)))
    propError = true;
    fprintf(1,'WARNING: Co-Blade may have failed during the structural analysis. Continuing anyways...\n');
else
    propError = false;
end

%% Write a BModes input file, execute BModes, and read the output file    
if N_MODES >= 1 && ~propError
    writeInpFileBModes(SIM, ...
                       ANLS, ...
                       BLADE, ...
                       elas_tw, ...
                       iner_tw, ...
                       mass_den, ...
                       flapIner_cm, ...
                       edgeIner_cm, ...
                       flapEI_sc, ...
                       edgeEI_sc, ...
                       tor_stff, ...
                       axial_stff, ...
                       cm_offst, ...
                       sc_offst, ...
                       tc_offst);

    % execute BModes
    cd(SIM.outputDir);
    [status result] = system(['BModes ' SIM.case '_BModes.bmi']);
    cd(SIM.rootDir);

    % read BModes output file
    Modes = readOutFileBModes(SIM, ANLS);
else
    if N_MODES >= 1 && propError
        fprintf(1,'WARNING: skipping BModes analysis due to possible error in stuctural properties. Continuing anyways...\n');
    end
    Modes.natFreq    = Inf; % create placeholder values
    Modes.span_loc   = [];
    Modes.flap_disp  = [];
    Modes.flap_slope = [];
    Modes.lag_disp   = [];
    Modes.lag_slope  = [];
    Modes.twist      = [];
end                

end % function structAnalysis

