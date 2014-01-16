function [xCapSt_inb, ...        
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
          t_web_core] = assignDesignVars(xo, OPT, BLADE, WEB, z_oub, z_CP)


%% re-assign some structure variable names (for convenience)
NUM_CP   = OPT.NUM_CP;
INB_STN  = OPT.INB_STN;
TRAN_STN = OPT.TRAN_STN;
OUB_STN  = OPT.OUB_STN;
NUM_SEC  = BLADE.NUM_SEC;
zSec     = BLADE.zSec;
pitAxis  = BLADE.pitAxis;
NUM_WEBS = WEB.NUM_WEBS;

%% assign the elemtents of vector xo to meaningful variable names
w_cap_inb   = xo(1);
w_cap_oub   = xo(2);
rootThick   = xo(3);
xo_t_pan    = reshape(xo(4:3+NUM_CP*5), NUM_CP, 5);
xo_t_web    = reshape(xo(4+NUM_CP*5:end), 2, 2);

% web chordwise locations
xCapSt_inb  = pitAxis(INB_STN) - w_cap_inb/2;
xCapEnd_inb = pitAxis(INB_STN) + w_cap_inb/2;
xCapSt_oub  = pitAxis(OUB_STN) - w_cap_oub/2;
xCapEnd_oub = pitAxis(OUB_STN) + w_cap_oub/2;
if NUM_WEBS > 1
    inbChLoc = linspace(xCapSt_inb, xCapEnd_inb, NUM_WEBS)';
    oubChLoc = linspace(xCapSt_oub, xCapEnd_oub, NUM_WEBS)';
else
    inbChLoc = (xCapSt_inb + xCapEnd_inb)/2;
    oubChLoc = (xCapSt_oub + xCapEnd_oub)/2;
end 

% root build-up thickness
m = -rootThick / (zSec(TRAN_STN) - zSec(INB_STN));    % slope
b = rootThick - m*zSec(INB_STN);                      % intercept
t_blade_root                   = zeros(NUM_SEC, 1);   % thickness of root build-up material
t_blade_root(1:INB_STN)        = rootThick;
t_blade_root(INB_STN:TRAN_STN) = m*zSec(INB_STN:TRAN_STN) + b;    

% blade skin thickness
t_blade_skin                   = zeros(NUM_SEC, 1);
t_blade_skin(1:TRAN_STN)       = xo_t_pan(1,1);
t_blade_skin(TRAN_STN:OUB_STN) = interp1(z_CP, xo_t_pan(:,1), z_oub);
t_blade_skin(OUB_STN:end)      = xo_t_pan(end,1);

% spar cap lamina thicknesses
t_cap_uni                   = zeros(NUM_SEC, 1);   % spar cap uni-directional thickness
t_cap_uni(INB_STN:TRAN_STN) = interp1(zSec([INB_STN, TRAN_STN]), [0; xo_t_pan(1,2)], zSec(INB_STN:TRAN_STN)); 
t_cap_uni(TRAN_STN:OUB_STN) = interp1(z_CP, xo_t_pan(:,2), z_oub);

t_cap_core                   = zeros(NUM_SEC, 1);   % spar cap core thickness
t_cap_core(INB_STN:TRAN_STN) = interp1(zSec([INB_STN, TRAN_STN]), [0; xo_t_pan(1,3)], zSec(INB_STN:TRAN_STN)); 
t_cap_core(TRAN_STN:OUB_STN) = interp1(z_CP, xo_t_pan(:,3), z_oub);

% leading edge panel lamina thicknesses
t_lep_core                   = zeros(NUM_SEC, 1);   % leading edge panel core thickness
t_lep_core(INB_STN:TRAN_STN) = interp1(zSec([INB_STN, TRAN_STN]), [0; xo_t_pan(1,4)], zSec(INB_STN:TRAN_STN)); 
t_lep_core(TRAN_STN:OUB_STN) = interp1(z_CP, xo_t_pan(:,4), z_oub);

% trailing edge panel lamina thicknesses
t_tep_core                   = zeros(NUM_SEC, 1);   % trailing edge panel core thickness
t_tep_core(INB_STN:TRAN_STN) = interp1(zSec([INB_STN, TRAN_STN]), [0; xo_t_pan(1,5)], zSec(INB_STN:TRAN_STN)); 
t_tep_core(TRAN_STN:OUB_STN) = interp1(z_CP, xo_t_pan(:,5), z_oub);

% web lamina thicknesses
t_web_skin = zeros(NUM_SEC, 1); % web skin thickness
t_web_skin(INB_STN:OUB_STN) = interp1(zSec([INB_STN, OUB_STN]), xo_t_web(:,1), zSec(INB_STN:OUB_STN));

t_web_core = zeros(NUM_SEC, 1); % web core thickness
t_web_core(INB_STN:OUB_STN) = interp1(zSec([INB_STN, OUB_STN]), xo_t_web(:,2), zSec(INB_STN:OUB_STN));

end % function assignDesignVars

