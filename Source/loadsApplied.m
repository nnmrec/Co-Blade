function AppLoads = loadsApplied(zSec, ...
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
                                 CENTRIF)

% vector of gravitational acceleration w.r.t. the global coordinate system
grav = [0; 0; -GRAV];      

% vectors of self weight force per unit length (span) w.r.t. the global coordinate system
if SELF_WEIGHT
    p_sw = grav * mass_den';
else
    p_sw = zeros(3, NUM_SEC);
end

% vectors of buoyancy force per unit length (span) w.r.t. the global coordinate system
if BUOYANCY
    p_b = -grav * (FLUID_DEN.*A_disp');
else
    p_b = zeros(3, NUM_SEC);
end

% vectors of net self weight force per unit length (span) w.r.t. the global coordinate system
p_w  = p_sw + p_b;

% transform the net buoyancy force vectors from the global coordinate system to the blade coordinate system
p_w  = Tran.GB*p_w; % net buoyancy force vectors w.r.t. the blade coordinate system
px_w = p_w(1,:)'; 
py_w = p_w(2,:)';
pz_w = p_w(3,:)';

% centrifugal forces per unit length w.r.t. the shaft coordinate system
if CENTRIF
    rSwept = (HUB_RAD + zSec').*cosd(PRE_CONE); % swept radius of a coned rotor
    p_c    = [zeros(1,NUM_SEC);
              zeros(1,NUM_SEC);
              rSwept .* mass_den' .* (ROT_SPD*pi/30)^2];
else
    p_c = zeros(3, NUM_SEC);
end

% transform the centrifugal force vectors from the shaft to the blade coordinate system    
p_c  = Tran.SB*p_c;
px_c = p_c(1,:)'; 
py_c = p_c(2,:)';
pz_c = p_c(3,:)';

%% Collect the output
AppLoads.px_a = px_a;   % x-component of aerodynamic force per unit length (span) w.r.t. the blade coordinate system
AppLoads.py_a = py_a;   % y-component of aerodynamic force per unit length (span) w.r.t. the blade coordinate system
AppLoads.qz_a = qz_a;   % z-component of aerodynamic moment per unit length (span) w.r.t. the blade coordinate system
AppLoads.px_w = px_w;   % x-component of net self weight force per unit length (span) w.r.t. the blade coordinate system
AppLoads.py_w = py_w;   % y-component of net self weight force per unit length (span) w.r.t. the blade coordinate system
AppLoads.pz_w = pz_w;   % z-component of net self weight force per unit length (span) w.r.t. the blade coordinate system
AppLoads.px_c = px_c;   % x-component of centrifugal force per unit length (span) w.r.t. the blade coordinate system
AppLoads.py_c = py_c;   % y-component of centrifugal force per unit length (span) w.r.t. the blade coordinate system
AppLoads.pz_c = pz_c;   % z-component of centrifugal force per unit length (span) w.r.t. the blade coordinate system

end % function loadsApplied

