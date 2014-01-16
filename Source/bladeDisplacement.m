function Disp = bladeDisplacement(zSec, ...
                                  zFrac, ...
                                  NUM_SEC, ...
                                  ResLoads, ...
                                  flapI_mc, ...
                                  edgeI_mc, ...
                                  EIx, ...
                                  EIy, ...
                                  EIxy, ...
                                  axial_stff, ...
                                  tor_stff, ...
                                  Tran, ...
                                  DISP_CF)

Vz     = ResLoads.Vz;
Mx     = ResLoads.Mx;
My     = ResLoads.My;
dMz_dz = ResLoads.dMz_dz;
                              
%% Compute correction factors for the displacements (correction factors for displacement of tapered cantilever beams)     
if DISP_CF                              
    cf = dispCorrectionFactor(zFrac, NUM_SEC, flapI_mc, edgeI_mc);
else
    cf = ones(NUM_SEC, 1);
end

%% Compute the displacement of the centroidal (tension center) axis by integration of the 2nd order ODE's

d2uo_dz2 = cf .* ( My.*EIx + Mx.*EIxy)./(EIx.*EIy - EIxy.^2);
d2vo_dz2 = cf .* (-Mx.*EIy - My.*EIxy)./(EIx.*EIy - EIxy.^2);
dwo_dz   = Vz ./ axial_stff;

% numerically compute the indefinite integrals
duo_dz = cumtrapzf(zSec, d2uo_dz2);
dvo_dz = cumtrapzf(zSec, d2vo_dz2);
dtz_dz = cumtrapzf(zSec, dMz_dz);
% add the constant of integration, from the appropriate boundary condition
duo_dz = duo_dz - duo_dz(1);
dvo_dz = dvo_dz - dvo_dz(1);
dtz_dz = dtz_dz - dtz_dz(end);

% numerically compute the indefinite integrals
uo = cumtrapzf(zSec, duo_dz);
vo = cumtrapzf(zSec, dvo_dz);
wo = cumtrapzf(zSec, dwo_dz);
tz = cumtrapzf(zSec, dtz_dz ./ tor_stff);
% add the constant of integration, from the appropriate boundary condition
uo = uo - uo(1);
vo = vo - vo(1);
wo = wo - wo(1);
tz = tz - tz(1);

%% Compute the rigid body rotations of the cross sections
N    = NUM_SEC;
angY = zeros(N,1);
angX = zeros(N,1);
for i = 2:N
    dx      = uo(i) - uo(i-1);
    dy      = vo(i) - vo(i-1);
    dz      = (zSec(i) + wo(i)) - (zSec(i-1) + wo(i-1));
    angY(i) = atan2(dx, dz);
    angX(i) = atan2(dy, sqrt(dx^2 + dz^2)); 
end
Rz  = cell(N,1);
Rxy = cell(N,1);
for i = 1:N
    % rotation matrices
    RotZ = [  cos(tz(i)),    -sin(tz(i)),              0;
              sin(tz(i)),     cos(tz(i)),              0;
                       0,              0,              1]; 
    RotY = [ cos(angY(i)),             0,   sin(angY(i));
                        0,             1,              0;
            -sin(angY(i)),             0,   cos(angY(i))];
    RotX = [            1,             0,              0;
                        0, cos(-angX(i)), -sin(-angX(i));
                        0, sin(-angX(i)),  cos(-angX(i))];
    % this rotation matrix rotates the cross section so it is orthogonal 
    % to the centriodal (tension center) axis in accordance with the Euler-Bernoulli assumption
    Rz{i}  = RotZ;
    Rxy{i} = RotY * RotX;
end

%% Compute the component of displacment in the direction of the tower
disp_G     = Tran.BG * [uo'; vo'; wo']; % displacement vector w.r.t. the global coordinate system
tipDeflect = disp_G(1,:)';              % in the global coordinate system the x-component is in the direction of the tower

%% Collect the output
Disp.uo         = uo;
Disp.vo         = vo; 
Disp.wo         = wo;
Disp.tz         = tz;
Disp.Rz         = Rz;
Disp.Rxy        = Rxy;
Disp.tipDeflect = tipDeflect(end);  % only take the displacment at the tip of the blade, in the direction of the tower

end % function bladeDisplacement



