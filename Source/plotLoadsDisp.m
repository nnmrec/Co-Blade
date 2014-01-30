function plotLoadsDisp(SIM, BLADE, AF, Coord, AppLoads, Disp, StrProps, OUT)
                       
%% re-assign some structure variable names (for convenience)
zSec         = BLADE.zSec;
NUM_SEC      = BLADE.NUM_SEC;
INTERP_AF    = BLADE.INTERP_AF;
N_AF         = BLADE.N_AF;
pitAxis      = BLADE.pitAxis;
chord        = BLADE.chord;
aeroTwst     = BLADE.aeroTwst;
NormAFcoords = AF.NormAFcoords; 

%%
figTitle = [SIM.case ' Blade Loads and Displacement'];
fig = figure('name', figTitle, ...
             'color', 'white', ...
             'units','normalized',...
             'outerposition',[0.1 0.1 0.8 0.8]);

if OUT.PLOT_GBL_SYS    
    % translation vector from the blade origin to the global origin
    t_vec = Coord.O.B - Coord.O.G;
    
    % transformation matrix from the blade coord. sys. to the global coord. sys.
    T_mat = Coord.Tran.BG;
    
    % create arrows that indicate the global origin and orientation
    hold on
    plot3(Coord.O.G(1), Coord.O.G(2), Coord.O.G(3), '.k')
    plot3([Coord.O.G(1) Coord.O.G(1)+Coord.U.G(1)], [Coord.O.G(2) Coord.O.G(2)+Coord.U.G(2)], [Coord.O.G(3) Coord.O.G(3)+Coord.U.G(3)],'-r')
    plot3([Coord.O.G(1) Coord.O.G(1)+Coord.V.G(1)], [Coord.O.G(2) Coord.O.G(2)+Coord.V.G(2)], [Coord.O.G(3) Coord.O.G(3)+Coord.V.G(3)],'-b')
    plot3([Coord.O.G(1) Coord.O.G(1)+Coord.W.G(1)], [Coord.O.G(2) Coord.O.G(2)+Coord.W.G(2)], [Coord.O.G(3) Coord.O.G(3)+Coord.W.G(3)],'-g')
    
    az = -60;
    el = 30;  
else
    t_vec = [0; 0; 0];   
    T_mat = [1, 0, 0;
             0, 1, 0;    
             0, 0, 1];  
    az = 50;
    el = 20;    
end


switch INTERP_AF 
    case 'none'
        
        % we will have to interpolate the coordinates anyways to create the
        % surface plot.  The surf() function requires square matrices as
        % inputs, so we need to make all the coordinate matrices the same
        % size at each blade station.
        
        nAF = zeros(NUM_SEC, 1);
        for i = 1:NUM_SEC
            nAF(i) = numel(NormAFcoords(i).x);            
        end
        N_AF = max(nAF);
        
        for i = 1:NUM_SEC
            
            oldX           = NormAFcoords(i).x;
            oldY           = NormAFcoords(i).y;
            [unused i_TE]  = min(abs(oldX - 1)); % find the index of the trailing edge 
            oldXup         = oldX(1:i_TE);
            oldYup         = oldY(1:i_TE);
            oldXlo         = [oldX(i_TE:end); 0];
            oldYlo         = [oldY(i_TE:end); 0];
            
            x_cUp  = cosspace(0, 1, N_AF/2 + 1, 'both');
            x_cLo  = flipud( x_cUp );
            newYup = interp1(oldXup, oldYup, x_cUp);
            newYlo = interp1(oldXlo, oldYlo, x_cLo);
            NormAFcoords(i).x = [x_cUp; x_cLo(2:end-1)];
            NormAFcoords(i).y = [newYup; newYlo(2:end-1)];
        
        end
         
    otherwise
        
end

X  = zeros(NUM_SEC, N_AF+1);
Y  = zeros(NUM_SEC, N_AF+1);
Z  = zeros(NUM_SEC, N_AF+1);
XX = zeros(NUM_SEC, N_AF+1);
YY = zeros(NUM_SEC, N_AF+1);
ZZ = zeros(NUM_SEC, N_AF+1);
for i = 1:NUM_SEC    
    x_c    = NormAFcoords(i).x;  % normalized x-coordinate of airfoil nodes
    y_c    = NormAFcoords(i).y;  % normalized y-coordinate of airfoil nodes
    [x y]  = scaleAirfoilCoords([x_c; x_c(1)], [y_c; y_c(1)], pitAxis(i), chord(i), aeroTwst(i));   
    X(i,:) = x;
    Y(i,:) = y;
    Z(i,:) = zSec(i).*ones(1, N_AF+1);

    r = T_mat * [x'; y'; zSec(i).*ones(1, N_AF+1)] + t_vec * ones(1, N_AF+1);   % rotation and translation between coordinate systems

    XX(i,:) = r(1,:);
    YY(i,:) = r(2,:);
    ZZ(i,:) = r(3,:);    
end

%%
if OUT.PLOT_GBL_SYS
    % create plot w.r.t. the global coordinate system
    for i = 1:NUM_SEC
        hold on
        plot3(XX(i,:), YY(i,:), ZZ(i,:), '-k') % plot the airfoil perimeter
    end
    plot3(XX(:,1),        YY(:,1),        ZZ(:,1),        '-k')    % plot the leading edge
    plot3(XX(:,N_AF/2+1), YY(:,N_AF/2+1), ZZ(:,N_AF/2+1), '-k')    % plot the trailing edge
    s = surf(XX, YY, ZZ);
    set(s,'EdgeColor','none','FaceColor',[0.8 0.8 0.8],'FaceAlpha',0.4)
    view(az, el)
    xlabel('x (m)')
    ylabel('y (m)')
    zlabel('z (m)')
    
else
    % create plot w.r.t. the blade coordinate system
    for i = 1:NUM_SEC
        hold on
        plot3(ZZ(i,:), XX(i,:), YY(i,:), '-k') % plot the airfoil perimeter
    end
    plot3(ZZ(:,1),        XX(:,1),        YY(:,1),        '-k')    % plot the leading edge
    plot3(ZZ(:,N_AF/2+1), XX(:,N_AF/2+1), YY(:,N_AF/2+1), '-k')    % plot the trailing edge
    s = surf(ZZ, XX, YY);
    set(s,'EdgeColor','none','FaceColor',[0.8 0.8 0.8],'FaceAlpha',0.4)
    view(az, el)
    xlabel('z (m)')
    ylabel('x (m)')
    zlabel('y (m)')
    
end
     
%%
if OUT.PLOT_F_BLD
    % determine the scale factor for the length of the arrows
    maxLength   = 0.25 * zSec(end);
    maxForce    = max([ max(AppLoads.px_a), ...
                        max(AppLoads.py_a), ...
                        max(AppLoads.px_w), ...
                        max(AppLoads.py_w), ...
                        max(AppLoads.pz_w), ...
                        max(AppLoads.px_c), ...
                        max(AppLoads.py_c), ...
                        max(AppLoads.pz_c) ]);
    sf = maxLength / maxForce;

    xyz_pa = T_mat * [zeros(NUM_SEC,1)'; zeros(NUM_SEC,1)'; zSec'] + t_vec * ones(1, NUM_SEC);
    xx_pa   = xyz_pa(1,:);
    yy_pa   = xyz_pa(2,:);
    zz_pa   = xyz_pa(3,:);

    xyz_cm = T_mat * [StrProps.x_cm'; StrProps.y_cm'; zSec'] + t_vec * ones(1, NUM_SEC);
    xx_cm   = xyz_cm(1,:);
    yy_cm   = xyz_cm(2,:);
    zz_cm   = xyz_cm(3,:);

    uvw_a = T_mat * [AppLoads.px_a'; AppLoads.py_a'; zeros(NUM_SEC, 1)'] + t_vec * ones(1, NUM_SEC);
    u_a   = uvw_a(1,:) .* sf;
    v_a   = uvw_a(2,:) .* sf;
    w_a   = uvw_a(3,:) .* sf;

    uvw_w = T_mat * [AppLoads.px_w'; AppLoads.py_w'; AppLoads.pz_w'] + t_vec * ones(1, NUM_SEC);
    u_w   = uvw_w(1,:) .* sf;
    v_w   = uvw_w(2,:) .* sf;
    w_w   = uvw_w(3,:) .* sf;

    uvw_c = T_mat * [AppLoads.px_c'; AppLoads.py_c'; AppLoads.pz_c'] + t_vec * ones(1, NUM_SEC);
    u_c   = uvw_c(1,:) .* sf;
    v_c   = uvw_c(2,:) .* sf;
    w_c   = uvw_c(3,:) .* sf;
    
    if OUT.PLOT_GBL_SYS      
        quiver3(xx_pa, yy_pa, zz_pa, u_a, v_a, w_a, 0, 'b');   % aerodynamic forces
        quiver3(xx_cm, yy_cm, zz_cm, u_w, v_w, w_w, 0, 'g');   % self weight & buoynacy forces
        quiver3(xx_cm, yy_cm, zz_cm, u_c, v_c, w_c, 0, 'r');   % centrifugal forces          
    else
        quiver3(zz_pa, xx_pa, yy_pa, w_a, u_a, v_a, 0, 'b');   % aerodynamic forces
        quiver3(zz_cm, xx_cm, yy_cm, w_w, u_w, v_w, 0, 'g');   % self weight & buoynacy forces
        quiver3(zz_cm, xx_cm, yy_cm, w_c, u_c, v_c, 0, 'r');   % centrifugal forces
    end
    
    h = plot3(0,0,0,'b', ...
              0,0,0,'g', ...
              0,0,0,'r');
    legend(h, 'Aerodynamic Forces','Net Weight Forces','Centrifugal Forces', 'Location', 'best')
    
    axis image
end

%%
if OUT.PLOT_DISP_BLD
    % plot the displaced blade geometry
    
    t_sc = [StrProps.x_sc, StrProps.y_sc, zeros(NUM_SEC, 1)]';                    % translation of the shear center
    t_tc = [Disp.uo + StrProps.x_tc, Disp.vo + StrProps.y_tc, zSec + Disp.wo]';   % translation of the displaced tension center
    
    dX = zeros(NUM_SEC, N_AF+1);
    dY = zeros(NUM_SEC, N_AF+1);
    dZ = zeros(NUM_SEC, N_AF+1);
    for i = 1:NUM_SEC     
        
        % translate the displaced shear center of every cross section to
        % the origin of the coordinate system, and then rotate the cross
        % sections about the z-axis (accounts for the torsional twist)
        r1 = Disp.Rz{i} * [X(i,:) - StrProps.x_sc(i); Y(i,:) - StrProps.y_sc(i); 0.*Z(i,:)] + t_sc(:,i) * ones(1, N_AF+1);
                
        % next, rotate the cross section while its centroid is centered at the origin
        % (i.e. x = y = z = 0), and then translate the cross section back to the
        % correct (x,y,z) position and centroid offset and displacement offset
        r2 = Disp.Rxy{i} * [r1(1,:) - StrProps.x_tc(i); r1(2,:) - StrProps.y_tc(i); 0.*Z(i,:)] + t_tc(:,i) * ones(1, N_AF+1);
        
        r2 = T_mat * r2 + t_vec * ones(1, N_AF+1);
        
        dX(i,:) = r2(1,:);
        dY(i,:) = r2(2,:);
        dZ(i,:) = r2(3,:);
        
    end
         
    if OUT.PLOT_GBL_SYS
        for i = 1:NUM_SEC    
            hold on
            plot3(dX(i,:), dY(i,:), dZ(i,:), '-r') % plot the displaced airfoil perimeter
        end
        plot3(dX(:,1),        dY(:,1),        dZ(:,1),        '-r')    % plot the displaced leading edge
        plot3(dX(:,N_AF/2+1), dY(:,N_AF/2+1), dZ(:,N_AF/2+1), '-r')    % plot the displaced trailing edge 
        s = surf(dX, dY, dZ);
    else
        for i = 1:NUM_SEC    
            hold on
            plot3(dZ(i,:), dX(i,:), dY(i,:), '-r') % plot the displaced airfoil perimeter
        end
        plot3(dZ(:,1),        dX(:,1),        dY(:,1),        '-r')    % plot the displaced leading edge
        plot3(dZ(:,N_AF/2+1), dX(:,N_AF/2+1), dY(:,N_AF/2+1), '-r')    % plot the displaced trailing edge 
        s = surf(dZ, dX, dY);
    end
      
    set(s,'EdgeColor','none','FaceColor',[0.8 0.8 0.8],'FaceAlpha',0.4)
    axis image
end

if OUT.SAVE_PLOTS
    savePlots(fig, figTitle, OUT.SAVE_FIG_FMT)
end

end % function plotLoadsDisp
