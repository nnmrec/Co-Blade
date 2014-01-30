function plotSurfaceData(SIM, OPT, BLADE, WEB, AF, stn1, stn2, Panel, PlotSurf, LaminaSS) 

%% re-assign some structure variable names (for convenience)
INB_STN          = OPT.INB_STN;
OUB_STN          = OPT.OUB_STN;
NUM_SEC          = BLADE.NUM_SEC;
zSec             = BLADE.zSec;
aeroTwst         = BLADE.aeroTwst;
chord            = BLADE.chord;
pitAxis          = BLADE.pitAxis;
NormAFcoords     = AF.NormAFcoords;
OrigNormAFcoords = AF.OrigNormAFcoords;
INTERP_AF        = BLADE.INTERP_AF;
N_AF             = BLADE.N_AF;
WEB_NODES        = WEB.WEB_NODES;                    
           
%%                     
if PlotSurf.LE_TE || PlotSurf.CrossSec

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

        r = [x'; y'; zSec(i).*ones(1, N_AF+1)];  

        XX(i,:) = r(1,:);
        YY(i,:) = r(2,:);
        ZZ(i,:) = r(3,:);    
    end

    % % create plot w.r.t. the blade coordinate system
    if PlotSurf.CrossSec
        for i = stn1:stn2
            hold on
            plot3(ZZ(i,:), XX(i,:), YY(i,:), '-k') % plot the airfoil perimeter
        end
    end
    if PlotSurf.LE_TE
        hold on
        plot3(ZZ(stn1:stn2,1),        XX(stn1:stn2,1),        YY(stn1:stn2,1),        '-k')    % plot the leading edge
        plot3(ZZ(stn1:stn2,N_AF/2+1), XX(stn1:stn2,N_AF/2+1), YY(stn1:stn2,N_AF/2+1), '-k')    % plot the trailing edge
    end
    
    if PlotSurf.Top || PlotSurf.Bot || PlotSurf.Web
        % don't plot this
        colorbar('delete')
    else
        hold on
        s = surf(ZZ, XX, YY);
    
        xlabel('z (m)')
        ylabel('x (m)')
        zlabel('y (m)')

        set(s,'EdgeColor','none','FaceColor',[0.8 0.8 0.8],'FaceAlpha',0.4)
    end
end

%% now plot the surface patches

% initialize a colormap, but it will get changed later
cmap   = colormap('jet'); % know thy enemy
cmtype = '';

% compute properties at midpoints between stations
zSec_midPt     = (    zSec(1:NUM_SEC-1) +     zSec(2:NUM_SEC))./2;
pitAxis_midPt  = ( pitAxis(1:NUM_SEC-1) +  pitAxis(2:NUM_SEC))./2;
chord_midPt    =    (chord(1:NUM_SEC-1) +    chord(2:NUM_SEC))./2;
aeroTwst_midPt = (aeroTwst(1:NUM_SEC-1) + aeroTwst(2:NUM_SEC))./2;

xWebNode = defineWebNodes(BLADE, WEB);
                                            
% need to get normalized airfoil coordinates (with identical x_c values) at each station
NormAFcoords = interpAirfoilCoords(OrigNormAFcoords, BLADE);

PanelSurf(NUM_SEC,1) = struct('Top',[],'Bot',[],'Web',[]);
for i = stn1:stn2
 
    if PlotSurf.Top % plot Top surfaces
        % airfoil coordinates at current station    
        [unused i_TE] = min(abs(NormAFcoords(i).x - 1));    % index of the trailing edge airfoil node
        x_c = NormAFcoords(i).x(1:i_TE);    % top surface only

        for n = eval(PlotSurf.TopPanel{i})
            % get x-y coordinates for exterior surface of panel
            xp     = Panel(i).Top.x{n};
            yp     = Panel(i).Top.y{n};
            xPanel = xp(1:numel(xp)/2); % only take the exterior surface of the panel coordinates
            yPanel = yp(1:numel(yp)/2);

            % unscale the coordinates
            [x_cPanel y_cPanel] = unscaleAirfoilCoords(xPanel, yPanel, pitAxis(i), chord(i), aeroTwst(i));

            if i == 1 % first station
                y_cOub = interp1(x_c, NormAFcoords(i+1).y(1:i_TE), x_cPanel);
                y_cMid = (y_cPanel + y_cOub)./2;

                % now scale the coordinates
                [x1 y1] = scaleAirfoilCoords(x_cPanel, y_cPanel, pitAxis(i), chord(i), aeroTwst(i));
                [x2 y2] = scaleAirfoilCoords(x_cPanel, y_cMid, pitAxis_midPt(i), chord_midPt(i), aeroTwst_midPt(i));
                z1      = zSec(i)       .* ones(numel(x1),1);
                z2      = zSec_midPt(i) .* ones(numel(x2),1);
            elseif i == NUM_SEC % last station
                y_cInb = interp1(x_c, NormAFcoords(i-1).y(1:i_TE), x_cPanel);
                y_cMid = (y_cPanel + y_cInb)./2;

                % now scale the coordinates
                [x1 y1] = scaleAirfoilCoords(x_cPanel, y_cPanel, pitAxis(i), chord(i), aeroTwst(i));
                [x2 y2] = scaleAirfoilCoords(x_cPanel, y_cMid, pitAxis_midPt(i-1), chord_midPt(i-1), aeroTwst_midPt(i-1));
                z1      = zSec(i)         .* ones(numel(x1),1);
                z2      = zSec_midPt(i-1) .* ones(numel(x2),1);
            else % inbetween first and last stations
                y_cInb = interp1(x_c, NormAFcoords(i-1).y(1:i_TE), x_cPanel);
                y_cOub = interp1(x_c, NormAFcoords(i+1).y(1:i_TE), x_cPanel);
                y_cMidInb = (y_cPanel + y_cInb)./2;
                y_cMidOub = (y_cPanel + y_cOub)./2;

                % now scale the coordinates
                [x1 y1] = scaleAirfoilCoords(x_cPanel, y_cMidInb, pitAxis_midPt(i-1), chord_midPt(i-1), aeroTwst_midPt(i-1));
                [x2 y2] = scaleAirfoilCoords(x_cPanel, y_cPanel, pitAxis(i), chord(i), aeroTwst(i));
                [x3 y3] = scaleAirfoilCoords(x_cPanel, y_cMidOub, pitAxis_midPt(i), chord_midPt(i), aeroTwst_midPt(i));
                z1      = zSec_midPt(i-1) .* ones(numel(x1),1);
                z2      = zSec(i)         .* ones(numel(x2),1);
                z3      = zSec_midPt(i)   .* ones(numel(x3),1);
            end

            switch PlotSurf.Data
                case 'Panel Data'
                    % get the property for the selected panel
                    [c cbTitle cmtype] = getPanelProperty(Panel(i).Top, n, PlotSurf.panelData);
                    c                  = ones(numel(x1),1) .* c;
                case 'Layer Data'
                    % get the property for the selected layer
                    iLayer             = eval(PlotSurf.TopLayer(i).layer{n});
                    [c cbTitle cmtype] = getLayerProperty(LaminaSS(i).Top, n, iLayer, PlotSurf.layerData);
                otherwise
                    error('ERROR: unrecognized choice for Panel or Lamina data');
            end
            
            if i == 1 || i == NUM_SEC
                % get the verticies for the patch polygons in 3D
                x1_vert = zeros(3,numel(x1)-1);
                y1_vert = zeros(3,numel(x1)-1);
                z1_vert = zeros(3,numel(x1)-1);
                c1_vert = zeros(3,numel(x1)-1);
                x2_vert = zeros(3,numel(x1)-1);
                y2_vert = zeros(3,numel(x1)-1);
                z2_vert = zeros(3,numel(x1)-1);
                c2_vert = zeros(3,numel(x1)-1);
                for k = 1:numel(x1)-1
                   x1_vert(:,k)   = [x1(k); x2(k);   x1(k+1)]; 
                   y1_vert(:,k)   = [y1(k); y2(k);   y1(k+1)]; 
                   z1_vert(:,k)   = [z1(k); z2(k);   z1(k+1)]; 
                   c1_vert(:,k)   = [ c(k);  c(k);    c(k+1)];
                   x2_vert(:,k+1) = [x2(k); x2(k+1); x1(k+1)];
                   y2_vert(:,k+1) = [y2(k); y2(k+1); y1(k+1)];
                   z2_vert(:,k+1) = [z2(k); z2(k+1); z1(k+1)];   
                   c2_vert(:,k+1) = [ c(k);  c(k+1);  c(k+1)];
                end
                PanelSurf(i).Top.x{n} = [x1_vert x2_vert];
                PanelSurf(i).Top.y{n} = [y1_vert y2_vert];
                PanelSurf(i).Top.z{n} = [z1_vert z2_vert];
                PanelSurf(i).Top.c{n} = [c1_vert c2_vert];
            else
                % get the verticies for the patch polygons in 3D
                x1_vert = zeros(3,numel(x1)-1);
                y1_vert = zeros(3,numel(x1)-1);
                z1_vert = zeros(3,numel(x1)-1);
                c1_vert = zeros(3,numel(x1)-1);
                x2_vert = zeros(3,numel(x1)-1);
                y2_vert = zeros(3,numel(x1)-1);
                z2_vert = zeros(3,numel(x1)-1);
                c2_vert = zeros(3,numel(x1)-1);
                x3_vert = zeros(3,numel(x1)-1);
                y3_vert = zeros(3,numel(x1)-1);
                z3_vert = zeros(3,numel(x1)-1);
                c3_vert = zeros(3,numel(x1)-1);
                x4_vert = zeros(3,numel(x1)-1);
                y4_vert = zeros(3,numel(x1)-1);
                z4_vert = zeros(3,numel(x1)-1);
                c4_vert = zeros(3,numel(x1)-1);
                for k = 1:numel(x1)-1
                   x1_vert(:,k)   = [x1(k); x2(k);   x1(k+1)]; 
                   y1_vert(:,k)   = [y1(k); y2(k);   y1(k+1)]; 
                   z1_vert(:,k)   = [z1(k); z2(k);   z1(k+1)];
                   c1_vert(:,k)   = [ c(k);  c(k);    c(k+1)];
                   x2_vert(:,k+1) = [x2(k); x2(k+1); x1(k+1)];
                   y2_vert(:,k+1) = [y2(k); y2(k+1); y1(k+1)];
                   z2_vert(:,k+1) = [z2(k); z2(k+1); z1(k+1)];
                   c2_vert(:,k+1) = [ c(k);  c(k+1);  c(k+1)];
                   x3_vert(:,k)   = [x2(k); x3(k);   x2(k+1)]; 
                   y3_vert(:,k)   = [y2(k); y3(k);   y2(k+1)]; 
                   z3_vert(:,k)   = [z2(k); z3(k);   z2(k+1)]; 
                   c3_vert(:,k)   = [ c(k);  c(k);    c(k+1)];
                   x4_vert(:,k+1) = [x3(k); x3(k+1); x2(k+1)];
                   y4_vert(:,k+1) = [y3(k); y3(k+1); y2(k+1)];
                   z4_vert(:,k+1) = [z3(k); z3(k+1); z2(k+1)];  
                   c4_vert(:,k+1) = [ c(k);  c(k+1);  c(k+1)];
                end
                PanelSurf(i).Top.x{n} = [x1_vert x2_vert x3_vert x4_vert];
                PanelSurf(i).Top.y{n} = [y1_vert y2_vert y3_vert y4_vert];
                PanelSurf(i).Top.z{n} = [z1_vert z2_vert z3_vert z4_vert];
                PanelSurf(i).Top.c{n} = [c1_vert c2_vert c3_vert c4_vert];              
                
            end

            x = PanelSurf(i).Top.x{n};
            y = PanelSurf(i).Top.y{n};
            z = PanelSurf(i).Top.z{n};
            c = PanelSurf(i).Top.c{n};
            patch(z, x, y, c, 'EdgeColor', 'none')
        end
    end
    
    if PlotSurf.Bot % plot Bottom surfaces
        % airfoil coordinates at current station    
        [unused i_TE] = min(abs(NormAFcoords(i).x - 1));    % index of the trailing edge airfoil node
        x_c = [NormAFcoords(i).x(i_TE:end); 0];    % bottom surface only

        for n = eval(PlotSurf.BotPanel{i})
            % get x-y coordinates for exterior surface of panel
            xp     = Panel(i).Bot.x{n};
            yp     = Panel(i).Bot.y{n};
            xPanel = xp(1:numel(xp)/2); % only take the exterior surface of the panel coordinates
            yPanel = yp(1:numel(yp)/2);

            % unscale the coordinates
            [x_cPanel y_cPanel] = unscaleAirfoilCoords(xPanel, yPanel, pitAxis(i), chord(i), aeroTwst(i));

            if i == 1 % first station
                y_cOub = interp1(x_c, [NormAFcoords(i+1).y(i_TE:end); 0], x_cPanel);
                y_cMid = (y_cPanel + y_cOub)./2;

                % now scale the coordinates
                [x1 y1] = scaleAirfoilCoords(x_cPanel, y_cPanel, pitAxis(i), chord(i), aeroTwst(i));
                [x2 y2] = scaleAirfoilCoords(x_cPanel, y_cMid, pitAxis_midPt(i), chord_midPt(i), aeroTwst_midPt(i));
                z1      = zSec(i)       .* ones(numel(x1),1);
                z2      = zSec_midPt(i) .* ones(numel(x2),1);
            elseif i == NUM_SEC % last station
                y_cInb = interp1(x_c, [NormAFcoords(i-1).y(i_TE:end); 0], x_cPanel);
                y_cMid = (y_cPanel + y_cInb)./2;

                % now scale the coordinates
                [x1 y1] = scaleAirfoilCoords(x_cPanel, y_cPanel, pitAxis(i), chord(i), aeroTwst(i));
                [x2 y2] = scaleAirfoilCoords(x_cPanel, y_cMid, pitAxis_midPt(i-1), chord_midPt(i-1), aeroTwst_midPt(i-1));
                z1      = zSec(i)         .* ones(numel(x1),1);
                z2      = zSec_midPt(i-1) .* ones(numel(x2),1);
            else % inbetween first and last stations
                y_cInb = interp1(x_c, [NormAFcoords(i-1).y(i_TE:end); 0], x_cPanel);
                y_cOub = interp1(x_c, [NormAFcoords(i+1).y(i_TE:end); 0], x_cPanel);
                y_cMidInb = (y_cPanel + y_cInb)./2;
                y_cMidOub = (y_cPanel + y_cOub)./2;

                % now scale the coordinates
                [x1 y1] = scaleAirfoilCoords(x_cPanel, y_cMidInb, pitAxis_midPt(i-1), chord_midPt(i-1), aeroTwst_midPt(i-1));
                [x2 y2] = scaleAirfoilCoords(x_cPanel, y_cPanel, pitAxis(i), chord(i), aeroTwst(i));
                [x3 y3] = scaleAirfoilCoords(x_cPanel, y_cMidOub, pitAxis_midPt(i), chord_midPt(i), aeroTwst_midPt(i));
                z1      = zSec_midPt(i-1) .* ones(numel(x1),1);
                z2      = zSec(i)         .* ones(numel(x2),1);
                z3      = zSec_midPt(i)   .* ones(numel(x3),1);
            end

            switch PlotSurf.Data
                case 'Panel Data'
                    % get the property for the selected panel
                    [c cbTitle cmtype] = getPanelProperty(Panel(i).Bot, n, PlotSurf.panelData);
                    c                  = ones(numel(x1),1) .* c;
                case 'Layer Data'
                    % get the property for the selected layer
                    iLayer             = eval(PlotSurf.BotLayer(i).layer{n});
                    [c cbTitle cmtype] = getLayerProperty(LaminaSS(i).Bot, n, iLayer, PlotSurf.layerData);
                otherwise
                    error('ERROR: unrecognized choice for Panel or Lamina data');
            end
            
            if i == 1 || i == NUM_SEC
                % get the verticies for the patch polygons in 3D
                x1_vert = zeros(3,numel(x1)-1);
                y1_vert = zeros(3,numel(x1)-1);
                z1_vert = zeros(3,numel(x1)-1);
                c1_vert = zeros(3,numel(x1)-1);
                x2_vert = zeros(3,numel(x1)-1);
                y2_vert = zeros(3,numel(x1)-1);
                z2_vert = zeros(3,numel(x1)-1);
                c2_vert = zeros(3,numel(x1)-1);
                for k = 1:numel(x1)-1
                   x1_vert(:,k)   = [x1(k); x2(k);   x1(k+1)]; 
                   y1_vert(:,k)   = [y1(k); y2(k);   y1(k+1)]; 
                   z1_vert(:,k)   = [z1(k); z2(k);   z1(k+1)]; 
                   c1_vert(:,k)   = [ c(k);  c(k);    c(k+1)];
                   x2_vert(:,k+1) = [x2(k); x2(k+1); x1(k+1)];
                   y2_vert(:,k+1) = [y2(k); y2(k+1); y1(k+1)];
                   z2_vert(:,k+1) = [z2(k); z2(k+1); z1(k+1)];   
                   c2_vert(:,k+1) = [ c(k);  c(k+1);  c(k+1)];
                end
                PanelSurf(i).Bot.x{n} = [x1_vert x2_vert];
                PanelSurf(i).Bot.y{n} = [y1_vert y2_vert];
                PanelSurf(i).Bot.z{n} = [z1_vert z2_vert];
                PanelSurf(i).Bot.c{n} = [c1_vert c2_vert];
            else
                % get the verticies for the patch polygons in 3D
                x1_vert = zeros(3,numel(x1)-1);
                y1_vert = zeros(3,numel(x1)-1);
                z1_vert = zeros(3,numel(x1)-1);
                c1_vert = zeros(3,numel(x1)-1);
                x2_vert = zeros(3,numel(x1)-1);
                y2_vert = zeros(3,numel(x1)-1);
                z2_vert = zeros(3,numel(x1)-1);
                c2_vert = zeros(3,numel(x1)-1);
                x3_vert = zeros(3,numel(x1)-1);
                y3_vert = zeros(3,numel(x1)-1);
                z3_vert = zeros(3,numel(x1)-1);
                c3_vert = zeros(3,numel(x1)-1);
                x4_vert = zeros(3,numel(x1)-1);
                y4_vert = zeros(3,numel(x1)-1);
                z4_vert = zeros(3,numel(x1)-1);
                c4_vert = zeros(3,numel(x1)-1);
                for k = 1:numel(x1)-1
                   x1_vert(:,k)   = [x1(k); x2(k);   x1(k+1)]; 
                   y1_vert(:,k)   = [y1(k); y2(k);   y1(k+1)]; 
                   z1_vert(:,k)   = [z1(k); z2(k);   z1(k+1)];
                   c1_vert(:,k)   = [ c(k);  c(k);    c(k+1)];
                   x2_vert(:,k+1) = [x2(k); x2(k+1); x1(k+1)];
                   y2_vert(:,k+1) = [y2(k); y2(k+1); y1(k+1)];
                   z2_vert(:,k+1) = [z2(k); z2(k+1); z1(k+1)];
                   c2_vert(:,k+1) = [ c(k);  c(k+1);  c(k+1)];
                   x3_vert(:,k)   = [x2(k); x3(k);   x2(k+1)]; 
                   y3_vert(:,k)   = [y2(k); y3(k);   y2(k+1)]; 
                   z3_vert(:,k)   = [z2(k); z3(k);   z2(k+1)]; 
                   c3_vert(:,k)   = [ c(k);  c(k);    c(k+1)];
                   x4_vert(:,k+1) = [x3(k); x3(k+1); x2(k+1)];
                   y4_vert(:,k+1) = [y3(k); y3(k+1); y2(k+1)];
                   z4_vert(:,k+1) = [z3(k); z3(k+1); z2(k+1)];  
                   c4_vert(:,k+1) = [ c(k);  c(k+1);  c(k+1)];
                end
                PanelSurf(i).Bot.x{n} = [x1_vert x2_vert x3_vert x4_vert];
                PanelSurf(i).Bot.y{n} = [y1_vert y2_vert y3_vert y4_vert];
                PanelSurf(i).Bot.z{n} = [z1_vert z2_vert z3_vert z4_vert];
                PanelSurf(i).Bot.c{n} = [c1_vert c2_vert c3_vert c4_vert];              
                
            end

            x = PanelSurf(i).Bot.x{n};
            y = PanelSurf(i).Bot.y{n};
            z = PanelSurf(i).Bot.z{n};
            c = PanelSurf(i).Bot.c{n};
            patch(z, x, y, c, 'EdgeColor', 'none')
        end
    end
end

if PlotSurf.Web % plot Web surfaces
    
    for i = stn1:stn2
      
        % get airfoil coordinates on top and bottom
        [unused i_TE] = min(abs(NormAFcoords(i).x - 1));    % index of the trailing edge airfoil node
        x_cTop =  NormAFcoords(i).x(1:i_TE);    
        y_cTop =  NormAFcoords(i).y(1:i_TE);    
        x_cBot = [NormAFcoords(i).x(i_TE:end); 0];    
        y_cBot = [NormAFcoords(i).y(i_TE:end); 0];

        for n = eval(PlotSurf.WebPanel{i})

            y_cwTop       = interp1(x_cTop, y_cTop, xWebNode(i,n));
            y_cwBot       = interp1(x_cBot, y_cBot, xWebNode(i,n));
            [xwTop ywTop] = scaleAirfoilCoords(xWebNode(i,n), y_cwTop, pitAxis(i), chord(i), aeroTwst(i));
            [xwBot ywBot] = scaleAirfoilCoords(xWebNode(i,n), y_cwBot, pitAxis(i), chord(i), aeroTwst(i));

            if i == INB_STN % first station
                y_cwTopOub = interp1(x_cTop,  NormAFcoords(i+1).y(1:i_TE),       xWebNode(i+1,n));
                y_cwBotOub = interp1(x_cBot, [NormAFcoords(i+1).y(i_TE:end); 0], xWebNode(i+1,n));
                
                [xwTopOub ywTopOub] = scaleAirfoilCoords(xWebNode(i+1,n), y_cwTopOub, pitAxis(i+1), chord(i+1), aeroTwst(i+1));
                [xwBotOub ywBotOub] = scaleAirfoilCoords(xWebNode(i+1,n), y_cwBotOub, pitAxis(i+1), chord(i+1), aeroTwst(i+1));
                
                x1 = linspace(xwTop,xwBot,WEB_NODES)';
                y1 = linspace(ywTop,ywBot,WEB_NODES)';
                z1 = zSec(i) .* ones(numel(x1), 1);
                x2 = linspace((xwTop + xwTopOub)./2, (xwBot + xwBotOub)./2, WEB_NODES)';
                y2 = linspace((ywTop + ywTopOub)./2, (ywBot + ywBotOub)./2, WEB_NODES)';
                z2 = zSec_midPt(i) .* ones(numel(x2), 1);
            elseif i == OUB_STN % last station
                y_cwTopInb = interp1(x_cTop,  NormAFcoords(i-1).y(1:i_TE),       xWebNode(i-1,n));
                y_cwBotInb = interp1(x_cBot, [NormAFcoords(i-1).y(i_TE:end); 0], xWebNode(i-1,n));

                [xwTopInb ywTopInb] = scaleAirfoilCoords(xWebNode(i-1,n), y_cwTopInb, pitAxis(i-1), chord(i-1), aeroTwst(i-1));
                [xwBotInb ywBotInb] = scaleAirfoilCoords(xWebNode(i-1,n), y_cwBotInb, pitAxis(i-1), chord(i-1), aeroTwst(i-1));
                
                x1 = linspace((xwTop + xwTopInb)./2, (xwBot + xwBotInb)./2, WEB_NODES)';
                y1 = linspace((ywTop + ywTopInb)./2, (ywBot + ywBotInb)./2, WEB_NODES)';
                z1 = zSec_midPt(i-1) .* ones(numel(x1), 1);
                x2 = linspace(xwTop,xwBot,WEB_NODES)';
                y2 = linspace(ywTop,ywBot,WEB_NODES)';
                z2 = zSec(i) .* ones(numel(x2), 1);
            else % inbetween first and last stations
                y_cwTopInb = interp1(x_cTop,  NormAFcoords(i-1).y(1:i_TE),       xWebNode(i-1,n));
                y_cwBotInb = interp1(x_cBot, [NormAFcoords(i-1).y(i_TE:end); 0], xWebNode(i-1,n));
                 
                y_cwTopOub = interp1(x_cTop,  NormAFcoords(i+1).y(1:i_TE),       xWebNode(i+1,n));
                y_cwBotOub = interp1(x_cBot, [NormAFcoords(i+1).y(i_TE:end); 0], xWebNode(i+1,n));
                
                [xwTopInb ywTopInb] = scaleAirfoilCoords(xWebNode(i-1,n), y_cwTopInb, pitAxis(i-1), chord(i-1), aeroTwst(i-1));
                [xwBotInb ywBotInb] = scaleAirfoilCoords(xWebNode(i-1,n), y_cwBotInb, pitAxis(i-1), chord(i-1), aeroTwst(i-1));
                [xwTopOub ywTopOub] = scaleAirfoilCoords(xWebNode(i+1,n), y_cwTopOub, pitAxis(i+1), chord(i+1), aeroTwst(i+1));
                [xwBotOub ywBotOub] = scaleAirfoilCoords(xWebNode(i+1,n), y_cwBotOub, pitAxis(i+1), chord(i+1), aeroTwst(i+1));

                x1 = linspace((xwTop + xwTopInb)./2, (xwBot + xwBotInb)./2, WEB_NODES)';
                y1 = linspace((ywTop + ywTopInb)./2, (ywBot + ywBotInb)./2, WEB_NODES)';
                z1 = zSec_midPt(i-1) .* ones(numel(x1), 1);
                x2 = linspace(xwTop,xwBot,WEB_NODES)';
                y2 = linspace(ywTop,ywBot,WEB_NODES)';
                z2 = zSec(i) .* ones(numel(x2), 1);
                x3 = linspace((xwTop + xwTopOub)./2, (xwBot + xwBotOub)./2, WEB_NODES)';
                y3 = linspace((ywTop + ywTopOub)./2, (ywBot + ywBotOub)./2, WEB_NODES)';
                z3 = zSec_midPt(i) .* ones(numel(x3), 1);
            end
            
            switch PlotSurf.Data
                case 'Panel Data'
                    % get the property for the selected panel
                    [c cbTitle cmtype] = getPanelProperty(Panel(i).Web, n, PlotSurf.panelData);
                    c                  = ones(numel(x1),1) .* c;
                case 'Layer Data'
                    % get the property for the selected layer
                    iLayer             = eval(PlotSurf.WebLayer(i).layer{n});
                    [c cbTitle cmtype] = getLayerProperty(LaminaSS(i).Web, n, iLayer, PlotSurf.layerData);
                otherwise
                    error('ERROR: unrecognized choice for Panel or Lamina data');
            end

            if i == INB_STN || i == OUB_STN
                % get the verticies for the patch polygons in 3D
                x1_vert = zeros(3,numel(x1)-1);
                y1_vert = zeros(3,numel(x1)-1);
                z1_vert = zeros(3,numel(x1)-1);
                c1_vert = zeros(3,numel(x1)-1);
                x2_vert = zeros(3,numel(x1)-1);
                y2_vert = zeros(3,numel(x1)-1);
                z2_vert = zeros(3,numel(x1)-1);
                c2_vert = zeros(3,numel(x1)-1);
                for k = 1:numel(x1)-1
                   x1_vert(:,k)   = [x1(k); x2(k);   x1(k+1)]; 
                   y1_vert(:,k)   = [y1(k); y2(k);   y1(k+1)]; 
                   z1_vert(:,k)   = [z1(k); z2(k);   z1(k+1)]; 
                   c1_vert(:,k)   = [ c(k);  c(k);    c(k+1)];
                   x2_vert(:,k+1) = [x2(k); x2(k+1); x1(k+1)];
                   y2_vert(:,k+1) = [y2(k); y2(k+1); y1(k+1)];
                   z2_vert(:,k+1) = [z2(k); z2(k+1); z1(k+1)];   
                   c2_vert(:,k+1) = [ c(k);  c(k+1);  c(k+1)];
                end
                PanelSurf(i).Web.x{n} = [x1_vert x2_vert];
                PanelSurf(i).Web.y{n} = [y1_vert y2_vert];
                PanelSurf(i).Web.z{n} = [z1_vert z2_vert];
                PanelSurf(i).Web.c{n} = [c1_vert c2_vert];
            else
                % get the verticies for the patch polygons in 3D
                x1_vert = zeros(3,numel(x1)-1);
                y1_vert = zeros(3,numel(x1)-1);
                z1_vert = zeros(3,numel(x1)-1);
                c1_vert = zeros(3,numel(x1)-1);
                x2_vert = zeros(3,numel(x1)-1);
                y2_vert = zeros(3,numel(x1)-1);
                z2_vert = zeros(3,numel(x1)-1);
                c2_vert = zeros(3,numel(x1)-1);
                x3_vert = zeros(3,numel(x1)-1);
                y3_vert = zeros(3,numel(x1)-1);
                z3_vert = zeros(3,numel(x1)-1);
                c3_vert = zeros(3,numel(x1)-1);
                x4_vert = zeros(3,numel(x1)-1);
                y4_vert = zeros(3,numel(x1)-1);
                z4_vert = zeros(3,numel(x1)-1);
                c4_vert = zeros(3,numel(x1)-1);
                for k = 1:numel(x1)-1
                   x1_vert(:,k)   = [x1(k); x2(k);   x1(k+1)]; 
                   y1_vert(:,k)   = [y1(k); y2(k);   y1(k+1)]; 
                   z1_vert(:,k)   = [z1(k); z2(k);   z1(k+1)];
                   c1_vert(:,k)   = [ c(k);  c(k);    c(k+1)];
                   x2_vert(:,k+1) = [x2(k); x2(k+1); x1(k+1)];
                   y2_vert(:,k+1) = [y2(k); y2(k+1); y1(k+1)];
                   z2_vert(:,k+1) = [z2(k); z2(k+1); z1(k+1)];
                   c2_vert(:,k+1) = [ c(k);  c(k+1);  c(k+1)];
                   x3_vert(:,k)   = [x2(k); x3(k);   x2(k+1)]; 
                   y3_vert(:,k)   = [y2(k); y3(k);   y2(k+1)]; 
                   z3_vert(:,k)   = [z2(k); z3(k);   z2(k+1)]; 
                   c3_vert(:,k)   = [ c(k);  c(k);    c(k+1)];
                   x4_vert(:,k+1) = [x3(k); x3(k+1); x2(k+1)];
                   y4_vert(:,k+1) = [y3(k); y3(k+1); y2(k+1)];
                   z4_vert(:,k+1) = [z3(k); z3(k+1); z2(k+1)];  
                   c4_vert(:,k+1) = [ c(k);  c(k+1);  c(k+1)];
                end
                PanelSurf(i).Web.x{n} = [x1_vert x2_vert x3_vert x4_vert];
                PanelSurf(i).Web.y{n} = [y1_vert y2_vert y3_vert y4_vert];
                PanelSurf(i).Web.z{n} = [z1_vert z2_vert z3_vert z4_vert];
                PanelSurf(i).Web.c{n} = [c1_vert c2_vert c3_vert c4_vert];
            end

            x = PanelSurf(i).Web.x{n};
            y = PanelSurf(i).Web.y{n};
            z = PanelSurf(i).Web.z{n};
            c = PanelSurf(i).Web.c{n};
            patch(z, x, y, c, 'EdgeColor', 'none')

        end
    end
end
  
if PlotSurf.Top || PlotSurf.Bot || PlotSurf.Web
    cb = colorbar('location','NorthOutside');
    set(get(cb,'xlabel'), 'String', cbTitle);
end
xlabel('z (m)')
ylabel('x (m)')
zlabel('y (m)')
axis image
alpha(1)

switch cmtype
    case 'seq_grayscale'
        cmap = colormap(flipud(colormap('bone')));
        colormap(cmap)
    case 'seq_color'
%         cmap = colormap(cbrewer('seq', 'YlOrRd', 9));
%         colormap(cmap)

        nLevels = 11;
        cmap    = colormap(bipolar(2*nLevels, 0.51));
        colormap(cmap)
        cabs = max(abs(get(gca, 'CLim')));
        set(gca, 'CLim',  [-cabs cabs]);
        set(cb,  'XTick', linspace(0, cabs, nLevels+1));
        set(cb,  'XLim',  [0 cabs]);

    case 'div_color'
        cmx = max(abs(get(gca, 'CLim')));
        set(gca, 'CLim', [-cmx cmx]);
        cmap = colormap(bipolar(21, 0.51));
        colormap(cmap)
    otherwise
        colormap(cmap) % just use the initiated map
end

end % function plotSurfaceData

function [c cbTitle cmtype] = getPanelProperty(Input, n, dataType)

switch dataType
    case 't:      thickness'
        c       = Input.t(n) .* 1000;
        cbTitle = 'Panel Thickness, t [mm]';
        cmtype  = 'seq_grayscale';        
    case 'E_eff:  effective Young''s modulus'
        c       = Input.E_eff(n) ./ 1e9;
        cbTitle = 'Panel Effective Young''s Modulus, E_{eff} [GPa]';
        cmtype  = 'seq_grayscale';
    case 'G_eff:  effective shear modulus'
        c       = Input.G_eff(n) ./ 1e9;
        cbTitle = 'Panel Effective Shear Modulus, G_{eff} [GPa]';
        cmtype  = 'seq_grayscale';
    case 'nu_eff: effective Poisson ratio'
        c       = Input.nu_eff(n);
        cbTitle = 'Panel Effective Poisson Ratio, nu_{eff}';
        cmtype  = 'seq_grayscale';
    otherwise
        error('ERROR: unrecognized choice for Panel Data.')
end

end % function getPanelProperty

function [c cbTitle cmtype] = getLayerProperty(Input, n, m, dataType)

switch dataType
    case 'e_11_maxabs: normal strain'
        c       = Input(n).e_11_maxabs{m} .* 1e6;
        cbTitle = 'normal strain, \epsilon_{11} [\mu-m/m]';
        cmtype  = 'div_color';
    case 'e_22_maxabs: transverse strain'
        c       = Input(n).e_22_maxabs{m} .* 1e6;
        cbTitle = 'transverse strain, \epsilon_{22} [\mu-m/m]';
        cmtype  = 'div_color';
    case 'e_12_maxabs: shear strain'
        c       = Input(n).e_12_maxabs{m} .* 1e6;
        cbTitle = 'shear strain, | \gamma_{12} | [\mu-rad]';
        cmtype  = 'seq_color';
    case 's_11_maxabs: normal stress'
        c       = Input(n).s_11_maxabs{m} ./ 1e6;
        cbTitle = 'normal stress, \sigma_{11} [MPa]';
        cmtype  = 'div_color';
    case 's_22_maxabs: transverse stress'
        c       = Input(n).s_22_maxabs{m} ./ 1e6;
        cbTitle = 'transverse stress, \sigma_{22} [MPa]';
        cmtype  = 'div_color';
    case 's_12_maxabs: shear stress'
        c       = Input(n).s_12_maxabs{m} ./ 1e6;
        cbTitle = 'shear stress, | \tau_{12} | [MPa]';
        cmtype  = 'seq_color';
    case 's_11_fc_T:   max stress failure criteria (normal tension)'
        c       = Input(n).s_11_fc_T{m};
        cbTitle = 'max stress failure criteria (normal tension), \sigma_{11} / \sigma_{11, fT}';
        cmtype  = 'seq_color';
    case 's_11_fc_C:   max stress failure criteria (normal compression)'
        c       = Input(n).s_11_fc_C{m};
        cbTitle = 'max stress failure criteria (normal compression), \sigma_{11} / \sigma_{11, fC}';
        cmtype  = 'seq_color';
    case 's_22_fc_T:   max stress failure criteria (transverse tension)'
        c       = Input(n).s_22_fc_T{m};
        cbTitle = 'max stress failure criteria (transverse tension), \sigma_{22} / \sigma_{22, fT}';
        cmtype  = 'seq_color';
    case 's_22_fc_C:   max stress failure criteria (transverse compression)'
        c       = Input(n).s_22_fc_C{m};
        cbTitle = 'max stress failure criteria (transverse compression), \sigma_{22} / \sigma_{22, fC}';
        cmtype  = 'seq_color';
    case 's_12_fc_S:   max stress failure criteria (shear)'
        c       = Input(n).s_12_fc_S{m};
        cbTitle = 'max stress failure criteria (shear), |\tau_{12}| / \tau_{y}';
        cmtype  = 'seq_color';
    otherwise
        error('ERROR: unrecognized choice for Layer Data.')
end

end % function getLayerProperty