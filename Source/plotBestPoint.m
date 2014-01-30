function plotBestPoint(fig, x_current, OPT, BLADE, WEB, AF, z_oub, z_CP)

%% re-assign some structure variable names (for convenience)
NUM_CP       = OPT.NUM_CP;
INB_STN      = OPT.INB_STN;
TRAN_STN     = OPT.TRAN_STN;
OUB_STN      = OPT.OUB_STN;
NUM_SEC      = BLADE.NUM_SEC;
zSec         = BLADE.zSec;
BLD_LENGTH   = BLADE.BLD_LENGTH;
aeroTwst     = BLADE.aeroTwst;
chord        = BLADE.chord;
pitAxis      = BLADE.pitAxis;
WEB_NODES    = WEB.WEB_NODES;
nWebs        = WEB.nWebs;
NormAFcoords = AF.NormAFcoords;

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
 t_web_core] = assignDesignVars(x_current, OPT, BLADE, WEB, z_oub, z_CP);

xo_t_pan = reshape(x_current(4:3+NUM_CP*5), NUM_CP, 5);
xo_t_web = reshape(x_current(4+NUM_CP*5:end), 2, 2);
    
%% plot the panel thickness
set(0,'CurrentFigure',fig);
clf

%% plot the blade geometry
subplot(2,5,[6 10])

xCapSt        = planformLine(zSec, chord, pitAxis, xCapSt_inb,  xCapSt_oub,  INB_STN, OUB_STN);
xCapEnd       = planformLine(zSec, chord, pitAxis, xCapEnd_inb, xCapEnd_oub, INB_STN, OUB_STN);
xsec_node_st  = xCapSt  ./ chord + pitAxis; 
xsec_node_end = xCapEnd ./ chord + pitAxis; 

WEB.inbChLoc  = inbChLoc;
WEB.oubChLoc  = oubChLoc;
xWebNode      = defineWebNodes(BLADE, WEB);

xNodeTop = cell(NUM_SEC, 1);
xNodeBot = cell(NUM_SEC, 1);
nPanels  = zeros(NUM_SEC, 1);
for i = 1:NUM_SEC
    % on the top and bottom, assign the number of segments and chordwise locations of the segment boundaries
    if i <= INB_STN || i > OUB_STN
        % only 1 sector
        xNodeTop{i} = [0; 1];
        xNodeBot{i} = [0; 1];
    else
        % 3 sectors
        xNodeTop{i} = [0; xsec_node_st(i); xsec_node_end(i); 1];
        xNodeBot{i} = [0; xsec_node_st(i); xsec_node_end(i); 1];
    end       
    nPanels(i) = numel(xNodeTop{i}) - 1;
end

Panel(NUM_SEC,1) = struct('Top',[],'Bot',[],'Web',[]);
webLocs          = cell(NUM_SEC,1);
xw_top           = cell(NUM_SEC, 1);                    % x-coordinates of the middle of the shear webs on the top surface
yw_top           = cell(NUM_SEC, 1);                    % y-coordinates of the middle of the shear webs on the top surface
xw_bot           = cell(NUM_SEC, 1);                    % x-coordinates of the middle of the shear webs on the bottom surface
yw_bot           = cell(NUM_SEC, 1);                    % y-coordinates of the middle of the shear webs on the bottom surface
for i = 1:NUM_SEC
    % web nodes
    wl         = xWebNode(i,:);
    webLocs{i} = wl(~isnan(wl))'; % an empty array if no webs exist at this station

    % airfoil node coordinates 
    x             = NormAFcoords(i).x;  % normalized x-coordinate of airfoil nodes
    y             = NormAFcoords(i).y;  % normalized y-coordinate of airfoil nodes
    [unused i_TE] = min(abs(x - 1));    % index of the trailing edge airfoil node

    % airfoil coordinates with spar cap chordwise locations embedded
    x_top = consolidator([x(1:i_TE); xNodeTop{i}], [] ,[], 1e-4);                   % a sorted array in ascending order
    y_top = interp1(x(1:i_TE), y(1:i_TE), x_top);
    x_bot = consolidator([flipud([x(i_TE:end); 0]); xNodeBot{i}], [] ,[], 1e-4);	% a sorted array in ascending order
    y_bot = interp1([x(i_TE:end); 0], [y(i_TE:end); 0], x_bot);

    i_sNodesTop = findClosest(x_top, xNodeTop{i});      % indices of the sector nodes on the airfoil top 
    i_wNodesTop = findClosest(x_top, webLocs{i});       % indices of the web nodes on the airfoil top 
    i_sNodesBot = findClosest(x_bot, xNodeBot{i});      % indices of the sector nodes on the airfoil bottom
    i_wNodesBot = findClosest(x_bot, webLocs{i});       % indices of the web nodes on the airfoil top 

    % Now scale, translate, and rotate the all coordinates
    [x_top y_top] = scaleAirfoilCoords(x_top, y_top, pitAxis(i), chord(i), aeroTwst(i));
    [x_bot y_bot] = scaleAirfoilCoords(x_bot, y_bot, pitAxis(i), chord(i), aeroTwst(i));

    % endpoints of the shear webs
    xw_top{i} = x_top(i_wNodesTop);
    yw_top{i} = y_top(i_wNodesTop);
    xw_bot{i} = x_bot(i_wNodesBot);
    yw_bot{i} = y_bot(i_wNodesBot);

    for n = 1:nPanels(i)
        % determine the total thickness of the panel
        if     i <= INB_STN
            % blade root build-up
            tPly  = sum([t_blade_skin(i);
                         t_blade_root(i) * 2;
                         t_blade_skin(i)]);       
        elseif i > INB_STN && i < TRAN_STN
            % transition from root build-up
            if     n == 1  % leading edge panel
                tPly  = sum([t_blade_skin(i);
                             t_blade_root(i);
                             t_lep_core(i) * 2;
                             t_blade_root(i);
                             t_blade_skin(i)]);
            elseif n == 2  % spar cap panel
                tPly  = sum([t_blade_skin(i);
                             t_blade_root(i);
                             t_cap_uni(i);
                             t_cap_core(i) * 2;
                             t_cap_uni(i);
                             t_blade_root(i);
                             t_blade_skin(i)]);
            elseif n == 3  % trailing edge panel
                tPly  = sum([t_blade_skin(i);
                             t_blade_root(i);
                             t_tep_core(i) * 2;
                             t_blade_root(i);
                             t_blade_skin(i)]);
            end
        elseif i >= TRAN_STN && i <= OUB_STN  
            % root material no longer exists
            if     n == 1  % leading edge panel
                tPly  = sum([t_blade_skin(i);
                             t_lep_core(i) * 2;
                             t_blade_skin(i)]);
            elseif n == 2  % spar cap panel
                tPly  = sum([t_blade_skin(i);
                             t_cap_uni(i);
                             t_cap_core(i) * 2;
                             t_cap_uni(i);
                             t_blade_skin(i)]);
            elseif n == 3  % trailing edge panel
                tPly  = sum([t_blade_skin(i);
                             t_tep_core(i) * 2;
                             t_blade_skin(i)]);
            end
        elseif i > OUB_STN
            % blade tip
            tPly  = t_blade_skin(i) * 2;
        end % if i <= INB_STN 

        % airfoil nodes of the panel
        xTop = x_top(i_sNodesTop(n):i_sNodesTop(n+1));
        yTop = y_top(i_sNodesTop(n):i_sNodesTop(n+1)); 
        xBot = x_bot(i_sNodesBot(n):i_sNodesBot(n+1));
        yBot = y_bot(i_sNodesBot(n):i_sNodesBot(n+1));

        % x-y coordinates around the perimeter of the panel, forming a non-closed polygon 
        ang_top   = atan2( finiteDiff(yTop), finiteDiff(xTop) );
        ang_bot   = atan2( finiteDiff(yBot), finiteDiff(xBot) );
        xPoly_top = [xTop; flipud( xTop + tPly.*sin(ang_top) )];
        yPoly_top = [yTop; flipud( yTop - tPly.*cos(ang_top) )];
        xPoly_bot = [xBot; flipud( xBot - tPly.*sin(ang_bot) )];
        yPoly_bot = [yBot; flipud( yBot + tPly.*cos(ang_bot) )];

        Panel(i).Top.x{n}    = xPoly_top;
        Panel(i).Top.y{n}    = yPoly_top;
        Panel(i).Top.z{n}    = zSec(i) .* ones(size(xPoly_top));
        Panel(i).Top.nPanels = nPanels(i);
        Panel(i).Bot.x{n}    = xPoly_bot;
        Panel(i).Bot.y{n}    = yPoly_bot;
        Panel(i).Bot.z{n}    = zSec(i) .* ones(size(xPoly_bot));
        Panel(i).Bot.nPanels = nPanels(i);           
    end

    [x_web y_web] = webCoords(Panel(i), ...     % computes the x-y coordinates of the web mid-lines
                              xw_top{i}, ...
                              yw_top{i}, ...
                              xw_bot{i}, ...
                              yw_bot{i}, ...
                              nWebs(i), ...
                              WEB_NODES);

    for n = 1:nWebs(i)
        % x-y coordinates around the perimeter of the web panel, forming a non-closed polygon
        tWeb  = 2*(t_web_skin(i) + t_web_core(i));
        alfa  = atan2( finiteDiff(y_web{n}), finiteDiff(x_web{n}) );
        beta  = alfa  + pi/2;
        x1    = x_web{n} + (tWeb/2).*cos(beta);
        y1    = y_web{n} + (tWeb/2).*sin(beta);
        beta  = alfa  - pi/2;
        x2    = x_web{n} + (tWeb/2).*cos(beta);
        y2    = y_web{n} + (tWeb/2).*sin(beta);
        Panel(i).Web.x{n} = [x1; flipud(x2)];
        Panel(i).Web.y{n} = [y1; flipud(y2)];
        Panel(i).Web.z{n} = zSec(i) .* ones(size(Panel(i).Web.x{n}));
    end       

end

for i = 1:NUM_SEC
    % create a patch for the polygons of the top surface
    for n = 1:nPanels(i)
        x = Panel(i).Top.x{n};
        y = Panel(i).Top.y{n};
        z = Panel(i).Top.z{n};
        c = getColor(i,n,INB_STN,TRAN_STN,OUB_STN);
        patch(z, x, y, c, 'LineWidth', 0.01, 'FaceAlpha', 1)
    end

    % create a patch for the polygons of the bottom surface
    for n = 1:nPanels(i)
        x = Panel(i).Bot.x{n};
        y = Panel(i).Bot.y{n};
        z = Panel(i).Bot.z{n};
        c = getColor(i,n,INB_STN,TRAN_STN,OUB_STN);
        patch(z, x, y, c, 'LineWidth', 0.01, 'FaceAlpha', 1)
    end

    % create a patch for the polygons of the webs
    for n = 1:nWebs(i)
        x = Panel(i).Web.x{n};
        y = Panel(i).Web.y{n};
        z = Panel(i).Web.z{n};
        c = 0;
        patch(z, x, y, c, 'LineWidth', 0.01, 'FaceAlpha', 1)
    end
end

xlabel('z [m]')
ylabel('x [m]')
zlabel('y [m]')
view(40,8)
axis image
% alpha(0.75);

%%
% spar cap
subplot(2,5,1)
hold on
plot(zSec, 2000 .* t_blade_root, ':r')   	% root material
plot(zSec, 2000 .* t_blade_skin, '-b')   	% skin material
plot(zSec, 2000 .* t_cap_uni,    '.-k')     % uni-directional material
plot(zSec, 2000 .* t_cap_core,   'd-m')     % cap core material
plot(zSec(INB_STN), 2000 .* x_current(3), 'ok', 'MarkerFaceColor', 'g');
plot(z_CP, 2000 .* xo_t_pan(:,1), 'ok', 'MarkerFaceColor', 'g');
plot(z_CP, 2000 .* xo_t_pan(:,2), 'ok', 'MarkerFaceColor', 'g');
plot(z_CP, 2000 .* xo_t_pan(:,3), 'ok', 'MarkerFaceColor', 'g');    
title('Spar Caps')
xlabel('z [m]')
ylabel('t [mm]')
xlim([0 BLD_LENGTH])
box on

% leading edge panel
subplot(2,5,2)
hold on
plot(zSec, 2000 .* t_blade_root, ':r')      % root material
plot(zSec, 2000 .* t_blade_skin, '-b')      % skin material
plot(zSec, 2000 .* t_lep_core,   'o-m')     % lep core material
plot(zSec(INB_STN), 2000 .* x_current(3), 'ok', 'MarkerFaceColor', 'g');
plot(z_CP, 2000 .* xo_t_pan(:,1), 'ok', 'MarkerFaceColor', 'g');
plot(z_CP, 2000 .* xo_t_pan(:,4), 'ok', 'MarkerFaceColor', 'g');
title('Leading Edge Panels')
xlabel('z [m]')
xlim([0 BLD_LENGTH])
box on

% trailing edge panel
subplot(2,5,3)
hold on
plot(zSec, 2000 .* t_blade_root, ':r')   	% root material
plot(zSec, 2000 .* t_blade_skin, '-b')   	% skin material
plot(zSec, 2000 .* t_tep_core,   's-m')     % tep core material
plot(zSec(INB_STN), 2000 .* x_current(3), 'ok', 'MarkerFaceColor', 'g');
plot(z_CP, 2000 .* xo_t_pan(:,1), 'ok', 'MarkerFaceColor', 'g');
plot(z_CP, 2000 .* xo_t_pan(:,5), 'ok', 'MarkerFaceColor', 'g');
title('Trailing Edge Panels')
xlabel('z [m]')
xlim([0 BLD_LENGTH])
box on

% webs 
subplot(2,5,4)
hold on 
plot(Inf, Inf, ':r')    % plot fake data point, so we can create the legend entry blade-root
plot(Inf, Inf, '-b')    % plot fake data point, so we can create the legend entry blade-shell
plot(Inf, Inf, '.-k')   % plot fake data point, so we can create the legend entry spar-uni
plot(Inf, Inf, 'd-m')	% plot fake data point, so we can create the legend entry spar-core
plot(Inf, Inf, 'o-m') 	% plot fake data point, so we can create the legend entry LEP-core
plot(Inf, Inf, 's-m') 	% plot fake data point, so we can create the legend entry TEP-core
plot(zSec, 2000 .* t_web_skin, '+-b')      	 % web shell material
plot(zSec, 2000 .* t_web_core, 'x-m')        % web core material
plot(zSec([INB_STN OUB_STN]), 2000 .* xo_t_web(:,1), 'ok', 'MarkerFaceColor', 'g');
legend('blade-root', 'blade-shell', 'spar-uni', 'spar-core', ...
       'LEP-core', 'TEP-core', 'web-shell', 'web-core', 'control points', ...
       'Location','EastOutside')
plot(zSec([INB_STN OUB_STN]), 2000 .* xo_t_web(:,2), 'ok', 'MarkerFaceColor', 'g');
title('Shear Webs')
xlabel('z [m]')
xlim([0 BLD_LENGTH])
box on

end % function plotBestPoint

function c = getColor(i,n,INB_STN,TRAN_STN,OUB_STN)

    color = linspace(0,1,9);

    % determine the total thickness of the panel
    if     i <= INB_STN
        % blade root build-up
        c = color(9);       
    elseif i > INB_STN && i < TRAN_STN
        % transition from root build-up
        if     n == 1  % leading edge panel
            c = color(6);
        elseif n == 2  % spar cap panel
            c = color(8);
        elseif n == 3  % trailing edge panel
            c = color(4);
        end
    elseif i >= TRAN_STN && i <= OUB_STN  
        % root material no longer exists
        if     n == 1  % leading edge panel
            c = color(5);
        elseif n == 2  % spar cap panel
            c = color(7);
        elseif n == 3  % trailing edge panel
            c = color(3);
        end
    elseif i > OUB_STN
        % blade tip
        c = color(1);
    end % if i <= INB_STN
    
end % function getColor

