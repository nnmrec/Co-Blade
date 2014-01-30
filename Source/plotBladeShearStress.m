function plotBladeShearStress(SIM, BLADE, WEB, OUT, Panel, ShearS)

figTitle = [SIM.case ' beam shear stress'];
fig = figure('name', figTitle, ...
             'color', 'white', ...
             'units','normalized',...
             'outerposition',[0.1 0.1 0.8 0.8]);

convStress = 1 / 1e6;  	% conversion factor for stresses

for i = 1:BLADE.NUM_SEC
    % create a patch for the polygons of the top surface
    for n = 1:Panel(i).Top.nPanels
        x = Panel(i).Top.x{n};
        y = Panel(i).Top.y{n};
        z = Panel(i).Top.z{n};
        d = ShearS(i).Top.stress_zs{n} .* convStress;
        c = [d; flipud(d)];
        patch(z, x, y, abs(c), 'LineWidth', 0.01) 
    end
    
    % create a patch for the polygons of the bottom surface
    for n = 1:Panel(i).Bot.nPanels
        x = Panel(i).Bot.x{n};
        y = Panel(i).Bot.y{n};
        z = Panel(i).Bot.z{n};
        d = ShearS(i).Bot.stress_zs{n} .* convStress;
        c = [d; flipud(d)];
        patch(z, x, y, abs(c), 'LineWidth', 0.01) 
    end
    
    % create a patch for the polygons of the shear webs
    for n = 1:WEB.nWebs(i)
        x = Panel(i).Web.x{n};
        y = Panel(i).Web.y{n};
        z = Panel(i).Web.z{n};
        d = ShearS(i).Web.stress_zs{n} .* convStress;
        c = [d; flipud(d)];
        patch(z, x, y, abs(c), 'LineWidth', 0.01)   
    end
    
end

xlabel('z (m)')
ylabel('x (m)')
zlabel('y (m)')
cb = colorbar('location','NorthOutside');
set(get(cb,'xlabel'), 'String', 'shear stress, | \tau_{zs} | (MPa)');
view(50,20)
axis image
alpha(1)
% levels    = 100;
% rgb_blue  = [59 76 192];
% rgb_red   = [180 4 38];
% colormap(diverging_map(linspace(0.5,1,levels), rgb_blue, rgb_red));
% cmx = max(abs(get(gca, 'CLim')));
% set(gca, 'CLim', [0 cmx]);
% cmap = colormap(cbrewer('seq', 'YlOrRd', 9));
% colormap(cmap)

nLevels = 11;
cmap    = colormap(bipolar(2*nLevels, 0.51));
colormap(cmap)
cabs = max(abs(get(gca, 'CLim')));
set(gca, 'CLim',  [-cabs cabs]);
set(cb,  'XTick', linspace(0, cabs, nLevels+1));
set(cb,  'XLim',  [0 cabs]);

set(fig, 'renderer', 'painters')
% set(fig, 'renderer', 'opengl')
% set(fig, 'renderer', 'zbuffer')

if OUT.SAVE_PLOTS
    savePlots(fig, figTitle, OUT.SAVE_FIG_FMT)
end

end % function plotBladeShearStress


