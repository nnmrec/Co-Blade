function plotBladeYModulus(SIM, BLADE, WEB, OUT, Panel)

figTitle = [SIM.case ' effective Young''s modulus'];
fig = figure('name', figTitle, ...
             'color', 'white', ...
             'units','normalized',...
             'outerposition',[0.1 0.1 0.8 0.8]);

convModulus = 1 / 1e9;  	% conversion factor for elastic modulus

for i = 1:BLADE.NUM_SEC
    % create a patch for the polygons of the top surface
    for n = 1:Panel(i).Top.nPanels
        x = Panel(i).Top.x{n};
        y = Panel(i).Top.y{n};
        z = Panel(i).Top.z{n};
        c = Panel(i).Top.E_eff(n) * convModulus;
        patch(z, x, y, c, 'LineWidth', 0.01)
    end
    
    % create a patch for the polygons of the bottom surface
    for n = 1:Panel(i).Bot.nPanels
        x = Panel(i).Bot.x{n};
        y = Panel(i).Bot.y{n};
        z = Panel(i).Bot.z{n};
        c = Panel(i).Bot.E_eff(n) * convModulus;
        patch(z, x, y, c, 'LineWidth', 0.01)
    end
    
    % create a patch for the polygons of the shear webs
    for n = 1:WEB.nWebs(i)
        x = Panel(i).Web.x{n};
        y = Panel(i).Web.y{n};
        z = Panel(i).Web.z{n};
        c = Panel(i).Web.E_eff(n) * convModulus;
        patch(z, x, y, c, 'LineWidth', 0.01)    
    end
    
end

xlabel('z (m)')
ylabel('x (m)')
zlabel('y (m)')
cb = colorbar('location','NorthOutside');
set(get(cb,'xlabel'), 'String', 'effective Young''s modulus, E_{eff} (GPa)');
view(50,20)
axis image
alpha(1)
colormap(flipud(colormap('gray')));


set(fig, 'renderer', 'opengl')

if OUT.SAVE_PLOTS
    savePlots(fig, figTitle, OUT.SAVE_FIG_FMT)
end

end % function plotBladeYModulus

