function plotLoads(SIM, BLADE, AppLoads, ResLoads, OUT)

convForce  = 1/1000; % convert N to kN 
convMoment = 1/1000; % convert N-m to kN-m 

%%
if OUT.PLOT_APPLOADS
    figTitle = [SIM.case ' applied loads'];
    fig = figure('name', figTitle, ...
                 'color', 'white', ...
                 'units','normalized',...
                 'outerposition',[0.1 0.1 0.8 0.8]);
    
    subplot(4,1,1);
    hold on
    plot(BLADE.zSec, AppLoads.px_a .* convForce, 'o-b')
    plot(BLADE.zSec, AppLoads.px_w .* convForce, 's-g')
    plot(BLADE.zSec, AppLoads.px_c .* convForce, '.-r')
    ylabel('p_x (kN/m)')
    box on
    
    legend('aerodynamic loads', ...
           'net weight loads', ...
           'centrifugal loads', ...
           'Location', 'EastOutside');
    pos = get(gca, 'Position');
    
    subplot(4,1,2);
    hold on
    plot(BLADE.zSec, AppLoads.py_a .* convForce, 'o-b')
    plot(BLADE.zSec, AppLoads.py_w .* convForce, 's-g')
    plot(BLADE.zSec, AppLoads.py_c .* convForce, '.-r')
    ylabel('p_y (kN/m)')
    box on
    p = get(gca, 'Position');
    set(gca, 'Position', [p(1) p(2) pos(3) p(4)]);    
    
    subplot(4,1,3)
    hold on
    plot(BLADE.zSec, AppLoads.pz_w .* convForce, 's-g')
    plot(BLADE.zSec, AppLoads.pz_c .* convForce, '.-r')
    ylabel('p_z (kN/m)')
    box on
    p = get(gca, 'Position');
    set(gca, 'Position', [p(1) p(2) pos(3) p(4)]);  
    
    subplot(4,1,4)
    hold on
    plot(BLADE.zSec, AppLoads.qz_a .* convMoment, 'o-b')
    xlabel('z (m)')
    ylabel('q_z (kN-m/m)')
    box on   
    p = get(gca, 'Position');
    set(gca, 'Position', [p(1) p(2) pos(3) p(4)]);  
    
    if OUT.SAVE_PLOTS
        savePlots(fig, figTitle, OUT.SAVE_FIG_FMT)
    end
end

%%
if OUT.PLOT_RESLOADS
    figTitle = [SIM.case ' resultant loads'];
    fig = figure('name', figTitle, ...
                 'color', 'white', ...
                 'units','normalized',...
                 'outerposition',[0.1 0.1 0.8 0.8]);
    
    subplot(2,1,1);
    hold on
    plot(BLADE.zSec, ResLoads.Vx .* convForce, 'o-b')
    plot(BLADE.zSec, ResLoads.Vy .* convForce, 's-g')
    plot(BLADE.zSec, ResLoads.Vz .* convForce, '.-r')
    ylabel('resultant force (kN)')
    xlabel('z (m)')
    box on
    
    legend('V_x', ...
           'V_y', ...
           'V_z', ...
           'Location', 'best')
       
    subplot(2,1,2);
    hold on
    plot(BLADE.zSec, ResLoads.Mx .* convMoment, 'o-b')
    plot(BLADE.zSec, ResLoads.My .* convMoment, 's-g')
    plot(BLADE.zSec, ResLoads.Mz .* convMoment, '.-r')
    ylabel('resultant moment (kN-m)')
    xlabel('z (m)')
    box on
    
    legend('M_x', ...
           'M_y', ...
           'M_z', ...
           'Location', 'best')
       
    if OUT.SAVE_PLOTS
        savePlots(fig, figTitle, OUT.SAVE_FIG_FMT)
    end
end

end % function plotLoads

