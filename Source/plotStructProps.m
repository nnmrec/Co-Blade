function plotStructProps(SIM, BLADE, StrProps, OUT)

if OUT.PLOT_MASS_DEN
    figTitle = [SIM.case ' mass density'];
    fig = figure('name', figTitle, ...
                 'color', 'white', ...
                 'units','normalized',...
                 'outerposition',[0.1 0.1 0.8 0.8]); 
             
    plot(BLADE.zSec, StrProps.mass_den,  'o-k')
    ylabel('mass density (kg/m)')
    xlabel('z (m)')
    box on
    
    if OUT.SAVE_PLOTS
        savePlots(fig, figTitle, OUT.SAVE_FIG_FMT)
    end
end
%% 
if OUT.PLOT_PRIN_ANG
    figTitle = [SIM.case ' principal angles'];
    fig = figure('name', figTitle, ...
                 'color', 'white', ...
                 'units','normalized',...
                 'outerposition',[0.1 0.1 0.8 0.8]);
             
    hold on
    plot(BLADE.zSec, StrProps.iner_tw, 's-r')
    plot(BLADE.zSec, StrProps.cent_tw, 'o-b')
    plot(BLADE.zSec, StrProps.elas_tw, '.-g')
    ylabel('flapwise principal angle (deg)')
    xlabel('z (m)')
    box on

    legend('w.r.t. mass center', ...
           'w.r.t. tension center', ...
           'w.r.t. shear center')
       
    if OUT.SAVE_PLOTS
        savePlots(fig, figTitle, OUT.SAVE_FIG_FMT)
    end
end
%%
if OUT.PLOT_AT_STFF
    figTitle = [SIM.case ' axial and torsional stiffness'];
    fig = figure('name', figTitle, ...
                 'color', 'white', ...
                 'units','normalized',...
                 'outerposition',[0.1 0.1 0.8 0.8]);   
    
    subplot(2,1,1)
    plot(BLADE.zSec, StrProps.axial_stff,  'o-r')
    ylabel('axial stiffness (N)')
    xlabel('z (m)')
    box on
    
    subplot(2,1,2)
    plot(BLADE.zSec, StrProps.tor_stff,  'o-b')
    ylabel('torsional stiffness (N-m^2)')
    xlabel('z (m)')
    box on 
    
    if OUT.SAVE_PLOTS
        savePlots(fig, figTitle, OUT.SAVE_FIG_FMT)
    end
end
%%
if OUT.PLOT_BSTFF
    figTitle = [SIM.case ' bending stiffness'];
    fig = figure('name', figTitle, ...
                 'color', 'white', ...
                 'units','normalized',...
                 'outerposition',[0.1 0.1 0.8 0.8]);
              
    subplot(3,1,1)
    hold on
    plot(BLADE.zSec, StrProps.flapEI_cm, 's-r')
    plot(BLADE.zSec, StrProps.flapEI_tc, 'o-b')
    plot(BLADE.zSec, StrProps.flapEI_sc, '.-g')
    ylabel('EI_{flap} (N-m^2)')
    xlabel('z (m)')
    box on
    legend('w.r.t. mass center', ...
           'w.r.t. tension center', ...
           'w.r.t. shear center', ...
           'Location', 'EastOutside')
    pos = get(gca, 'Position');
    
    subplot(3,1,2)
    hold on
    plot(BLADE.zSec, StrProps.edgeEI_cm, 's-r')
    plot(BLADE.zSec, StrProps.edgeEI_tc, 'o-b')
    plot(BLADE.zSec, StrProps.edgeEI_sc, '.-g')
    ylabel('EI_{edge} (N-m^2)')
    xlabel('z (m)')
    box on
    legend('w.r.t. mass center', ...
           'w.r.t. tension center', ...
           'w.r.t. shear center', ...
           'Location', 'EastOutside') 
    p = get(gca, 'Position');
    set(gca, 'Position', [p(1) p(2) pos(3) p(4)]);
    
    subplot(3,1,3)
    hold on
    plot(BLADE.zSec, StrProps.EIx, 's-r')
    plot(BLADE.zSec, StrProps.EIy, 'o-b')
    plot(BLADE.zSec, StrProps.EIxy, '.-g')
    ylabel('EI (N-m^2)')
    xlabel('z (m)')
    box on
    legend('(EI)_x', ...
           '(EI)_y', ...
           '(EI)_{xy}', ...
           'Location', 'EastOutside')  
    p = get(gca, 'Position');
    set(gca, 'Position', [p(1) p(2) pos(3) p(4)]);
    
    if OUT.SAVE_PLOTS
        savePlots(fig, figTitle, OUT.SAVE_FIG_FMT)
    end
end
%%
if OUT.PLOT_INER
    figTitle = [SIM.case ' mass moment of inertia'];
    fig = figure('name', figTitle, ...
                 'color', 'white', ...
                 'units','normalized',...
                 'outerposition',[0.1 0.1 0.8 0.8]);
              
    subplot(3,1,1)
    hold on
    plot(BLADE.zSec, StrProps.flapIner_cm, 's-r')
    plot(BLADE.zSec, StrProps.flapIner_tc, 'o-b')
    plot(BLADE.zSec, StrProps.flapIner_sc, '.-g')
    ylabel('mI_{flap} (kg-m)')
    xlabel('z (m)')
    box on
    legend('w.r.t. mass center', ...
           'w.r.t. tension center', ...
           'w.r.t. shear center', ...
           'Location', 'EastOutside')
    pos = get(gca, 'Position');
    
    subplot(3,1,2)
    hold on
    plot(BLADE.zSec, StrProps.edgeIner_cm, 's-r')
    plot(BLADE.zSec, StrProps.edgeIner_tc, 'o-b')
    plot(BLADE.zSec, StrProps.edgeIner_sc, '.-g')
    ylabel('mI_{edge} (kg-m)')
    xlabel('z (m)')
    box on
    legend('w.r.t. mass center', ...
           'w.r.t. tension center', ...
           'w.r.t. shear center', ...
           'Location', 'EastOutside') 
    p = get(gca, 'Position');
    set(gca, 'Position', [p(1) p(2) pos(3) p(4)]);
    
    subplot(3,1,3)
    hold on
    plot(BLADE.zSec, StrProps.mIx, 's-r')
    plot(BLADE.zSec, StrProps.mIy, 'o-b')
    plot(BLADE.zSec, StrProps.mIxy, '.-g')
    ylabel('mI (kg-m)')
    xlabel('z (m)')
    box on
    legend('mI_x', ...
           'mI_y', ...
           'mI_{xy}', ...
           'Location', 'EastOutside') 
    p = get(gca, 'Position');
    set(gca, 'Position', [p(1) p(2) pos(3) p(4)]);
    
    if OUT.SAVE_PLOTS
        savePlots(fig, figTitle, OUT.SAVE_FIG_FMT)
    end
end
%%
if OUT.PLOT_CENTERS
    figTitle = [SIM.case ' center of mass, tension, shear'];
    fig = figure('name', figTitle, ...
                 'color', 'white', ...
                 'units','normalized',...
                 'outerposition',[0.1 0.1 0.8 0.8]);
              
    subplot(3,1,1)
    hold on
    plot(BLADE.zSec, StrProps.x_cm, 's-r')
    plot(BLADE.zSec, StrProps.x_tc, 'o-b')
    plot(BLADE.zSec, StrProps.x_sc, '.-g')
    ylabel('x (m)')
    xlabel('z (m)')
    box on
    legend('mass center', ...
           'tension center', ...
           'shear center', ...
           'Location', 'EastOutside')
    pos = get(gca, 'Position');
    
    subplot(3,1,2)
    hold on
    plot(BLADE.zSec, StrProps.y_cm, 's-r')
    plot(BLADE.zSec, StrProps.y_tc, 'o-b')
    plot(BLADE.zSec, StrProps.y_sc, '.-g')
    ylabel('y (m)')
    xlabel('z (m)')
    box on
    p = get(gca, 'Position');
    set(gca, 'Position', [p(1) p(2) pos(3) p(4)]);
    
    subplot(3,1,3)
    hold on
    plot(BLADE.zSec, StrProps.cm_offst, 's-r')
    plot(BLADE.zSec, StrProps.tc_offst, 'o-b')
    plot(BLADE.zSec, StrProps.sc_offst, '.-g')
    ylabel('chord offset (x/c)')
    xlabel('z (m)')
    box on  
    p = get(gca, 'Position');
    set(gca, 'Position', [p(1) p(2) pos(3) p(4)]);
    
    if OUT.SAVE_PLOTS
        savePlots(fig, figTitle, OUT.SAVE_FIG_FMT)
    end
end
                    
end %function plotStructProps