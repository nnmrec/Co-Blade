function plotModes(SIM, ANLS, Modes, OUT)

cmap = hsv(ANLS.N_MODES);  % creates a N_MODES-by-3 set of colors from the HSV colormap

if isinf(Modes.natFreq)
    fprintf(1,'Error: The BModes analysis failed, cannot create output plots.\n');
    return
end

% plot modal displacements
if OUT.PLOT_MODE_D
    figTitle = [SIM.case ' flapwise modal displacement'];
    fig = figure('name', figTitle, ...
                 'color', 'white', ...
                 'units','normalized',...
                 'outerposition',[0.1 0.1 0.8 0.8]);
    
    for n = 1:ANLS.N_MODES 
        hold on
        plot(Modes.span_loc{n}, Modes.flap_disp{n}, '.-', 'Color', cmap(n,:));
        leg{n} = ['mode ' num2str(n) ': f_N = ' num2str(Modes.natFreq(n),'%4.1f') ' Hz'];
    end
    xlabel('span location, (z/R)')
    ylabel('flapwise modal displacement, (-)')
    legend(leg,'Location','Best')
    box on
    if OUT.SAVE_PLOTS
        savePlots(fig, figTitle, OUT.SAVE_FIG_FMT)
    end
    %%
    figTitle = [SIM.case ' edgewise modal displacement'];
    fig = figure('name', figTitle, ...
                 'color', 'white', ...
                 'units','normalized',...
                 'outerposition',[0.1 0.1 0.8 0.8]);
    
    for n = 1:ANLS.N_MODES 
        hold on
        plot(Modes.span_loc{n}, Modes.lag_disp{n}, '.-', 'Color', cmap(n,:));
        leg{n} = ['mode ' num2str(n) ': f_N = ' num2str(Modes.natFreq(n),'%4.1f') ' Hz'];
    end
    xlabel('span location, (z/R)')
    ylabel('edgewise modal displacement, (-)')
    legend(leg,'Location','Best')
    box on
    if OUT.SAVE_PLOTS
        savePlots(fig, figTitle, OUT.SAVE_FIG_FMT)
    end
    %%
    figTitle = [SIM.case ' modal twist'];
    fig = figure('name', figTitle, ...
                 'color', 'white', ...
                 'units','normalized',...
                 'outerposition',[0.1 0.1 0.8 0.8]);
    
    for n = 1:ANLS.N_MODES 
        hold on
        plot(Modes.span_loc{n}, Modes.twist{n}, '.-', 'Color', cmap(n,:));
        leg{n} = ['mode ' num2str(n) ': f_N = ' num2str(Modes.natFreq(n),'%4.1f') ' Hz'];
    end
    xlabel('span location, (z/R)')
    ylabel('modal twist, (rad)')
    legend(leg,'Location','Best')
    box on
    if OUT.SAVE_PLOTS
        savePlots(fig, figTitle, OUT.SAVE_FIG_FMT)
    end
end

% plot modal slopes
if OUT.PLOT_MODE_S
    figTitle = [SIM.case ' flapwise modal slope'];
    fig = figure('name', figTitle, ...
                 'color', 'white', ...
                 'units','normalized',...
                 'outerposition',[0.1 0.1 0.8 0.8]);
    
    for n = 1:ANLS.N_MODES 
        hold on
        plot(Modes.span_loc{n}, Modes.flap_slope{n}, '.-', 'Color', cmap(n,:));
        leg{n} = ['mode ' num2str(n) ': f_N = ' num2str(Modes.natFreq(n),'%4.1f') ' Hz'];
    end
    xlabel('span location, (z/R)')
    ylabel('flapwise modal slope, (rad)')
    legend(leg,'Location','Best')
    box on
    if OUT.SAVE_PLOTS
        savePlots(fig, figTitle, OUT.SAVE_FIG_FMT)
    end
    %%
    figTitle = [SIM.case ' edgewise modal slope'];
    fig = figure('name', figTitle, ...
                 'color', 'white', ...
                 'units','normalized',...
                 'outerposition',[0.1 0.1 0.8 0.8]);
    
    for n = 1:ANLS.N_MODES 
        hold on
        plot(Modes.span_loc{n}, Modes.lag_slope{n}, '.-', 'Color', cmap(n,:));
        leg{n} = ['mode ' num2str(n) ': f_N = ' num2str(Modes.natFreq(n),'%4.1f') ' Hz'];
    end
    xlabel('span location, (z/R)')
    ylabel('edgewise modal slope, (rad)')
    legend(leg,'Location','Best')
    box on
    if OUT.SAVE_PLOTS
        savePlots(fig, figTitle, OUT.SAVE_FIG_FMT)
    end
    
end
  
end % function plotModes

