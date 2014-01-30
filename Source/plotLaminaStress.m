function plotLaminaStress(SIM, BLADE, WEB, MATS, Panel, LaminaSS, OUT)

%% re-assign some structure variable names (for convenience)
NUM_SEC    = BLADE.NUM_SEC;
BLD_LENGTH = BLADE.BLD_LENGTH;
aeroTwst   = BLADE.aeroTwst;
nWebs      = WEB.nWebs;
matID      = MATS.matID;    
matName    = MATS.matName; 

%%
convStress = 1e-6;          % conversion factor for stress (converts Pa to MPa)
d2r        = pi/180;        % convert degrees to radians

matID_used = [];            % material identification numbers for materials that were actually used
for i = 1:NUM_SEC
    matID_t = unique( cat(1, Panel(i).Top.matID{:}) );
    matID_b = unique( cat(1, Panel(i).Bot.matID{:}) );
    if isempty(Panel(i).Web.matID)
        matID_w = [];
    else
        matID_w = unique( cat(1, Panel(i).Web.matID{:}) );
    end
    
    matID_used = unique( [matID_used; matID_t; matID_b; matID_w] );
end

nMats = numel(matID);   % number of materials listed in the materials file
cmap  = lines(nMats);	% creates a nMats-by-3 set of colors from the lines colormap
az    = 45;             % azimuth of view for the 3D plots
el    = 10;             % elevation of the view for 3D plots

%% 1st principal stresses
if OUT.PLOT_S11
    figTitle = [SIM.case ' lamina 1st principal stresses'];
    fig = figure('name', figTitle, ...
                 'color', 'white', ...
                 'units','normalized',...
                 'outerposition',[0.1 0.1 0.8 0.8]);

    % create axes
    ha    = zeros(4,1);
    ha(1) = axes('position', [0.10  0.1  0.2  0.8]);
    ha(2) = axes('position', [0.35  0.1  0.2  0.8]);
    ha(3) = axes('position', [0.60  0.1  0.2  0.8]);
    ha(4) = axes('position', [0.8   0.1  0.2  0.8]);

    for i = 1:NUM_SEC

        c_vec = [cos(aeroTwst(i)*d2r); sin(aeroTwst(i)*d2r)];           % unit vector pointing in direction of chordline towards the TE

        % top panels
        subplot(ha(1))
        le    = [Panel(i).Top.xs{1}(1);     Panel(i).Top.ys{1}(1)];     % position vector of the LE
        te    = [Panel(i).Top.xs{end}(end); Panel(i).Top.ys{end}(end)]; % position vector of the TE
        c     = norm(te - le);                                          % "chord length"
        for j = 1:Panel(i).Top.nPanels  
            r   = [Panel(i).Top.xs{j} - le(1), Panel(i).Top.ys{j} - le(2)];
            x_c = (r * c_vec) ./ c;
            z_L = Panel(i).Top.zs{j} ./ BLD_LENGTH;

            nLam = Panel(i).Top.nLam(j);
            for k = 1:nLam
                hold on
                plot3(z_L, x_c, LaminaSS(i).Top(j).s_11_maxabs{k} .* convStress, 'Color', cmap(Panel(i).Top.matID{j}(k),:))
            end
        end

        % bottom panels
        subplot(ha(2))
        le    = [Panel(i).Bot.xs{1}(1);     Panel(i).Bot.ys{1}(1)];     % position vector of the LE
        te    = [Panel(i).Bot.xs{end}(end); Panel(i).Bot.ys{end}(end)]; % position vector of the TE
        c     = norm(te - le);                                          % "chord length"
        for j = 1:Panel(i).Bot.nPanels  
            r   = [Panel(i).Bot.xs{j} - le(1), Panel(i).Bot.ys{j} - le(2)];
            x_c = (r * c_vec) ./ c;
            z_L = Panel(i).Bot.zs{j} ./ BLD_LENGTH;

            nLam = Panel(i).Bot.nLam(j);
            for k = 1:nLam
                hold on
                plot3(z_L, x_c, LaminaSS(i).Bot(j).s_11_maxabs{k} .* convStress, 'Color', cmap(Panel(i).Bot.matID{j}(k),:))
            end
        end

        % web panels
        subplot(ha(3))
        for j = 1:nWebs(i)

            w   = Panel(i).Web.s{j}(end);
            s_w = Panel(i).Web.s{j} ./ w;
            z_L = Panel(i).Web.zs{j} ./ BLD_LENGTH;

            nLam = Panel(i).Web.nLam(j);
            for k = 1:nLam
                hold on
                plot3(z_L, s_w, LaminaSS(i).Web(j).s_11_maxabs{k} .* convStress, 'Color', cmap(Panel(i).Web.matID{j}(k),:))
            end
        end

    end

    % set axis labels and scale
    ttl = {'Top Panels','Bottom Panels','Web Panels'};
    for i = 1:3
        xlabel(ha(i), 'z/L')
        if i == 3
            ylabel(ha(i), 's/L_w')
        else
            ylabel(ha(i), 'x/c')
        end
        if i == 1
            zlabel(ha(i), '\sigma_{11}, (MPa)')
        end
        xlim(ha(i), [0 1])
        ylim(ha(i), [0 1])
        title(ha(i), ttl(i));
        view(ha(i), az,el)    
    end

    % create legend
    subplot(ha(4))
    for i = 1:numel(matID_used)
       hold on
       plot(0, 0, 'Color', cmap(matID_used(i),:));  % just need to create dummy data points, so a legend can be created
    end
    legend(matName(matID_used), 'Location', 'West')
    axis off
    if OUT.SAVE_PLOTS
        savePlots(fig, figTitle, OUT.SAVE_FIG_FMT)
    end

end
%% 2nd principal stresses
if OUT.PLOT_S22
    figTitle = [SIM.case ' lamina 2nd principal stresses'];
    fig = figure('name', figTitle, ...
                 'color', 'white', ...
                 'units','normalized',...
                 'outerposition',[0.1 0.1 0.8 0.8]);

    % create axes
    ha    = zeros(4,1);
    ha(1) = axes('position', [0.10  0.1  0.2  0.8]);
    ha(2) = axes('position', [0.35  0.1  0.2  0.8]);
    ha(3) = axes('position', [0.60  0.1  0.2  0.8]);
    ha(4) = axes('position', [0.8   0.1  0.2  0.8]);

    for i = 1:NUM_SEC

        c_vec = [cos(aeroTwst(i)*d2r); sin(aeroTwst(i)*d2r)];           % unit vector pointing in direction of chordline towards the TE

        % top panels
        subplot(ha(1))
        le    = [Panel(i).Top.xs{1}(1);     Panel(i).Top.ys{1}(1)];     % position vector of the LE
        te    = [Panel(i).Top.xs{end}(end); Panel(i).Top.ys{end}(end)]; % position vector of the TE
        c     = norm(te - le);                                          % "chord length"
        for j = 1:Panel(i).Top.nPanels  
            r   = [Panel(i).Top.xs{j} - le(1), Panel(i).Top.ys{j} - le(2)];
            x_c = (r * c_vec) ./ c;
            z_L = Panel(i).Top.zs{j} ./ BLD_LENGTH;

            nLam = Panel(i).Top.nLam(j);
            for k = 1:nLam
                hold on
                plot3(z_L, x_c, LaminaSS(i).Top(j).s_22_maxabs{k} .* convStress, 'Color', cmap(Panel(i).Top.matID{j}(k),:))
            end
        end

        % bottom panels
        subplot(ha(2))
        le    = [Panel(i).Bot.xs{1}(1);     Panel(i).Bot.ys{1}(1)];     % position vector of the LE
        te    = [Panel(i).Bot.xs{end}(end); Panel(i).Bot.ys{end}(end)]; % position vector of the TE
        c     = norm(te - le);                                          % "chord length"
        for j = 1:Panel(i).Bot.nPanels  
            r   = [Panel(i).Bot.xs{j} - le(1), Panel(i).Bot.ys{j} - le(2)];
            x_c = (r * c_vec) ./ c;
            z_L = Panel(i).Bot.zs{j} ./ BLD_LENGTH;

            nLam = Panel(i).Bot.nLam(j);
            for k = 1:nLam
                hold on
                plot3(z_L, x_c, LaminaSS(i).Bot(j).s_22_maxabs{k} .* convStress, 'Color', cmap(Panel(i).Bot.matID{j}(k),:))
            end
        end

        % web panels
        subplot(ha(3))
        for j = 1:nWebs(i)

            w   = Panel(i).Web.s{j}(end);
            s_w = Panel(i).Web.s{j} ./ w;
            z_L = Panel(i).Web.zs{j} ./ BLD_LENGTH;

            nLam = Panel(i).Web.nLam(j);
            for k = 1:nLam
                hold on
                plot3(z_L, s_w, LaminaSS(i).Web(j).s_22_maxabs{k} .* convStress, 'Color', cmap(Panel(i).Web.matID{j}(k),:))
            end
        end

    end

    % set axis labels and scale
    ttl = {'Top Panels','Bottom Panels','Web Panels'};
    for i = 1:3
        xlabel(ha(i), 'z/L')
        if i == 3
            ylabel(ha(i), 's/L_w')
        else
            ylabel(ha(i), 'x/c')
        end
        if i == 1
            zlabel(ha(i), '\sigma_{22}, (MPa)')
        end
        xlim(ha(i), [0 1])
        ylim(ha(i), [0 1])
        title(ha(i), ttl(i));
        view(ha(i), az,el)    
    end

    % create legend
    subplot(ha(4))
    for i = 1:numel(matID_used)
       hold on
       plot(0, 0, 'Color', cmap(matID_used(i),:));  % just need to create dummy data points, so a legend can be created
    end
    legend(matName(matID_used), 'Location', 'West')
    axis off
    if OUT.SAVE_PLOTS
        savePlots(fig, figTitle, OUT.SAVE_FIG_FMT)
    end

end
%% principal shear stresses
if OUT.PLOT_S12
    figTitle = [SIM.case ' lamina principal shear stresses'];
    fig = figure('name', figTitle, ...
                 'color', 'white', ...
                 'units','normalized',...
                 'outerposition',[0.1 0.1 0.8 0.8]);

    % create axes
    ha    = zeros(4,1);
    ha(1) = axes('position', [0.10  0.1  0.2  0.8]);
    ha(2) = axes('position', [0.35  0.1  0.2  0.8]);
    ha(3) = axes('position', [0.60  0.1  0.2  0.8]);
    ha(4) = axes('position', [0.8   0.1  0.2  0.8]);

    for i = 1:NUM_SEC

        c_vec = [cos(aeroTwst(i)*d2r); sin(aeroTwst(i)*d2r)];           % unit vector pointing in direction of chordline towards the TE

        % top panels
        subplot(ha(1))
        le    = [Panel(i).Top.xs{1}(1);     Panel(i).Top.ys{1}(1)];     % position vector of the LE
        te    = [Panel(i).Top.xs{end}(end); Panel(i).Top.ys{end}(end)]; % position vector of the TE
        c     = norm(te - le);                                          % "chord length"
        for j = 1:Panel(i).Top.nPanels  
            r   = [Panel(i).Top.xs{j} - le(1), Panel(i).Top.ys{j} - le(2)];
            x_c = (r * c_vec) ./ c;
            z_L = Panel(i).Top.zs{j} ./ BLD_LENGTH;

            nLam = Panel(i).Top.nLam(j);
            for k = 1:nLam
                hold on
                plot3(z_L, x_c, LaminaSS(i).Top(j).s_12_maxabs{k} .* convStress, 'Color', cmap(Panel(i).Top.matID{j}(k),:))
            end
        end

        % bottom panels
        subplot(ha(2))
        le    = [Panel(i).Bot.xs{1}(1);     Panel(i).Bot.ys{1}(1)];     % position vector of the LE
        te    = [Panel(i).Bot.xs{end}(end); Panel(i).Bot.ys{end}(end)]; % position vector of the TE
        c     = norm(te - le);                                          % "chord length"
        for j = 1:Panel(i).Bot.nPanels  
            r   = [Panel(i).Bot.xs{j} - le(1), Panel(i).Bot.ys{j} - le(2)];
            x_c = (r * c_vec) ./ c;
            z_L = Panel(i).Bot.zs{j} ./ BLD_LENGTH;

            nLam = Panel(i).Bot.nLam(j);
            for k = 1:nLam
                hold on
                plot3(z_L, x_c, LaminaSS(i).Bot(j).s_12_maxabs{k} .* convStress, 'Color', cmap(Panel(i).Bot.matID{j}(k),:))
            end
        end

        % web panels
        subplot(ha(3))
        for j = 1:nWebs(i)

            w   = Panel(i).Web.s{j}(end);
            s_w = Panel(i).Web.s{j} ./ w;
            z_L = Panel(i).Web.zs{j} ./ BLD_LENGTH;

            nLam = Panel(i).Web.nLam(j);
            for k = 1:nLam
                hold on
                plot3(z_L, s_w, LaminaSS(i).Web(j).s_12_maxabs{k} .* convStress, 'Color', cmap(Panel(i).Web.matID{j}(k),:))
            end
        end

    end

    % set axis labels and scale
    ttl = {'Top Panels','Bottom Panels','Web Panels'};
    for i = 1:3
        xlabel(ha(i), 'z/L')
        if i == 3
            ylabel(ha(i), 's/L_w')
        else
            ylabel(ha(i), 'x/c')
        end
        if i == 1
            zlabel(ha(i), '|\tau_{12}|, (MPa)')
        end
        xlim(ha(i), [0 1])
        ylim(ha(i), [0 1])
        title(ha(i), ttl(i));
        view(ha(i), az,el)    
    end

    % create legend
    subplot(ha(4))
    for i = 1:numel(matID_used)
       hold on
       plot(0, 0, 'Color', cmap(matID_used(i),:));  % just need to create dummy data points, so a legend can be created
    end
    legend(matName(matID_used), 'Location', 'West')
    axis off
    if OUT.SAVE_PLOTS
        savePlots(fig, figTitle, OUT.SAVE_FIG_FMT)
    end
    
end

end % function plotLaminaStress
