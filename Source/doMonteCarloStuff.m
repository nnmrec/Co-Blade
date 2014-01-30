function MC = doMonteCarloStuff(SIM, ANLS, ENV, BLADE, WEB, MATS, AF, Coord)

%% Analysis portion of uncertian material properties...

% generate a set of random material properties, according to a normal PDF   
rnd_E11     = zeros(MATS.numMats, ANLS.SAMPLE_SIZE);
rnd_E22     = zeros(MATS.numMats, ANLS.SAMPLE_SIZE);
rnd_G12     = zeros(MATS.numMats, ANLS.SAMPLE_SIZE);
rnd_nu12    = zeros(MATS.numMats, ANLS.SAMPLE_SIZE);
rnd_density = zeros(MATS.numMats, ANLS.SAMPLE_SIZE);
for n = 1:MATS.numMats
    rnd_E11(n,:)     = randraw('norm', [MATS.E11(n),     ANLS.MAT_COV*MATS.E11(n)],     ANLS.SAMPLE_SIZE);
    rnd_E22(n,:)     = randraw('norm', [MATS.E22(n),     ANLS.MAT_COV*MATS.E22(n)],     ANLS.SAMPLE_SIZE);
    rnd_G12(n,:)     = randraw('norm', [MATS.G12(n),     ANLS.MAT_COV*MATS.G12(n)],     ANLS.SAMPLE_SIZE);
    rnd_nu12(n,:)    = randraw('norm', [MATS.nu12(n),    ANLS.MAT_COV*MATS.nu12(n)],    ANLS.SAMPLE_SIZE);
    rnd_density(n,:) = randraw('norm', [MATS.density(n), ANLS.MAT_COV*MATS.density(n)], ANLS.SAMPLE_SIZE);
end

% replace the last sample with the original mean values (so that any plots created by the user correspond to the nominal value, rather than some unknown value!)
% this should not affect any of the statistics if the sample size is large enough
rnd_E11(:,end)     = MATS.E11;
rnd_E22(:,end)     = MATS.E22;
rnd_G12(:,end)     = MATS.G12;
rnd_nu12(:,end)    = MATS.nu12;
rnd_density(:,end) = MATS.density;

% initialize random output variables
tipDeflect       = zeros(ANLS.SAMPLE_SIZE, 1);
bladeMass        = zeros(ANLS.SAMPLE_SIZE, 1);
maxStressFC_mat1 =  cell(ANLS.SAMPLE_SIZE, 1); % the maximum absolute stress (within the spatial domain of this material ID)
maxStressFC_mat2 =  cell(ANLS.SAMPLE_SIZE, 1);
maxStressFC_mat3 =  cell(ANLS.SAMPLE_SIZE, 1);
maxStressFC_mat4 =  cell(ANLS.SAMPLE_SIZE, 1);
maxStressFC_mat5 =  cell(ANLS.SAMPLE_SIZE, 1);
maxStressFC_mat6 =  cell(ANLS.SAMPLE_SIZE, 1);
maxStressFC_mat7 =  cell(ANLS.SAMPLE_SIZE, 1);
maxStressFC_mat8 =  cell(ANLS.SAMPLE_SIZE, 1);
for nSample = 1:ANLS.SAMPLE_SIZE 
    % modify the nominal material properties
    MATS.E11     =     rnd_E11(:,nSample);
    MATS.E22     =     rnd_E22(:,nSample);
    MATS.G12     =     rnd_G12(:,nSample);
    MATS.nu12    =    rnd_nu12(:,nSample);
    MATS.density = rnd_density(:,nSample);

    % Read the pre-defined laminate data from the laminate input files
    [SECNODES LamData] = readLaminateData(SIM, BLADE, WEB, MATS);

    % execute the structural analysis
    [Panel, StrProps, AppLoads, ResLoads, Disp, NormS, ShearS, Buckle, MidPlane, LaminaSS, Modes] ...
    = structAnalysis(SIM, ANLS, ENV, BLADE, WEB, MATS, LamData, AF, SECNODES, Coord);

    % collect random output variables
    tipDeflect(nSample)       = Disp.tipDeflect;
    bladeMass(nSample)        = StrProps.bladeMass;
    maxBuckleCrit(nSample)    = getExtremeBuckle(BLADE, Buckle, OPT);
    maxStressFC_mat1{nSample} = getExtremeStress(1, BLADE, Panel, LaminaSS);
    maxStressFC_mat2{nSample} = getExtremeStress(2, BLADE, Panel, LaminaSS);
    maxStressFC_mat3{nSample} = getExtremeStress(3, BLADE, Panel, LaminaSS);
    maxStressFC_mat4{nSample} = getExtremeStress(4, BLADE, Panel, LaminaSS);
    maxStressFC_mat5{nSample} = getExtremeStress(5, BLADE, Panel, LaminaSS);
    maxStressFC_mat6{nSample} = getExtremeStress(6, BLADE, Panel, LaminaSS);
    maxStressFC_mat7{nSample} = getExtremeStress(7, BLADE, Panel, LaminaSS);
    maxStressFC_mat8{nSample} = getExtremeStress(8, BLADE, Panel, LaminaSS);

    fprintf(1, 'Running Structural Analysis: material properties sample: %g of %g \n', nSample, ANLS.SAMPLE_SIZE);

end


%% now let's create some plots to make sense of all this...
nBins = 0.025*ANLS.SAMPLE_SIZE;

plotMonteCarloStressFC(SIM, nBins, ANLS.SAMPLE_SIZE, cell2mat(maxStressFC_mat1), 'MATERIAL: “NCT307-D1-E300” triaxial E-glass              COMPONENT: root build-up & transition')     % triaxial E-glass
plotMonteCarloStressFC(SIM, nBins, ANLS.SAMPLE_SIZE, cell2mat(maxStressFC_mat2), 'MATERIAL: "NB307-D1-7781-497A" double-bias E-glass       COMPONENT: blade shell')                    % double-bias E-glass
plotMonteCarloStressFC(SIM, nBins, ANLS.SAMPLE_SIZE, cell2mat(maxStressFC_mat7), 'MATERIAL: "NB307-D1-7781-497A" double-bias E-glass       COMPONENT: webs shell')                     % double-bias E-glass
plotMonteCarloStressFC(SIM, nBins, ANLS.SAMPLE_SIZE, cell2mat(maxStressFC_mat3), 'MATERIAL: "NCT307-D1-34-600" unidirectional carbon fiber COMPONENT: spar cap')                       % carbon fiber
plotMonteCarloStressFC(SIM, nBins, ANLS.SAMPLE_SIZE, cell2mat(maxStressFC_mat4), 'MATERIAL: "Corecell M-Foam M200" structural foam         COMPONENT: spar cap core')                  % core
plotMonteCarloStressFC(SIM, nBins, ANLS.SAMPLE_SIZE, cell2mat(maxStressFC_mat5), 'MATERIAL: "Corecell M-Foam M200" structural foam         COMPONENT: LEP core')                       % core
plotMonteCarloStressFC(SIM, nBins, ANLS.SAMPLE_SIZE, cell2mat(maxStressFC_mat6), 'MATERIAL: "Corecell M-Foam M200" structural foam         COMPONENT: TEP core')                       % core
plotMonteCarloStressFC(SIM, nBins, ANLS.SAMPLE_SIZE, cell2mat(maxStressFC_mat8), 'MATERIAL: "Corecell M-Foam M200" structural foam         COMPONENT: webs core')                      % core

plotMonteCarloBuckleFC(SIM, nBins, ANLS.SAMPLE_SIZE, maxBuckleCrit, '')


figTitle = [SIM.case ' Monte Carlo analysis'];
fig = figure('name',          figTitle, ...
             'color',         'white', ...
             'units',         'normalized',...
             'outerposition', [0.1 0.1 0.8 0.8]);
set(fig, 'color', 'white')

[count1, x1] = hist(tipDeflect, nBins);
bin_width1   = diff(x1(1:2));   
pdf1         = count1 / ANLS.SAMPLE_SIZE / bin_width1;     
z1           = (x1 - mean(tipDeflect))/std(tipDeflect); 
norm_pdf1    = std(tipDeflect) .* pdf1;
fit_p1       = fitdist(tipDeflect, 'normal');
xlim         = 1.2;
idx1         = linspace(min(0, min(tipDeflect)), max(xlim, max(tipDeflect)), ANLS.SAMPLE_SIZE);      
hold on     
plot(x1, norm_pdf1, '.r')
plot(idx1, std(tipDeflect) .* pdf('Normal', idx1, fit_p1.mu, fit_p1.sigma), '-r');
box on
xlabel('tip deflection penalty factor, p_7')
ylabel('normalized PDF')
fprintf(1, 'p_i     mean    std    cov    P_fail \n');
fprintf(1, '%s    %3.2f     %3.2f     %3.2f     %5.1f     \n', 'p7', mean(tipDeflect), std(tipDeflect), std(tipDeflect)/mean(tipDeflect), 100 * trapz(x1(x1>1), pdf1(x1>1)));

%%

figure
set(gcf, 'color', 'white')
hold on
plot_matID = 3;
% plot the values from the analytical normal distribution
x = linspace(min(rnd_E11(plot_matID,:)), max(rnd_E11(plot_matID,:)), 100);               
%     plot(x./1e9, y, 'b')
% plot the actual randomly generated values
[n,x] = hist(rnd_E11(plot_matID,:), nBins);
plot(x./1e9, n/ANLS.SAMPLE_SIZE/diff(x(1:2)), '.r')
plot(x./1e9, pdf('normal', x, mean(rnd_E11(plot_matID,:)), std(rnd_E11(plot_matID,:))), '-b')
xlabel('elastic modulus, E_{11} (GPa)')
ylabel('PDF')
box on

E11_mean = mean(rnd_E11(plot_matID,:));
E11_std  = std(rnd_E11(plot_matID,:));
E11_cov  = E11_std/E11_mean;
str = sprintf('mean = %0.2f (GPa), std = %0.2f (GPa), COV = %0.2f', E11_mean./1e9, E11_std./1e9, E11_cov);
title(str)

%% =========================================================        

tmp           = cell2mat(maxStressFC_mat3);
mat3_max_s_11 = [tmp.s_11_fc_C]';
%         
%         tipDeflect(nSample)     = Disp.tipDeflect;
%         bladeMass(nSample)      = trapzf(BLADE.zSec, StrProps.mass_den);
%         maxBuckleCrit(nSample)  = getExtremeBuckle(BLADE, Buckle, OPT);
%         maxStress_mat1{nSample} = getExtremeStress(1, BLADE, Panel, LaminaSS);
%         maxStress_mat2{nSample} = getExtremeStress(2, BLADE, Panel, LaminaSS);
%         maxStress_mat3{nSample} = getExtremeStress(3, BLADE, Panel, LaminaSS);
%         maxStress_mat4{nSample} = getExtremeStress(4, BLADE, Panel, LaminaSS);
%         maxStress_mat5{nSample} = getExtremeStress(5, BLADE, Panel, LaminaSS);
%         maxStress_mat6{nSample} = getExtremeStress(6, BLADE, Panel, LaminaSS);
%         maxStress_mat7{nSample} = getExtremeStress(7, BLADE, Panel, LaminaSS);
%         maxStress_mat8{nSample} = getExtremeStress(8, BLADE, Panel, LaminaSS);

%         
%         figure(21)
%         set(gcf, 'color', 'white')
%         [n,x] = hist(bladeMass, nBins);
%         pdf = n/ANLS.SAMPLE_SIZE/diff(x(1:2));
% %         plot(x./mean(bladeMass), pdf, '.r')
%         plot((x-mean(bladeMass))/std(bladeMass), pdf.*std(bladeMass), '.r')
%         xlabel('blade mass (m_{blade}), x/\mu_{x}')
%         ylabel('PDF')
%         box on 
%         bladeMass_mean = mean(bladeMass);
%         bladeMass_std  = std(bladeMass);
%         bladeMass_cov  = bladeMass_std/bladeMass_mean;
%         str = sprintf('mean = %0.2f (kg), std = %0.2f (kg), COV = %0.3f', bladeMass_mean, bladeMass_std, bladeMass_cov);
%         title(str)
%         xlim([0.4 2.4])
%         % =========================================================
%         figure(22)
%         set(gcf, 'color', 'white')
%         [n,x] = hist(tipDeflect, nBins);
%         pdf = n/ANLS.SAMPLE_SIZE/diff(x(1:2));
%         plot(x./mean(tipDeflect), pdf, '.r')
%         xlabel('tip deflection (\delta_{tip}), x/\mu_{x}')
%         ylabel('PDF')
%         box on 
%         tipDeflect_mean = mean(tipDeflect);
%         tipDeflect_std  = std(tipDeflect);
%         tipDeflect_cov  = tipDeflect_std/tipDeflect_mean;
%         str = sprintf('mean = %0.2f (m), std = %0.2f (m), COV = %0.3f', tipDeflect_mean, tipDeflect_std, tipDeflect_cov);
%         title(str)
%         xlim([0.4 2.4])
%         % =========================================================
%         figure(23)
%         set(gcf, 'color', 'white')
%         var = [maxBuckleCrit.TEP]';
%         [n,x] = hist(var, nBins);
%         pdf = n/ANLS.SAMPLE_SIZE/diff(x(1:2));
%         plot(x./mean(var), pdf, '.r')
%         xlabel('max buckling criteria in TEP (R_{b,TEP}), x/\mu_x')
%         ylabel('PDF')
%         box on 
%         var_mean = mean(var);
%         var_std  = std(var);
%         var_cov  = var_std/var_mean;
%         str = sprintf('mean = %0.2f, std = %0.2f, COV = %0.3f', var_mean, var_std, var_cov);
%         title(str)
%         xlim([0.4 2.4])
%         % =========================================================
%         figure(24)
%         set(gcf, 'color', 'white')
%         var = [maxBuckleCrit.webs]';
%         [n,x] = hist(var, nBins);
%         pdf = n/ANLS.SAMPLE_SIZE/diff(x(1:2));
%         plot(x./mean(var), pdf, '.r')
%         xlabel('max buckling criteria in webs (R_{b,webs}), x/\mu_x')
%         ylabel('PDF')
%         box on 
%         var_mean = mean(var);
%         var_std  = std(var);
%         var_cov  = var_std/var_mean;
%         str = sprintf('mean = %0.2f, std = %0.2f, COV = %0.3f', var_mean, var_std, var_cov);
%         title(str)
%         xlim([0.4 2.4])
%         % =========================================================
%         figure(25)
%         set(gcf, 'color', 'white')
%         tmp           = cell2mat(maxStress_mat3);
%         mat3_max_s_11 = [tmp.s_11]';
%         [n,x] = hist(mat3_max_s_11, nBins);
%         pdf = n/ANLS.SAMPLE_SIZE/diff(x(1:2));
%         plot(x/mean(mat3_max_s_11), pdf, '.r')
%         xlabel('max stress in material #3 (\sigma_{11,mat3}), x/\mu_x')
%         ylabel('PDF')
%         box on 
%         mat3_max_s_11_mean = mean(mat3_max_s_11);
%         mat3_max_s_11_std  = std(mat3_max_s_11);
%         mat3_max_s_11_cov  = mat3_max_s_11_std/mat3_max_s_11_mean;
%         str = sprintf('mean = %0.2f (MPa), std = %0.2f (MPa), COV = %0.3f', mat3_max_s_11_mean/1e6, mat3_max_s_11_std/1e6, mat3_max_s_11_cov);
%         title(str)
%         xlim([0.4 2.4])
%         % =========================================================
%         figure(26)
%         set(gcf, 'color', 'white')
%         tmp           = cell2mat(maxStress_mat2);
%         mat2_max_s_11 = [tmp.s_11]';
%         [n,x] = hist(mat2_max_s_11, nBins);
%         pdf = n/ANLS.SAMPLE_SIZE/diff(x(1:2));
%         plot(x/mean(mat2_max_s_11), pdf, '.r')
%         xlabel('max stress in material #2 (\sigma_{11,mat2}), x/\mu_x')
%         ylabel('PDF')
%         box on 
%         mat2_max_s_11_mean = mean(mat2_max_s_11);
%         mat2_max_s_11_std  = std(mat2_max_s_11);
%         mat2_max_s_11_cov  = mat2_max_s_11_std/mat2_max_s_11_mean;
%         str = sprintf('mean = %0.2f (MPa), std = %0.2f (MPa), COV = %0.3f', mat2_max_s_11_mean/1e6, mat2_max_s_11_std/1e6, mat2_max_s_11_cov);
%         title(str)
%         xlim([0.4 2.4])
%         % =========================================================
%         figure(27)
%         set(gcf, 'color', 'white')
%         tmp           = cell2mat(maxStress_mat2);
%         mat2_max_s_12 = [tmp.s_12]';
%         [n,x] = hist(mat2_max_s_12, nBins);
%         pdf = n/ANLS.SAMPLE_SIZE/diff(x(1:2));
%         plot(x/mean(mat2_max_s_12), pdf, '.r')
%         xlabel('max stress in material #2 (\tau_{12,mat2}), x/\mu_x')
%         ylabel('PDF')
%         box on 
%         mat2_max_s_12_mean = mean(mat2_max_s_12);
%         mat2_max_s_12_std  = std(mat2_max_s_12);
%         mat2_max_s_12_cov  = mat2_max_s_12_std/mat2_max_s_12_mean;
%         str = sprintf('mean = %0.2f (MPa), std = %0.2f (MPa), COV = %0.3f', mat2_max_s_12_mean/1e6, mat2_max_s_12_std/1e6, mat2_max_s_12_cov);
%         title(str)
%         xlim([0.4 2.4])
%         % =========================================================
%         figure(28)
%         set(gcf, 'color', 'white')
%         tmp           = cell2mat(maxStress_mat7);
%         mat7_max_s_12 = [tmp.s_12]';
%         [n,x] = hist(mat7_max_s_12, nBins);
%         pdf = n/ANLS.SAMPLE_SIZE/diff(x(1:2));
%         plot(x/mean(mat7_max_s_12), pdf, '.r')
%         xlabel('max stress in material #7 (\tau_{12,mat7}), x/\mu_x')
%         ylabel('PDF')
%         box on 
%         mat7_max_s_12_mean = mean(mat7_max_s_12);
%         mat7_max_s_12_std  = std(mat7_max_s_12);
%         mat7_max_s_12_cov  = mat7_max_s_12_std/mat7_max_s_12_mean;
%         str = sprintf('mean = %0.2f (MPa), std = %0.2f (MPa), COV = %0.3f', mat7_max_s_12_mean/1e6, mat7_max_s_12_std/1e6, mat7_max_s_12_cov);
%         title(str)
%         xlim([0.4 2.4])
%         %% =========================================================
%         %     figure(6) % Tornado diagram (I think the PDFs show more information)
%         %     set(gcf, 'color', 'white')
%         %       
%         %     var1 = (mat3_max_s_11 ./ mat3_max_s_11_mean);
%         %     var2 = bladeMass  ./ bladeMass_mean;
%         %     var3 = tipDeflect ./ tipDeflect_mean;
%         %     var4 = var ./ var_mean;
%         %     X    = [1; 2; 3; 4];
%         %     Y    = [var1'; var2'; var3'; var4'];
%         % 
%         %     hold on
%         %     plot(var1, 1.*ones(ANLS.SAMPLE_SIZE,1), 'or')
%         %     plot(var2, 2.*ones(ANLS.SAMPLE_SIZE,1), 'ob')
%         %     plot(var3, 3.*ones(ANLS.SAMPLE_SIZE,1), 'ok')
%         %     plot(var4, 4.*ones(ANLS.SAMPLE_SIZE,1), 'om')
%         %     ylim([0 5])
%         %     xlabel('x / \mu_x')
%         %     title(sprintf('sample size = %g',ANLS.SAMPLE_SIZE))
%         %     box on
%         %     
%         %     % Save label positions 
%         %     xtl = get(gca,'xticklabel'); 
%         %     xt  = get(gca,'xtick'); 
%         %     xl = get(gca,'xlabel'); 
%         %     yl = get(gca,'ylabel'); 
%         %     orig_label_units = get( xl, 'Units'); 
%         %     set([xl yl], 'Units','pixels'); 
%         %     xp = get( xl, 'Position'); 
%         %     yp = get( yl, 'Position'); 
%         %     % Generate the text object ticks (changing labels' Position): 
%         %     [hx,hy] = format_ticks(gca,cellstr(xtl),{'\sigma_{11,mat3}','m_{blade}','\delta_{tip}','R_{b,max}'},xt,[1 2 3 4]);
%         %     % Restore original label Position and Units properties: 
%         %     set( xl, 'Position', xp); 
%         %     set( yl, 'Position', yp); 
%         %     set([xl yl], 'Units', orig_label_units);    
% 
%         %% =========================================================

mat1 = [rnd_E11(1,:)', rnd_E22(1,:)', rnd_G12(1,:)', rnd_nu12(1,:)', rnd_density(1,:)'];
mat2 = [rnd_E11(2,:)', rnd_E22(2,:)', rnd_G12(2,:)', rnd_nu12(2,:)', rnd_density(2,:)'];
mat3 = [rnd_E11(3,:)', rnd_E22(3,:)', rnd_G12(3,:)', rnd_nu12(3,:)', rnd_density(3,:)'];
mat4 = [rnd_E11(4,:)', rnd_E22(4,:)', rnd_G12(4,:)', rnd_nu12(4,:)', rnd_density(4,:)'];
mat5 = [rnd_E11(5,:)', rnd_E22(5,:)', rnd_G12(5,:)', rnd_nu12(5,:)', rnd_density(5,:)'];
mat6 = [rnd_E11(6,:)', rnd_E22(6,:)', rnd_G12(6,:)', rnd_nu12(6,:)', rnd_density(6,:)'];
mat7 = [rnd_E11(7,:)', rnd_E22(7,:)', rnd_G12(7,:)', rnd_nu12(7,:)', rnd_density(7,:)'];
mat8 = [rnd_E11(8,:)', rnd_E22(8,:)', rnd_G12(8,:)', rnd_nu12(8,:)', rnd_density(8,:)'];


%% bladeMass vs. {E_11 for all mats}, {E_22 for all mats}, ... and {density for all mats}
%         figure
%         set(gcf, 'color', 'white')
%         
%         R1_mat1 = cov( [bladeMass, mat1] );
%         R1_mat2 = cov( [bladeMass, mat2] );
%         R1_mat3 = cov( [bladeMass, mat3] );
%         R1_mat4 = cov( [bladeMass, mat4] );
%         R1_mat5 = cov( [bladeMass, mat5] );
%         R1_mat6 = cov( [bladeMass, mat6] );
%         R1_mat7 = cov( [bladeMass, mat7] );
%         R1_mat8 = cov( [bladeMass, mat8] );
% 
%         rho_bladeMass_to_mat1 = R1_mat1(1, 2:6) ./ (std(bladeMass) .* std(mat1));
%         rho_bladeMass_to_mat2 = R1_mat2(1, 2:6) ./ (std(bladeMass) .* std(mat2));
%         rho_bladeMass_to_mat3 = R1_mat3(1, 2:6) ./ (std(bladeMass) .* std(mat3));
%         rho_bladeMass_to_mat4 = R1_mat4(1, 2:6) ./ (std(bladeMass) .* std(mat4));
%         rho_bladeMass_to_mat5 = R1_mat5(1, 2:6) ./ (std(bladeMass) .* std(mat5));
%         rho_bladeMass_to_mat6 = R1_mat6(1, 2:6) ./ (std(bladeMass) .* std(mat6));
%         rho_bladeMass_to_mat7 = R1_mat7(1, 2:6) ./ (std(bladeMass) .* std(mat7));
%         rho_bladeMass_to_mat8 = R1_mat8(1, 2:6) ./ (std(bladeMass) .* std(mat8));
% 
%         Y1 = [rho_bladeMass_to_mat1; ...
%               rho_bladeMass_to_mat2; ...
%               rho_bladeMass_to_mat3; ...
%               rho_bladeMass_to_mat4; ...
%               rho_bladeMass_to_mat5; ...
%               rho_bladeMass_to_mat6; ...
%               rho_bladeMass_to_mat7; ...
%               rho_bladeMass_to_mat8];
% 
%         bar(Y1,'grouped');
%         ylabel('blade mass, \rho_{m_{blade}}')
%         ylim([-1 1]);
%         %     legend('E_{11}','E_{22}','G_{12}','\nu_{12}','\rho','Location','EastOutside')
% 
%         % Save label positions 
%         xtl = get(gca,'xticklabel'); 
%         xt  = get(gca,'xtick');
%         ytl = get(gca,'yticklabel'); 
%         yt  = get(gca,'ytick'); 
%         xl  = get(gca,'xlabel'); 
%         yl  = get(gca,'ylabel'); 
%         orig_label_units = get( xl, 'Units'); 
%         set([xl yl], 'Units','pixels'); 
%         xp = get( xl, 'Position'); 
%         yp = get( yl, 'Position'); 
%         % Generate the text object ticks (changing labels' Position): 
%         [hx,hy] = format_ticks(gca,{'blade-root','blade-shell','spar-uni','spar-core','LEP-core','TEP-core','web-shell','web-core'},cellstr(ytl),xt,yt,45,0,0);
%         % Restore original label Position and Units properties: 
%         set( xl, 'Position', xp); 
%         set( yl, 'Position', yp); 
%         set([xl yl], 'Units', orig_label_units);

%% =========================================================
% tipDeflect vs. {E_11 for all mats}, {E_22 for all mats}, ... and {density for all mats}
% 
%         figure
%         set(gcf, 'color', 'white')
% 
%         R2_mat1 = cov( [tipDeflect, mat1] );
%         R2_mat2 = cov( [tipDeflect, mat2] );
%         R2_mat3 = cov( [tipDeflect, mat3] );
%         R2_mat4 = cov( [tipDeflect, mat4] );
%         R2_mat5 = cov( [tipDeflect, mat5] );
%         R2_mat6 = cov( [tipDeflect, mat6] );
%         R2_mat7 = cov( [tipDeflect, mat7] );
%         R2_mat8 = cov( [tipDeflect, mat8] );
% 
%         rho_tipDeflect_to_mat1 = R2_mat1(1, 2:6) ./ (std(tipDeflect) .* std(mat1));
%         rho_tipDeflect_to_mat2 = R2_mat2(1, 2:6) ./ (std(tipDeflect) .* std(mat2));
%         rho_tipDeflect_to_mat3 = R2_mat3(1, 2:6) ./ (std(tipDeflect) .* std(mat3));
%         rho_tipDeflect_to_mat4 = R2_mat4(1, 2:6) ./ (std(tipDeflect) .* std(mat4));
%         rho_tipDeflect_to_mat5 = R2_mat5(1, 2:6) ./ (std(tipDeflect) .* std(mat5));
%         rho_tipDeflect_to_mat6 = R2_mat6(1, 2:6) ./ (std(tipDeflect) .* std(mat6));
%         rho_tipDeflect_to_mat7 = R2_mat7(1, 2:6) ./ (std(tipDeflect) .* std(mat7));
%         rho_tipDeflect_to_mat8 = R2_mat8(1, 2:6) ./ (std(tipDeflect) .* std(mat8));
% 
%         Y2 = [rho_tipDeflect_to_mat1; ...
%               rho_tipDeflect_to_mat2; ...
%               rho_tipDeflect_to_mat3; ...
%               rho_tipDeflect_to_mat4; ...
%               rho_tipDeflect_to_mat5; ...
%               rho_tipDeflect_to_mat6; ...
%               rho_tipDeflect_to_mat7; ...
%               rho_tipDeflect_to_mat8];
% 
%         bar(Y2,'grouped');
%         ylabel('tip deflection, \rho_{\delta_{tip}}')
%         ylim([-1 1]);
%         %     legend('E_{11}','E_{22}','G_{12}','\nu_{12}','\rho','Location','EastOutside')
% 
%         % Save label positions 
%         xtl = get(gca,'xticklabel'); 
%         xt  = get(gca,'xtick');
%         ytl = get(gca,'yticklabel'); 
%         yt  = get(gca,'ytick'); 
%         xl = get(gca,'xlabel'); 
%         yl = get(gca,'ylabel'); 
%         orig_label_units = get( xl, 'Units'); 
%         set([xl yl], 'Units','pixels'); 
%         xp = get( xl, 'Position'); 
%         yp = get( yl, 'Position'); 
%         % Generate the text object ticks (changing labels' Position): 
%         [hx,hy] = format_ticks(gca,{'blade-root','blade-shell','spar-uni','spar-core','LEP-core','TEP-core','web-shell','web-core'},cellstr(ytl),xt,yt,45,0,0);
%         % Restore original label Position and Units properties: 
%         set( xl, 'Position', xp); 
%         set( yl, 'Position', yp); 
%         set([xl yl], 'Units', orig_label_units);

%% =========================================================
% maxBuckle in LEP vs. {E_11 for all mats}, {E_22 for all mats}, ... and {density for all mats} 
figure
set(gcf, 'color', 'white') 

var = [maxBuckleCrit.LEP]';

R3_mat1 = cov( [var, mat1] );
R3_mat2 = cov( [var, mat2] );
R3_mat3 = cov( [var, mat3] );
R3_mat4 = cov( [var, mat4] );
R3_mat5 = cov( [var, mat5] );
R3_mat6 = cov( [var, mat6] );
R3_mat7 = cov( [var, mat7] );
R3_mat8 = cov( [var, mat8] );

rho_maxBuckleCrit_to_mat1 = R3_mat1(1, 2:6) ./ (std(var) .* std(mat1));
rho_maxBuckleCrit_to_mat2 = R3_mat2(1, 2:6) ./ (std(var) .* std(mat2));
rho_maxBuckleCrit_to_mat3 = R3_mat3(1, 2:6) ./ (std(var) .* std(mat3));
rho_maxBuckleCrit_to_mat4 = R3_mat4(1, 2:6) ./ (std(var) .* std(mat4));
rho_maxBuckleCrit_to_mat5 = R3_mat5(1, 2:6) ./ (std(var) .* std(mat5));
rho_maxBuckleCrit_to_mat6 = R3_mat6(1, 2:6) ./ (std(var) .* std(mat6));
rho_maxBuckleCrit_to_mat7 = R3_mat7(1, 2:6) ./ (std(var) .* std(mat7));
rho_maxBuckleCrit_to_mat8 = R3_mat8(1, 2:6) ./ (std(var) .* std(mat8));

Y3 = [rho_maxBuckleCrit_to_mat1; ...
      rho_maxBuckleCrit_to_mat2; ...
      rho_maxBuckleCrit_to_mat3; ...
      rho_maxBuckleCrit_to_mat4; ...
      rho_maxBuckleCrit_to_mat5; ...
      rho_maxBuckleCrit_to_mat6; ...
      rho_maxBuckleCrit_to_mat7; ...
      rho_maxBuckleCrit_to_mat8];


bar(Y3,'grouped');
ylabel('max buckling criteria in LEP, \rho_{R_{b,LEP}}')
ylim([-1 1])
%     legend('E_{11}','E_{22}','G_{12}','\nu_{12}','\rho','Location','EastOutside')

% Save label positions 
xtl = get(gca,'xticklabel'); 
xt  = get(gca,'xtick');
ytl = get(gca,'yticklabel'); 
yt  = get(gca,'ytick'); 
xl = get(gca,'xlabel'); 
yl = get(gca,'ylabel'); 
orig_label_units = get( xl, 'Units'); 
set([xl yl], 'Units','pixels'); 
xp = get( xl, 'Position'); 
yp = get( yl, 'Position'); 
% Generate the text object ticks (changing labels' Position): 
[hx,hy] = format_ticks(gca,{'blade-root','blade-shell','spar-uni','spar-core','LEP-core','TEP-core','web-shell','web-core'},cellstr(ytl),xt,yt,45,0,0);
% Restore original label Position and Units properties: 
set( xl, 'Position', xp); 
set( yl, 'Position', yp); 
set([xl yl], 'Units', orig_label_units);

%% =========================================================
% maxBuckle in TEP vs. {E_11 for all mats}, {E_22 for all mats}, ... and {density for all mats} 
figure
set(gcf, 'color', 'white') 

var = [maxBuckleCrit.TEP]';

R3_mat1 = cov( [var, mat1] );
R3_mat2 = cov( [var, mat2] );
R3_mat3 = cov( [var, mat3] );
R3_mat4 = cov( [var, mat4] );
R3_mat5 = cov( [var, mat5] );
R3_mat6 = cov( [var, mat6] );
R3_mat7 = cov( [var, mat7] );
R3_mat8 = cov( [var, mat8] );

rho_maxBuckleCrit_to_mat1 = R3_mat1(1, 2:6) ./ (std(var) .* std(mat1));
rho_maxBuckleCrit_to_mat2 = R3_mat2(1, 2:6) ./ (std(var) .* std(mat2));
rho_maxBuckleCrit_to_mat3 = R3_mat3(1, 2:6) ./ (std(var) .* std(mat3));
rho_maxBuckleCrit_to_mat4 = R3_mat4(1, 2:6) ./ (std(var) .* std(mat4));
rho_maxBuckleCrit_to_mat5 = R3_mat5(1, 2:6) ./ (std(var) .* std(mat5));
rho_maxBuckleCrit_to_mat6 = R3_mat6(1, 2:6) ./ (std(var) .* std(mat6));
rho_maxBuckleCrit_to_mat7 = R3_mat7(1, 2:6) ./ (std(var) .* std(mat7));
rho_maxBuckleCrit_to_mat8 = R3_mat8(1, 2:6) ./ (std(var) .* std(mat8));

Y3 = [rho_maxBuckleCrit_to_mat1; ...
      rho_maxBuckleCrit_to_mat2; ...
      rho_maxBuckleCrit_to_mat3; ...
      rho_maxBuckleCrit_to_mat4; ...
      rho_maxBuckleCrit_to_mat5; ...
      rho_maxBuckleCrit_to_mat6; ...
      rho_maxBuckleCrit_to_mat7; ...
      rho_maxBuckleCrit_to_mat8];


bar(Y3,'grouped');
ylabel('max buckling criteria in TEP, \rho_{R_{b,TEP}}')
ylim([-1 1])
%     legend('E_{11}','E_{22}','G_{12}','\nu_{12}','\rho','Location','EastOutside')

% Save label positions 
xtl = get(gca,'xticklabel'); 
xt  = get(gca,'xtick');
ytl = get(gca,'yticklabel'); 
yt  = get(gca,'ytick'); 
xl = get(gca,'xlabel'); 
yl = get(gca,'ylabel'); 
orig_label_units = get( xl, 'Units'); 
set([xl yl], 'Units','pixels'); 
xp = get( xl, 'Position'); 
yp = get( yl, 'Position'); 
% Generate the text object ticks (changing labels' Position): 
[hx,hy] = format_ticks(gca,{'blade-root','blade-shell','spar-uni','spar-core','LEP-core','TEP-core','web-shell','web-core'},cellstr(ytl),xt,yt,45,0,0);
% Restore original label Position and Units properties: 
set( xl, 'Position', xp); 
set( yl, 'Position', yp); 
set([xl yl], 'Units', orig_label_units);

%% =========================================================
% maxBuckle in WEBS vs. {E_11 for all mats}, {E_22 for all mats}, ... and {density for all mats} 
%         figure
%         set(gcf, 'color', 'white') 
% 
%         var = [maxBuckleCrit.webs]';
% 
%         R4_mat1 = cov( [var, mat1] );
%         R4_mat2 = cov( [var, mat2] );
%         R4_mat3 = cov( [var, mat3] );
%         R4_mat4 = cov( [var, mat4] );
%         R4_mat5 = cov( [var, mat5] );
%         R4_mat6 = cov( [var, mat6] );
%         R4_mat7 = cov( [var, mat7] );
%         R4_mat8 = cov( [var, mat8] );
% 
%         rho_maxBuckleCrit_to_mat1 = R4_mat1(1, 2:6) ./ (std(var) .* std(mat1));
%         rho_maxBuckleCrit_to_mat2 = R4_mat2(1, 2:6) ./ (std(var) .* std(mat2));
%         rho_maxBuckleCrit_to_mat3 = R4_mat3(1, 2:6) ./ (std(var) .* std(mat3));
%         rho_maxBuckleCrit_to_mat4 = R4_mat4(1, 2:6) ./ (std(var) .* std(mat4));
%         rho_maxBuckleCrit_to_mat5 = R4_mat5(1, 2:6) ./ (std(var) .* std(mat5));
%         rho_maxBuckleCrit_to_mat6 = R4_mat6(1, 2:6) ./ (std(var) .* std(mat6));
%         rho_maxBuckleCrit_to_mat7 = R4_mat7(1, 2:6) ./ (std(var) .* std(mat7));
%         rho_maxBuckleCrit_to_mat8 = R4_mat8(1, 2:6) ./ (std(var) .* std(mat8));
% 
%         Y4 = [rho_maxBuckleCrit_to_mat1; ...
%               rho_maxBuckleCrit_to_mat2; ...
%               rho_maxBuckleCrit_to_mat3; ...
%               rho_maxBuckleCrit_to_mat4; ...
%               rho_maxBuckleCrit_to_mat5; ...
%               rho_maxBuckleCrit_to_mat6; ...
%               rho_maxBuckleCrit_to_mat7; ...
%               rho_maxBuckleCrit_to_mat8];
% 
% 
%         bar(Y4,'grouped');
%         ylabel('max buckling criteria in webs, \rho_{R_{b,webs}}')
%         ylim([-1 1])
%         %     legend('E_{11}','E_{22}','G_{12}','\nu_{12}','\rho','Location','EastOutside')
% 
%         % Save label positions 
%         xtl = get(gca,'xticklabel'); 
%         xt  = get(gca,'xtick');
%         ytl = get(gca,'yticklabel'); 
%         yt  = get(gca,'ytick'); 
%         xl = get(gca,'xlabel'); 
%         yl = get(gca,'ylabel'); 
%         orig_label_units = get( xl, 'Units'); 
%         set([xl yl], 'Units','pixels'); 
%         xp = get( xl, 'Position'); 
%         yp = get( yl, 'Position'); 
%         % Generate the text object ticks (changing labels' Position): 
%         [hx,hy] = format_ticks(gca,{'blade-root','blade-shell','spar-uni','spar-core','LEP-core','TEP-core','web-shell','web-core'},cellstr(ytl),xt,yt,45,0,0);
%         % Restore original label Position and Units properties: 
%         set( xl, 'Position', xp); 
%         set( yl, 'Position', yp); 
%         set([xl yl], 'Units', orig_label_units);

%% =========================================================
figure
set(gcf, 'color', 'white') 

tmp = cell2mat(maxStressFC_mat3);
var = [tmp.s_11_fc_C]';

R5_mat1 = cov( [var, mat1] );
R5_mat2 = cov( [var, mat2] );
R5_mat3 = cov( [var, mat3] );
R5_mat4 = cov( [var, mat4] );
R5_mat5 = cov( [var, mat5] );
R5_mat6 = cov( [var, mat6] );
R5_mat7 = cov( [var, mat7] );
R5_mat8 = cov( [var, mat8] );

rho_maxStress_to_mat1 = R5_mat1(1, 2:6) ./ (std(var) .* std(mat1));
rho_maxStress_to_mat2 = R5_mat2(1, 2:6) ./ (std(var) .* std(mat2));
rho_maxStress_to_mat3 = R5_mat3(1, 2:6) ./ (std(var) .* std(mat3));
rho_maxStress_to_mat4 = R5_mat4(1, 2:6) ./ (std(var) .* std(mat4));
rho_maxStress_to_mat5 = R5_mat5(1, 2:6) ./ (std(var) .* std(mat5));
rho_maxStress_to_mat6 = R5_mat6(1, 2:6) ./ (std(var) .* std(mat6));
rho_maxStress_to_mat7 = R5_mat7(1, 2:6) ./ (std(var) .* std(mat7));
rho_maxStress_to_mat8 = R5_mat8(1, 2:6) ./ (std(var) .* std(mat8));

Y5 = [rho_maxStress_to_mat1; ...
      rho_maxStress_to_mat2; ...
      rho_maxStress_to_mat3; ...
      rho_maxStress_to_mat4; ...
      rho_maxStress_to_mat5; ...
      rho_maxStress_to_mat6; ...
      rho_maxStress_to_mat7; ...
      rho_maxStress_to_mat8];


bar(Y5,'grouped');
ylabel('correlation of max compressive stress, \rho_{\sigma_{11}}')
title('MATERIAL: "NCT307-D1-34-600" unidirectional carbon fiber COMPONENT: spar cap')
ylim([-1 1])
%     legend('E_{11}','E_{22}','G_{12}','\nu_{12}','\rho','Location','EastOutside')

% Save label positions 
xtl = get(gca,'xticklabel'); 
xt  = get(gca,'xtick');
ytl = get(gca,'yticklabel'); 
yt  = get(gca,'ytick'); 
xl = get(gca,'xlabel'); 
yl = get(gca,'ylabel'); 
orig_label_units = get( xl, 'Units'); 
set([xl yl], 'Units','pixels'); 
xp = get( xl, 'Position'); 
yp = get( yl, 'Position'); 
% Generate the text object ticks (changing labels' Position): 
[hx,hy] = format_ticks(gca,{'blade-root','blade-shell','spar-uni','spar-core','LEP-core','TEP-core','web-shell','web-core'},cellstr(ytl),xt,yt,45,0,0);
% Restore original label Position and Units properties: 
set( xl, 'Position', xp); 
set( yl, 'Position', yp); 
set([xl yl], 'Units', orig_label_units);

%% =========================================================
figure
set(gcf, 'color', 'white') 

tmp = cell2mat(maxStressFC_mat2);
var = [tmp.s_11_fc_C]';

R6_mat1 = cov( [var, mat1] );
R6_mat2 = cov( [var, mat2] );
R6_mat3 = cov( [var, mat3] );
R6_mat4 = cov( [var, mat4] );
R6_mat5 = cov( [var, mat5] );
R6_mat6 = cov( [var, mat6] );
R6_mat7 = cov( [var, mat7] );
R6_mat8 = cov( [var, mat8] );

rho_maxStress_to_mat1 = R6_mat1(1, 2:6) ./ (std(var) .* std(mat1));
rho_maxStress_to_mat2 = R6_mat2(1, 2:6) ./ (std(var) .* std(mat2));
rho_maxStress_to_mat3 = R6_mat3(1, 2:6) ./ (std(var) .* std(mat3));
rho_maxStress_to_mat4 = R6_mat4(1, 2:6) ./ (std(var) .* std(mat4));
rho_maxStress_to_mat5 = R6_mat5(1, 2:6) ./ (std(var) .* std(mat5));
rho_maxStress_to_mat6 = R6_mat6(1, 2:6) ./ (std(var) .* std(mat6));
rho_maxStress_to_mat7 = R6_mat7(1, 2:6) ./ (std(var) .* std(mat7));
rho_maxStress_to_mat8 = R6_mat8(1, 2:6) ./ (std(var) .* std(mat8));

Y6 = [rho_maxStress_to_mat1; ...
      rho_maxStress_to_mat2; ...
      rho_maxStress_to_mat3; ...
      rho_maxStress_to_mat4; ...
      rho_maxStress_to_mat5; ...
      rho_maxStress_to_mat6; ...
      rho_maxStress_to_mat7; ...
      rho_maxStress_to_mat8];


bar(Y6,'grouped');
title('MATERIAL: "NB307-D1-7781-497A" double-bias E-glass       COMPONENT: blade-shell')
ylabel('correlation of max compressive stress, \rho_{\sigma_{11}}')
ylim([-1 1])
%     legend('E_{11}','E_{22}','G_{12}','\nu_{12}','\rho','Location','EastOutside')

% Save label positions 
xtl = get(gca,'xticklabel'); 
xt  = get(gca,'xtick');
ytl = get(gca,'yticklabel'); 
yt  = get(gca,'ytick'); 
xl = get(gca,'xlabel'); 
yl = get(gca,'ylabel'); 
orig_label_units = get( xl, 'Units'); 
set([xl yl], 'Units','pixels'); 
xp = get( xl, 'Position'); 
yp = get( yl, 'Position'); 
% Generate the text object ticks (changing labels' Position): 
[hx,hy] = format_ticks(gca,{'blade-root','blade-shell','spar-uni','spar-core','LEP-core','TEP-core','web-shell','web-core'},cellstr(ytl),xt,yt,45,0,0);
% Restore original label Position and Units properties: 
set( xl, 'Position', xp); 
set( yl, 'Position', yp); 
set([xl yl], 'Units', orig_label_units);

%% =========================================================
%         figure
%         set(gcf, 'color', 'white') 
% 
%         tmp = cell2mat(maxStressFC_mat2);
%         var = [tmp.s_12_fc_S]';
% 
%         R7_mat1 = cov( [var, mat1] );
%         R7_mat2 = cov( [var, mat2] );
%         R7_mat3 = cov( [var, mat3] );
%         R7_mat4 = cov( [var, mat4] );
%         R7_mat5 = cov( [var, mat5] );
%         R7_mat6 = cov( [var, mat6] );
%         R7_mat7 = cov( [var, mat7] );
%         R7_mat8 = cov( [var, mat8] );
% 
%         rho_maxStress_to_mat1 = R7_mat1(1, 2:6) ./ (std(var) .* std(mat1));
%         rho_maxStress_to_mat2 = R7_mat2(1, 2:6) ./ (std(var) .* std(mat2));
%         rho_maxStress_to_mat3 = R7_mat3(1, 2:6) ./ (std(var) .* std(mat3));
%         rho_maxStress_to_mat4 = R7_mat4(1, 2:6) ./ (std(var) .* std(mat4));
%         rho_maxStress_to_mat5 = R7_mat5(1, 2:6) ./ (std(var) .* std(mat5));
%         rho_maxStress_to_mat6 = R7_mat6(1, 2:6) ./ (std(var) .* std(mat6));
%         rho_maxStress_to_mat7 = R7_mat7(1, 2:6) ./ (std(var) .* std(mat7));
%         rho_maxStress_to_mat8 = R7_mat8(1, 2:6) ./ (std(var) .* std(mat8));
% 
%         Y7 = [rho_maxStress_to_mat1; ...
%               rho_maxStress_to_mat2; ...
%               rho_maxStress_to_mat3; ...
%               rho_maxStress_to_mat4; ...
%               rho_maxStress_to_mat5; ...
%               rho_maxStress_to_mat6; ...
%               rho_maxStress_to_mat7; ...
%               rho_maxStress_to_mat8];
% 
% 
%         bar(Y7,'grouped');
%         ylabel('max stress in material #2, \rho_{\tau_{12,mat2}}')
%         ylim([-1 1])
%         %     legend('E_{11}','E_{22}','G_{12}','\nu_{12}','\rho','Location','EastOutside')
% 
%         % Save label positions 
%         xtl = get(gca,'xticklabel'); 
%         xt  = get(gca,'xtick');
%         ytl = get(gca,'yticklabel'); 
%         yt  = get(gca,'ytick'); 
%         xl = get(gca,'xlabel'); 
%         yl = get(gca,'ylabel'); 
%         orig_label_units = get( xl, 'Units'); 
%         set([xl yl], 'Units','pixels'); 
%         xp = get( xl, 'Position'); 
%         yp = get( yl, 'Position'); 
%         % Generate the text object ticks (changing labels' Position): 
%         [hx,hy] = format_ticks(gca,{'blade-root','blade-shell','spar-uni','spar-core','LEP-core','TEP-core','web-shell','web-core'},cellstr(ytl),xt,yt,45,0,0);
%         % Restore original label Position and Units properties: 
%         set( xl, 'Position', xp); 
%         set( yl, 'Position', yp); 
%         set([xl yl], 'Units', orig_label_units);

%% =========================================================
%         figure
%         set(gcf, 'color', 'white') 
% 
%         tmp = cell2mat(maxStressFC_mat7);
%         var = [tmp.s_12_fc_S]';
% 
%         R8_mat1 = cov( [var, mat1] );
%         R8_mat2 = cov( [var, mat2] );
%         R8_mat3 = cov( [var, mat3] );
%         R8_mat4 = cov( [var, mat4] );
%         R8_mat5 = cov( [var, mat5] );
%         R8_mat6 = cov( [var, mat6] );
%         R8_mat7 = cov( [var, mat7] );
%         R8_mat8 = cov( [var, mat8] );
% 
%         rho_maxStress_to_mat1 = R8_mat1(1, 2:6) ./ (std(var) .* std(mat1));
%         rho_maxStress_to_mat2 = R8_mat2(1, 2:6) ./ (std(var) .* std(mat2));
%         rho_maxStress_to_mat3 = R8_mat3(1, 2:6) ./ (std(var) .* std(mat3));
%         rho_maxStress_to_mat4 = R8_mat4(1, 2:6) ./ (std(var) .* std(mat4));
%         rho_maxStress_to_mat5 = R8_mat5(1, 2:6) ./ (std(var) .* std(mat5));
%         rho_maxStress_to_mat6 = R8_mat6(1, 2:6) ./ (std(var) .* std(mat6));
%         rho_maxStress_to_mat7 = R8_mat7(1, 2:6) ./ (std(var) .* std(mat7));
%         rho_maxStress_to_mat8 = R8_mat8(1, 2:6) ./ (std(var) .* std(mat8));
% 
%         Y8 = [rho_maxStress_to_mat1; ...
%               rho_maxStress_to_mat2; ...
%               rho_maxStress_to_mat3; ...
%               rho_maxStress_to_mat4; ...
%               rho_maxStress_to_mat5; ...
%               rho_maxStress_to_mat6; ...
%               rho_maxStress_to_mat7; ...
%               rho_maxStress_to_mat8];
% 
% 
%         bar(Y8,'grouped');
%         ylabel('max stress in material #7, \rho_{\tau_{12,mat7}}')
%         ylim([-1 1])
%         legend('E_{11}','E_{22}','G_{12}','\nu_{12}','\rho','Location','EastOutside')
% 
%         % Save label positions 
%         xtl = get(gca,'xticklabel'); 
%         xt  = get(gca,'xtick');
%         ytl = get(gca,'yticklabel'); 
%         yt  = get(gca,'ytick'); 
%         xl = get(gca,'xlabel'); 
%         yl = get(gca,'ylabel'); 
%         orig_label_units = get( xl, 'Units'); 
%         set([xl yl], 'Units','pixels'); 
%         xp = get( xl, 'Position'); 
%         yp = get( yl, 'Position'); 
%         % Generate the text object ticks (changing labels' Position): 
%         [hx,hy] = format_ticks(gca,{'blade-root','blade-shell','spar-uni','spar-core','LEP-core','TEP-core','web-shell','web-core'},cellstr(ytl),xt,yt,45,0,0);
%         % Restore original label Position and Units properties: 
%         set( xl, 'Position', xp); 
%         set( yl, 'Position', yp); 
%         set([xl yl], 'Units', orig_label_units);

%% =========================================================
%         figure
%         set(gcf, 'color', 'white') 
% 
%         % maximum s_11_mat3 vs. {bladeMass, tipDeflect, maxBuckle}
%         R4_mat1 = cov( [mat3_max_s_11, mat1] );
%         R4_mat2 = cov( [mat3_max_s_11, mat2] );
%         R4_mat3 = cov( [mat3_max_s_11, mat3] );
%         R4_mat4 = cov( [mat3_max_s_11, mat4] );
%         R4_mat5 = cov( [mat3_max_s_11, mat5] );
%         R4_mat6 = cov( [mat3_max_s_11, mat6] );
%         R4_mat7 = cov( [mat3_max_s_11, mat7] );
%         R4_mat8 = cov( [mat3_max_s_11, mat8] );
% 
%         rho_mat3_max_s_11_to_mat1 = R4_mat1(1, 2:6) ./ (mat3_max_s_11_std .* std(mat1));
%         rho_mat3_max_s_11_to_mat2 = R4_mat2(1, 2:6) ./ (mat3_max_s_11_std .* std(mat2));
%         rho_mat3_max_s_11_to_mat3 = R4_mat3(1, 2:6) ./ (mat3_max_s_11_std .* std(mat3));
%         rho_mat3_max_s_11_to_mat4 = R4_mat4(1, 2:6) ./ (mat3_max_s_11_std .* std(mat4));
%         rho_mat3_max_s_11_to_mat5 = R4_mat5(1, 2:6) ./ (mat3_max_s_11_std .* std(mat5));
%         rho_mat3_max_s_11_to_mat6 = R4_mat6(1, 2:6) ./ (mat3_max_s_11_std .* std(mat6));
%         rho_mat3_max_s_11_to_mat7 = R4_mat7(1, 2:6) ./ (mat3_max_s_11_std .* std(mat7));
%         rho_mat3_max_s_11_to_mat8 = R4_mat8(1, 2:6) ./ (mat3_max_s_11_std .* std(mat8));
% 
%         Y4 = [rho_mat3_max_s_11_to_mat1; ...
%               rho_mat3_max_s_11_to_mat2; ...
%               rho_mat3_max_s_11_to_mat3; ...
%               rho_mat3_max_s_11_to_mat4; ...
%               rho_mat3_max_s_11_to_mat5; ...
%               rho_mat3_max_s_11_to_mat6; ...
%               rho_mat3_max_s_11_to_mat7; ...
%               rho_mat3_max_s_11_to_mat8];
% 
% 
%         bar(Y4,'grouped');
%         ylabel('\rho_{max stress in material #3}')
%         legend('E_{11}','E_{22}','G_{12}','\nu_{12}','\rho','Location','EastOutside')
% 
%         % Save label positions 
%         xtl = get(gca,'xticklabel'); 
%         xt  = get(gca,'xtick');
%         ytl = get(gca,'yticklabel'); 
%         yt  = get(gca,'ytick'); 
%         xl = get(gca,'xlabel'); 
%         yl = get(gca,'ylabel'); 
%         orig_label_units = get( xl, 'Units'); 
%         set([xl yl], 'Units','pixels'); 
%         xp = get( xl, 'Position'); 
%         yp = get( yl, 'Position'); 
%         % Generate the text object ticks (changing labels' Position): 
%         [hx,hy] = format_ticks(gca,{'blade-root','blade-shell','spar-uni','spar-core','LEP-core','TEP-core','web-shell','web-core'},cellstr(ytl),xt,yt,90,0,0,0);
%         % Restore original label Position and Units properties: 
%         set( xl, 'Position', xp); 
%         set( yl, 'Position', yp); 
%         set([xl yl], 'Units', orig_label_units);
end   

