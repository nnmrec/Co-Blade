function plotMonteCarloBuckleFC(iSIM, SIM, nBins, SAMPLE_SIZE, sample, title_str)

figTitle = [SIM.case{iSIM} ' Monte Carlo analysis'];
fig = figure('name',          figTitle, ...
             'color',         'white', ...
             'units',         'normalized',...
             'outerposition', [0.1 0.1 0.8 0.8]);
set(fig, 'color', 'white')
       
p1 = [sample.root]';
p2 = [sample.spar]';
p3 = [sample.LEP]';
p4 = [sample.TEP]';
p5 = [sample.webs]';
p6 = [sample.tip]';

%% compute basic statistics
mu1 = mean(p1);
mu2 = mean(p2);
mu3 = mean(p3);
mu4 = mean(p4);
mu5 = mean(p5);
mu6 = mean(p6);

sigma1 = std(p1);
sigma2 = std(p2);
sigma3 = std(p3);
sigma4 = std(p4);
sigma5 = std(p5);
sigma6 = std(p6);

cov1 = sigma1 / mu1;
cov2 = sigma2 / mu2;
cov3 = sigma3 / mu3;
cov4 = sigma4 / mu4;
cov5 = sigma5 / mu5;
cov6 = sigma6 / mu6;

%% compute the PDF (probability density function)
[count1, x1] = hist(p1, nBins);
[count2, x2] = hist(p2, nBins);
[count3, x3] = hist(p3, nBins);
[count4, x4] = hist(p4, nBins);
[count5, x5] = hist(p5, nBins);
[count6, x6] = hist(p6, nBins);

bin_width1 = diff(x1(1:2));
bin_width2 = diff(x2(1:2));
bin_width3 = diff(x3(1:2));
bin_width4 = diff(x4(1:2));
bin_width5 = diff(x5(1:2));
bin_width6 = diff(x6(1:2));

pdf1 = count1 / SAMPLE_SIZE / bin_width1;
pdf2 = count2 / SAMPLE_SIZE / bin_width2;
pdf3 = count3 / SAMPLE_SIZE / bin_width3;
pdf4 = count4 / SAMPLE_SIZE / bin_width4;
pdf5 = count5 / SAMPLE_SIZE / bin_width5;
pdf6 = count6 / SAMPLE_SIZE / bin_width6;

%% compute the normalized PDF
z1 = (x1 - mu1)/sigma1;
z2 = (x2 - mu2)/sigma2;
z3 = (x3 - mu3)/sigma3;
z4 = (x4 - mu4)/sigma4;
z5 = (x5 - mu5)/sigma5;
z6 = (x6 - mu6)/sigma6;

norm_pdf1 = sigma1 .* pdf1;
norm_pdf2 = sigma2 .* pdf2;
norm_pdf3 = sigma3 .* pdf3;
norm_pdf4 = sigma4 .* pdf4;
norm_pdf5 = sigma5 .* pdf5;
norm_pdf6 = sigma6 .* pdf6;

%% fit a probability distribution to the data
fit_p1 = fitdist(p1, 'normal');
fit_p2 = fitdist(p2, 'normal');
fit_p3 = fitdist(p3, 'normal');
fit_p4 = fitdist(p4, 'normal');
fit_p5 = fitdist(p5, 'normal');
fit_p6 = fitdist(p6, 'normal');

xlim = 1.2;
idx1 = linspace(min(0, min(p1)), max(xlim, max(p1)), SAMPLE_SIZE);
idx2 = linspace(min(0, min(p2)), max(xlim, max(p2)), SAMPLE_SIZE);
idx3 = linspace(min(0, min(p3)), max(xlim, max(p3)), SAMPLE_SIZE);
idx4 = linspace(min(0, min(p4)), max(xlim, max(p4)), SAMPLE_SIZE);
idx5 = linspace(min(0, min(p5)), max(xlim, max(p5)), SAMPLE_SIZE);
idx6 = linspace(min(0, min(p6)), max(xlim, max(p6)), SAMPLE_SIZE);

% plot the normalized PDFs
hold on
plot(x1, norm_pdf1, '.r')
plot(x2, norm_pdf2, 'ob')
plot(x3, norm_pdf3, 'xk')
plot(x4, norm_pdf4, '.m')
plot(x5, norm_pdf5, '^c')
plot(x6, norm_pdf6, '*g')

legend('root build-up', ...
       'spar cap', ...
       'LEP', ...
       'TEP', ...
       'webs', ...
       'blade tip')

plot(idx1, sigma1 .* pdf('Normal', idx1, fit_p1.mu, fit_p1.sigma), '-r');
plot(idx2, sigma2 .* pdf('Normal', idx2, fit_p2.mu, fit_p2.sigma), '-b');
plot(idx3, sigma3 .* pdf('Normal', idx3, fit_p3.mu, fit_p3.sigma), '-k');
plot(idx4, sigma4 .* pdf('Normal', idx4, fit_p4.mu, fit_p4.sigma), '-m');
plot(idx5, sigma5 .* pdf('Normal', idx5, fit_p5.mu, fit_p5.sigma), '-c');
plot(idx6, sigma6 .* pdf('Normal', idx6, fit_p6.mu, fit_p6.sigma), '-g');

title(title_str)
box on
xlabel('buckling penalty factor, p_6')
ylabel('normalized PDF')

%% calculate probability of failure
pFail1 = 100 * trapz(x1(x1>1), pdf1(x1>1));
pFail2 = 100 * trapz(x2(x2>1), pdf2(x2>1));
pFail3 = 100 * trapz(x3(x3>1), pdf3(x3>1));
pFail4 = 100 * trapz(x4(x4>1), pdf4(x4>1));
pFail5 = 100 * trapz(x5(x5>1), pdf5(x5>1));
pFail6 = 100 * trapz(x6(x6>1), pdf5(x6>1));

%% tabulate the data
disp(title_str)
fprintf(1, 'p_i     mean    std    cov    P_fail \n');
fprintf(1, '%s    %3.2f     %3.2f     %3.2f     %5.1f     \n', 'root build-up', mu1, sigma1, cov1, pFail1);
fprintf(1, '%s    %3.2f     %3.2f     %3.2f     %5.1f     \n', 'spar cap     ', mu2, sigma2, cov2, pFail2);
fprintf(1, '%s    %3.2f     %3.2f     %3.2f     %5.1f     \n', 'LEP          ', mu3, sigma3, cov3, pFail3);
fprintf(1, '%s    %3.2f     %3.2f     %3.2f     %5.1f     \n', 'TEP          ', mu4, sigma4, cov4, pFail4);
fprintf(1, '%s    %3.2f     %3.2f     %3.2f     %5.1f     \n', 'webs         ', mu5, sigma5, cov5, pFail5);
fprintf(1, '%s    %3.2f     %3.2f     %3.2f     %5.1f     \n', 'blade tip    ', mu6, sigma6, cov6, pFail6);

end

