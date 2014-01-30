function plotMonteCarloStressFC(SIM, nBins, SAMPLE_SIZE, sample, title_str)

figTitle = [SIM.case ' Monte Carlo analysis'];
fig = figure('name',          figTitle, ...
             'color',         'white', ...
             'units',         'normalized',...
             'outerposition', [0.1 0.1 0.8 0.8]);
set(fig, 'color', 'white')
       
p1 = [sample.s_11_fc_T]';
p2 = [sample.s_11_fc_C]';
p3 = [sample.s_22_fc_T]';
p4 = [sample.s_22_fc_C]';
p5 = [sample.s_12_fc_S]';

%% compute basic statistics
mu1 = mean(p1);
mu2 = mean(p2);
mu3 = mean(p3);
mu4 = mean(p4);
mu5 = mean(p5);

sigma1 = std(p1);
sigma2 = std(p2);
sigma3 = std(p3);
sigma4 = std(p4);
sigma5 = std(p5);

cov1 = sigma1 / mu1;
cov2 = sigma2 / mu2;
cov3 = sigma3 / mu3;
cov4 = sigma4 / mu4;
cov5 = sigma5 / mu5;

%% compute the PDF (probability density function)
[count1, x1] = hist(p1, nBins);
[count2, x2] = hist(p2, nBins);
[count3, x3] = hist(p3, nBins);
[count4, x4] = hist(p4, nBins);
[count5, x5] = hist(p5, nBins);

bin_width1 = diff(x1(1:2));
bin_width2 = diff(x2(1:2));
bin_width3 = diff(x3(1:2));
bin_width4 = diff(x4(1:2));
bin_width5 = diff(x5(1:2));

pdf1 = count1 / SAMPLE_SIZE / bin_width1;
pdf2 = count2 / SAMPLE_SIZE / bin_width2;
pdf3 = count3 / SAMPLE_SIZE / bin_width3;
pdf4 = count4 / SAMPLE_SIZE / bin_width4;
pdf5 = count5 / SAMPLE_SIZE / bin_width5;

%% compute the normalized PDF
z1 = (x1 - mu1)/sigma1;
z2 = (x2 - mu2)/sigma2;
z3 = (x3 - mu3)/sigma3;
z4 = (x4 - mu4)/sigma4;
z5 = (x5 - mu5)/sigma5;

norm_pdf1 = sigma1 .* pdf1;
norm_pdf2 = sigma2 .* pdf2;
norm_pdf3 = sigma3 .* pdf3;
norm_pdf4 = sigma4 .* pdf4;
norm_pdf5 = sigma5 .* pdf5;

%% fit a probability distribution to the data
fit_p1 = fitdist(p1, 'normal');
fit_p2 = fitdist(p2, 'normal');
fit_p3 = fitdist(p3, 'normal');
fit_p4 = fitdist(p4, 'normal');
fit_p5 = fitdist(p5, 'normal');
% fit_p1 = fitdist(p1, 'InverseGaussian');
% fit_p2 = fitdist(p2, 'InverseGaussian');
% fit_p3 = fitdist(p3, 'InverseGaussian');
% fit_p4 = fitdist(p4, 'InverseGaussian');
% fit_p5 = fitdist(p5, 'InverseGaussian');

% idx1 = linspace(min(p1), max(p1), SAMPLE_SIZE);
% idx2 = linspace(min(p2), max(p2), SAMPLE_SIZE);
% idx3 = linspace(min(p3), max(p3), SAMPLE_SIZE);
% idx4 = linspace(min(p4), max(p4), SAMPLE_SIZE);
% idx5 = linspace(min(p5), max(p5), SAMPLE_SIZE);
xlim = 1.2;
idx1 = linspace(min(0, min(p1)), max(xlim, max(p1)), SAMPLE_SIZE);
idx2 = linspace(min(0, min(p2)), max(xlim, max(p2)), SAMPLE_SIZE);
idx3 = linspace(min(0, min(p3)), max(xlim, max(p3)), SAMPLE_SIZE);
idx4 = linspace(min(0, min(p4)), max(xlim, max(p4)), SAMPLE_SIZE);
idx5 = linspace(min(0, min(p5)), max(xlim, max(p5)), SAMPLE_SIZE);

% plot the normalized PDFs
hold on
plot(x1, norm_pdf1, '.r')
plot(x2, norm_pdf2, 'ob')
plot(x3, norm_pdf3, 'xk')
plot(x4, norm_pdf4, '.m')
plot(x5, norm_pdf5, '^c')

legend('p_1: \sigma_{11,max} / \sigma_{11,fT}', ...
       'p_2: \sigma_{11,min} / \sigma_{11,fC}', ...
       'p_3: \sigma_{22,max} / \sigma_{22,yT}', ...
       'p_4: \sigma_{22,min} / \sigma_{22,yC}', ...
       'p_5: |\tau_{12}|,_{max} / \tau_{12,y}')

plot(idx1, sigma1 .* pdf('Normal', idx1, fit_p1.mu, fit_p1.sigma), '-r');
plot(idx2, sigma2 .* pdf('Normal', idx2, fit_p2.mu, fit_p2.sigma), '-b');
plot(idx3, sigma3 .* pdf('Normal', idx3, fit_p3.mu, fit_p3.sigma), '-k');
plot(idx4, sigma4 .* pdf('Normal', idx4, fit_p4.mu, fit_p4.sigma), '-m');
plot(idx5, sigma5 .* pdf('Normal', idx5, fit_p5.mu, fit_p5.sigma), '-c');
% plot(idx1, sigma1 .* pdf('inversegaussian', idx1, fit_p1.mu, fit_p1.lambda), '-r');
% plot(idx2, sigma2 .* pdf('inversegaussian', idx2, fit_p2.mu, fit_p2.lambda), '-b');
% plot(idx3, sigma3 .* pdf('inversegaussian', idx3, fit_p3.mu, fit_p3.lambda), '-k');
% plot(idx4, sigma4 .* pdf('inversegaussian', idx4, fit_p4.mu, fit_p4.lambda), '-m');
% plot(idx5, sigma5 .* pdf('inversegaussian', idx5, fit_p5.mu, fit_p5.lambda), '-c');

title(title_str)
box on
xlabel('maximum stress penalty factor, p_i')
ylabel('normalized PDF')

% check for correctness, integration of the normalized PDF should equal 1
% trapz(x1, pdf1)
% trapz(x2, pdf2)
% trapz(x3, pdf3)
% trapz(x4, pdf4)
% trapz(x5, pdf5)
% 
% trapz(z1, norm_pdf1)
% trapz(z2, norm_pdf2)
% trapz(z3, norm_pdf3)
% trapz(z4, norm_pdf4)
% trapz(z5, norm_pdf5)

%% calculate probability of failure
pFail1 = 100 * trapz(x1(x1>1), pdf1(x1>1));
pFail2 = 100 * trapz(x2(x2>1), pdf2(x2>1));
pFail3 = 100 * trapz(x3(x3>1), pdf3(x3>1));
pFail4 = 100 * trapz(x4(x4>1), pdf4(x4>1));
pFail5 = 100 * trapz(x5(x5>1), pdf5(x5>1));

%% tabulate the data
disp(title_str)
fprintf(1, 'p_i     mean    std    cov    P_fail \n');
fprintf(1, '%s    %3.2f     %3.2f     %3.2f     %5.1f     \n', 'p1', mu1, sigma1, cov1, pFail1);
fprintf(1, '%s    %3.2f     %3.2f     %3.2f     %5.1f     \n', 'p2', mu2, sigma2, cov2, pFail2);
fprintf(1, '%s    %3.2f     %3.2f     %3.2f     %5.1f     \n', 'p3', mu3, sigma3, cov3, pFail3);
fprintf(1, '%s    %3.2f     %3.2f     %3.2f     %5.1f     \n', 'p4', mu4, sigma4, cov4, pFail4);
fprintf(1, '%s    %3.2f     %3.2f     %3.2f     %5.1f     \n', 'p5', mu5, sigma5, cov5, pFail5);

end

