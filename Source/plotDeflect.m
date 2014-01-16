function plotDeflect(SIM, BLADE, Disp, OUT)

figTitle = [SIM.case ' centroidal deflections'];
fig = figure('name', figTitle, ...
             'color', 'white', ...
             'units','normalized',...
             'outerposition',[0.1 0.1 0.8 0.8]);

subplot(2,1,1);
hold on
plot(BLADE.zSec, Disp.uo, 'o-b')
plot(BLADE.zSec, Disp.vo, 's-g')
plot(BLADE.zSec, Disp.wo, '.-r')
ylabel('deflection (m)')
xlabel('z (m)')
box on

legend('u_o', ...
       'v_o', ...
       'w_o', ...
       'Location', 'NorthWest')

subplot(2,1,2);
hold on
plot(BLADE.zSec, Disp.tz .* 180/pi, '.-k')
ylabel('twist angle (deg)')
xlabel('z (m)')
box on

if OUT.SAVE_PLOTS
    savePlots(fig, figTitle, OUT.SAVE_FIG_FMT)
end
    
end % function plotDeflect