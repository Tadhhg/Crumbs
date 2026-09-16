function make_sensitivity_figure(SS, root)
%MAKE_SENSITIVITY_FIGURE  figures/fig9_sensitivity.png from the study output.
%   Split out of RUN_HEAT_SENSITIVITY so the figure can be redrawn from
%   results/heat_sensitivity.mat without repeating the 22 solves:
%       load Heatmodel/results/heat_sensitivity.mat
%       make_sensitivity_figure(SS, 'Heatmodel')

f = figure('Visible','off','Position',[100 100 1250 640]);
subplot(1,2,1);
nm = {SS.name};  P = [SS.P100];
[~,ord] = sort(P);
barh(P(ord));
set(gca,'YTick',1:numel(P),'YTickLabel',nm(ord),'FontSize',7, ...
        'TickLabelInterpreter','none');
xlabel('generator power for T_{max} = 100 \circC after 60 s  (W)');
grid on; title('Sensitivity of the operating power');
xline(P(1),'r--','LineWidth',1.2);

subplot(1,2,2);
keep = [1 2 4 5 9 12 16 20];
keep = keep(keep <= numel(SS));
co = lines(numel(keep)); lg = cell(1,numel(keep));
for q = 1:numel(keep)
    c = keep(q);
    loglog(SS(c).t, SS(c).Th, 'LineWidth',1.7,'Color',co(q,:)); hold on;
    lg{q} = SS(c).name;
end
grid on; xlabel('time  (s)');
ylabel('\Theta_{max}  (K per 150 W delivered)');
title('Peak temperature rise vs time');
legend(lg,'Location','southeast','FontSize',7,'Interpreter','none');
exportgraphics(f, fullfile(root,'figures','fig9_sensitivity.png'),'Resolution',160);
close(f);
end
