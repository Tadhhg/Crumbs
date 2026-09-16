function make_scenario_figures(S, M, o, cw)
%MAKE_SCENARIO_FIGURES  The two figures HEAT_SCENARIO produces.
%
%   figures/scenario_<name>_traces.png   drive, temperature vs time, depth
%                                        profiles, and a summary block
%   figures/scenario_<name>_maps.png     temperature maps through the run and
%                                        the thermal-dose map
%
%   The traces figure is the one to put a measurement next to: pass
%   'measured', [t T1 T2 ...] to HEAT_SCENARIO and the record is overlaid on
%   the predicted probe traces as dashed lines.

here = fileparts(mfilename('fullpath'));
fd = fullfile(fileparts(here),'figures');
if ~exist(fd,'dir'), mkdir(fd); end

t  = S.t(:).';          % fine grid: drive and traces
tp = S.t_map(:).';      % coarse grid: temperature maps and depth profiles
xm = S.x(:).'*1e3;  sm = -S.z(:).'*1e3;
Lant = 12;  Tb = M.T_blood;
np = size(S.probes,1);
lim = [43 50 60 100];

if o.save, vis = 'off'; else, vis = 'on'; end

% ===================================================== FIGURE 1: traces ===
f = figure('Visible',vis,'Position',[80 60 1180 760]);

% ---- drive ---------------------------------------------------------------
subplot(2,2,1);
stairs(t, S.P_drive, 'LineWidth',1.4); hold on;
yline(S.P_mean_generator,'r--','LineWidth',1.2, ...
      'Label',sprintf('mean %.2f W', S.P_mean_generator));
xline(S.t_stim,'k:','LineWidth',1.2,'Label','stimulation ends');
grid on; xlabel('time  (s)'); ylabel('generator forward power  (W)');
ylim([0 1.25*max(S.P_drive)+eps]);
if cw
    title(sprintf('Drive: %.3f W continuous for %.4g s', S.P_generator, S.t_stim));
else
    title(sprintf('Drive: %.3f W, %.0f %% duty at %.4g Hz, for %.4g s', ...
        S.P_generator, 100*S.duty, 1/S.period, S.t_stim));
end
xlim([0 t(end)]);

% ---- temperature vs time -------------------------------------------------
subplot(2,2,2);
h = gobjects(1,np+1);  lg = cell(1,np+1);
h(1) = plot(t, S.T_peak, 'k-', 'LineWidth',2); hold on;
lg{1} = 'peak in the block';
co = lines(max(np,1));
for r = 1:np
    h(r+1) = plot(t, S.T_probe(:,r), 'LineWidth',1.5, 'Color',co(r,:));
    lg{r+1} = sprintf('(%g, %g, %g) mm', 1e3*S.probes(r,1), 1e3*S.probes(r,2), ...
                      1e3*S.probes(r,3));
end
if ~isempty(o.measured)
    Mz = o.measured;
    for r = 1:min(np, size(Mz,2)-1)
        plot(Mz(:,1), Mz(:,r+1), '--', 'LineWidth',1.4, 'Color',co(r,:), ...
             'HandleVisibility','off');
    end
    text(0.02, 0.02, 'dashed = measured', 'Units','normalized', ...
         'FontAngle','italic');
end
for L = lim
    yline(L,':','Color',[.4 .4 .4],'HandleVisibility','off');
end
xline(S.t_stim,'k:','LineWidth',1.2,'HandleVisibility','off');
grid on; xlabel('time  (s)'); ylabel('temperature  (\circC)');
title('Temperature vs time  (dotted: 43, 50, 60, 100 \circC)');
legend(h, lg, 'Location','southeast','FontSize',8);
xlim([0 t(end)]);  ylim([Tb-1, max(S.T_max*1.08, 50)]);

% ---- depth profiles ------------------------------------------------------
subplot(2,2,3);
[~,ix6] = min(abs(S.x - 6e-3));
tsel = unique([S.t_stim*[0.1 0.3 1], tp(end)]);
tsel = tsel(tsel <= tp(end));
co2 = parula(numel(tsel)+1);  lg2 = cell(1,numel(tsel));
for q = 1:numel(tsel)
    [~,ii] = min(abs(tp - tsel(q)));
    plot(squeeze(S.T_plane(ix6,:,ii)), sm, 'LineWidth',1.7,'Color',co2(q,:));
    hold on;  lg2{q} = sprintf('t = %.3g s', tp(ii));
end
for L = lim, xline(L,':','Color',[.4 .4 .4],'HandleVisibility','off'); end
set(gca,'YDir','reverse'); grid on; ylim([0 30]);
xlabel('temperature  (\circC)'); ylabel('depth below the surface  (mm)');
title('Depth profile on the vertical x = 6 mm');
legend(lg2,'Location','southeast','FontSize',8);

% ---- summary block -------------------------------------------------------
subplot(2,2,4); axis off;
[~,ie] = min(abs(t - S.t_stim));
L = {};
L{end+1} = sprintf('\\bf%s\\rm', strrep(S.name,'_','\_'));
L{end+1} = '';
L{end+1} = sprintf('generator      %.3f W%s', S.P_generator, ...
    ternary(cw,' continuous', sprintf(', %.0f %% duty at %.4g Hz', 100*S.duty, 1/S.period)));
L{end+1} = sprintf('delivered      %.3f W while on  (\\eta = %.3f)', S.P_delivered, M.eta);
L{end+1} = sprintf('mean power     %.3f W generator', S.P_mean_generator);
L{end+1} = sprintf('stimulation    %.4g s     energy %.1f J', S.t_stim, S.energy_delivered);
L{end+1} = sprintf('surface        %s', ternary(strcmpi(o.model,'robin'), ...
    'convective, h = 1460 W/(m^2K)', 'held at 37 \circC (lower bound)'));
L{end+1} = '';
L{end+1} = sprintf('\\bfpeak %.2f \\circC\\rm at t = %.4g s', S.T_max, S.t_at_T_max);
L{end+1} = sprintf('   at (%.2f, %.2f, %.2f) mm', 1e3*S.hotspot);
L{end+1} = sprintf('peak at end of stimulation  %.2f \\circC', S.T_peak(ie));
L{end+1} = '';
for r = 1:np
    L{end+1} = sprintf('probe (%g, %g, %g) mm : %.2f \\circC at end, %.2f \\circC max', ...
        1e3*S.probes(r,1), 1e3*S.probes(r,2), 1e3*S.probes(r,3), ...
        S.T_probe(ie,r), max(S.T_probe(:,r))); %#ok<AGROW>
end
L{end+1} = '';
L{end+1} = sprintf('CEM43 \\geq 240 min : %s', zone_str(S.CEM43,240,S.x,S.z));
L{end+1} = sprintf('CEM43 \\geq  60 min : %s', zone_str(S.CEM43, 60,S.x,S.z));
L{end+1} = '';
if S.T_max > 100
    L{end+1} = '\bf*** above 100 \circC: the gel boils, the model does not\rm';
elseif S.T_max > 60
    L{end+1} = '* above 60 \circC: true rise is 10-15 % higher (Sec. 10K)';
elseif S.T_max > 50
    L{end+1} = '* above 50 \circC: check the Phytagel gel-sol transition';
else
    L{end+1} = 'nothing above 50 \circC -- constant properties are safe here';
end
text(0, 1, L, 'VerticalAlignment','top','FontName','FixedWidth','FontSize',8.5, ...
     'Interpreter','tex');

sgtitle(sprintf('Bench scenario -- %s surface model', o.model));
if o.save
    exportgraphics(f, fullfile(fd, sprintf('scenario_%s_traces.png', S.name)), ...
                   'Resolution',160);
    close(f);
end

% ======================================================= FIGURE 2: maps ===
f = figure('Visible',vis,'Position',[80 60 1250 640]);
tm = unique([S.t_stim*[0.05 0.2 0.5 1], tp(end)]);
tm = tm(tm <= tp(end));
tm = tm(1:min(5,numel(tm)));
cmax = max(max(S.T_max, 43.1), Tb+1);
% Blank the feed-exclusion sphere.  Inside it the field model puts a
% fictitious point charge (write-up Sec. 6.2), so the temperature there is not
% a prediction; leaving it on the map would dominate both the colour scale and
% the eye.  The blanked cells are drawn in grey and outlined.
msk = M.pmask_y0;
th = linspace(0,pi,120);  rex = 1e3*M.r_feed_excl;
cx = rex*cos(th);  cy = rex*sin(th);
for p = 1:numel(tm)
    [~,ii] = min(abs(tp - tm(p)));
    A = S.T_plane(:,:,ii);
    subplot(2,3,p);
    set(gca,'Color',[0.55 0.55 0.55]);
    imagesc(xm, sm, A.', 'AlphaData', double(msk.'));
    set(gca,'YDir','reverse','Color',[0.45 0.45 0.50]);
    clim([Tb cmax]); colormap(gca, hot); hold on;
    B = A;  B(~msk) = NaN;
    [c1,h1] = contour(xm, sm, B.', lim, 'w-','LineWidth',0.9);
    if ~isempty(c1), clabel(c1,h1,'Color','w','FontSize',7); end
    plot(cx, cy, 'w--', 'LineWidth',1);
    plot([0 Lant],[0 0],'c-','LineWidth',4);
    title(sprintf('t = %.3g s   (peak %.1f \\circC)', tp(ii), max(B(:))));
    if p >= 4, xlabel('x (mm)'); end
    if mod(p,3) == 1, ylabel('depth (mm)'); end
    if p == numel(tm)
        cb = colorbar; cb.Label.String = 'T (\circC)';
    end
    xlim([xm(1) xm(end)]); ylim([sm(1) sm(end)]);
end
subplot(2,3,6);
set(gca,'Color',[0.55 0.55 0.55]);
imagesc(xm, sm, log10(max(S.CEM43.',1e-3)), 'AlphaData', double(msk.'));
set(gca,'YDir','reverse','Color',[0.45 0.45 0.50]);
clim([-2 6]); colormap(gca, parula); cb = colorbar;
cb.Label.String = 'log_{10} CEM43  (min)';
hold on;
D = S.CEM43;  D(~msk) = NaN;
[c1,h1] = contour(xm, sm, D.', [60 240], 'w-','LineWidth',2);
if ~isempty(c1), clabel(c1,h1,'Color','w','FontSize',8); end
plot(cx, cy, 'w--', 'LineWidth',1);
plot([0 Lant],[0 0],'r-','LineWidth',4);
xlabel('x (mm)'); ylabel('depth (mm)');
xlim([xm(1) xm(end)]); ylim([sm(1) sm(end)]);
title('Thermal dose over the whole run');
sgtitle({sprintf(['Temperature on the y = 0 plane -- %s, %s surface ' ...
    '(cyan/red bar = the 12 mm monopole)'], strrep(S.name,'_','\_'), o.model), ...
    sprintf(['grey, inside the dashed %.0f mm arc: the feed singularity of the ' ...
    'field model -- not a prediction (Sec. 6.2)'], rex)});
if o.save
    exportgraphics(f, fullfile(fd, sprintf('scenario_%s_maps.png', S.name)), ...
                   'Resolution',160);
    close(f);
end
end

% =========================================================================
function s = zone_str(D, lev, x, z)
m = D >= lev;
if ~any(m(:)), s = 'none'; return; end
[ii,jj] = find(m);
s = sprintf('%.2f mm deep, %.1f mm long', max(-z(jj))*1e3, ...
            (max(x(ii))-min(x(ii)))*1e3);
end

function v = ternary(c,a,b), if c, v = a; else, v = b; end, end
