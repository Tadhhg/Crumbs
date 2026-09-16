function make_heat_figures(R, root)
%MAKE_HEAT_FIGURES  The eight figures of the thermal model.
%   Called by RUN_HEAT_MODEL; R is its output struct, root the Heatmodel
%   folder.  Everything is derived from the single stored solution Theta by
%   the exact power scaling, so no figure needs a new solve.

fd = fullfile(root,'figures');
xm = R.x(:).'*1e3;  sm = -R.z(:).'*1e3;         % x [mm], depth [mm]
Lant = R.Tc.L_ant*1e3;
Tb = R.T_blood;  eta = R.eta;  Pref = R.Pdel_ref;
Th_y0 = double(R.Th_y0);                        % [nx x nz x Nt]
t = R.t;

% K per WATT OF GENERATOR POWER -- the power-independent core result
KperW_ss = double(R.Th_ss_y0)*eta/Pref;
[~,i10] = min(abs(t-10));
KperW_10 = Th_y0(:,:,i10)*eta/Pref;

% the design power: peak reaches 100 degC after 60 s
[~,i60] = min(abs(t-60));
Pdes = Pref*(100-Tb)/(eta*R.Th_max(i60));
fprintf('\nfigures: design power = %.3f W generator (%.3f W delivered)\n', ...
        Pdes, Pdes*eta);
sdes = Pdes*eta/Pref;

% ---- Fig 1: temperature rise per watt -----------------------------------
f = figure('Visible','off','Position',[100 100 1150 420]);
pan = {KperW_10, KperW_ss};
ttl = {'after 10 s', 'steady state (t \rightarrow \infty)'};
for p = 1:2
    subplot(1,2,p);
    V = log10(max(pan{p}.',1e-4));
    imagesc(xm, sm, V); set(gca,'YDir','reverse');
    clim([-3 log10(max(KperW_ss(:)))]);
    colormap(turbo); cb = colorbar;
    cb.Label.String = 'log_{10} ( K per watt of generator power )';
    hold on;
    contour(xm, sm, pan{p}.', [0.1 0.3 1 3 10], 'w-','LineWidth',0.7);
    plot([0 Lant],[0 0],'w-','LineWidth',5); plot([0 Lant],[0 0],'k-','LineWidth',2);
    xlabel('x  (mm)'); if p==1, ylabel('depth below the surface  (mm)'); end
    title(['\Theta / P_{gen},  ' ttl{p}]);
end
sgtitle(['Temperature rise per watt of generator power, y = 0 plane ' ...
         '(black bar = the 12 mm monopole, held at 37 \circC by the blood)']);
exportgraphics(f, fullfile(fd,'fig1_rise_per_watt.png'),'Resolution',160); close(f);

% ---- Fig 2: temperature maps at the design power ------------------------
f = figure('Visible','off','Position',[100 100 1100 640]);
tshow = [1 10 60 600];
for p = 1:4
    [~,ii] = min(abs(t-tshow(p)));
    Tmap = Tb + sdes*Th_y0(:,:,ii).';
    subplot(2,2,p);
    imagesc(xm, sm, Tmap); set(gca,'YDir','reverse'); clim([37 100]);
    colormap(hot); hold on;
    [c1,h1] = contour(xm, sm, Tmap, [43 50 60 80], 'w-','LineWidth',0.9);
    clabel(c1,h1,'Color','w','FontSize',7);
    plot([0 Lant],[0 0],'c-','LineWidth',4);
    title(sprintf('t = %.3g s   (T_{max} = %.1f \\circC)', t(ii), Tb+sdes*R.Th_max(ii)));
    if p>2, xlabel('x (mm)'); end
    if mod(p,2)==1, ylabel('depth (mm)'); end
end
cb = colorbar('Position',[0.93 0.11 0.015 0.8]); cb.Label.String = 'T (\circC)';
sgtitle(sprintf(['Gel temperature, y = 0 plane, %.2f W generator ' ...
    '(%.2f W delivered) -- the power that just reaches 100 \\circC at 60 s'], ...
    Pdes, Pdes*eta));
exportgraphics(f, fullfile(fd,'fig2_T_maps.png'),'Resolution',160); close(f);

% ---- Fig 3: peak temperature vs time, several powers --------------------
f = figure('Visible','off','Position',[100 100 820 520]);
Pl = [1 2 5 10 20 50 100 200];
co = turbo(numel(Pl)+2); lg = cell(1,numel(Pl));
for p = 1:numel(Pl)
    semilogx(t(2:end), Tb + Pl(p)*eta/Pref*R.Th_max(2:end), ...
             'LineWidth',1.7,'Color',co(p+1,:)); hold on;
    lg{p} = sprintf('%g W', Pl(p));
end
for L = [43 50 60 100]
    yline(L,'k--','LineWidth',0.9,'Label',sprintf('%g \\circC',L), ...
          'LabelHorizontalAlignment','left');
end
grid on; set(gca,'YScale','log'); ylim([37 1e4]);
xlabel('time  (s)'); ylabel('peak temperature in the block  (\circC)');
title('Peak gel temperature vs time, by generator forward power');
legend(lg,'Location','northwest','NumColumns',2);
exportgraphics(f, fullfile(fd,'fig3_Tmax_vs_time.png'),'Resolution',160); close(f);

% ---- Fig 4: the safe-operating chart ------------------------------------
f = figure('Visible','off','Position',[100 100 820 520]);
lims = [43 50 60 100];  co = lines(4); lg = cell(1,4); hL = gobjects(1,4);
for p = 1:4
    Pcurve = Pref*(lims(p)-Tb)./(eta*R.Th_max(2:end));
    hL(p) = loglog(t(2:end), Pcurve, 'LineWidth',1.9,'Color',co(p,:)); hold on;
    Pss = Pref*(lims(p)-Tb)/(eta*R.Th_max_ss);
    yline(Pss,'--','Color',co(p,:),'LineWidth',1,'HandleVisibility','off');
    lg{p} = sprintf('%g \\circC  (steady limit %.2f W)', lims(p), Pss);
end
grid on; xlabel('exposure time  (s)');
ylabel('generator forward power that just reaches the limit  (W)');
title({'Safe-operating chart: how much power, for how long', ...
       'dashed = the steady-state limit that can be run indefinitely'});
legend(hL, lg, 'Location','southwest');
xlim([1e-2 t(end)]);
exportgraphics(f, fullfile(fd,'fig4_safe_power.png'),'Resolution',160); close(f);

% ---- Fig 5: depth profiles and the sub-surface peak ---------------------
f = figure('Visible','off','Position',[100 100 1100 460]);
[~,ix6] = min(abs(R.x - 6e-3));
subplot(1,2,1);
tp = [0.1 1 10 60 600];  co = lines(numel(tp)); lg = cell(1,numel(tp)+1);
for p = 1:numel(tp)
    [~,ii] = min(abs(t-tp(p)));
    plot(Tb + sdes*squeeze(Th_y0(ix6,:,ii)), sm, 'LineWidth',1.7,'Color',co(p,:));
    hold on; lg{p} = sprintf('t = %g s', tp(p));
end
plot(Tb + sdes*double(R.Th_ss_y0(ix6,:)), sm, 'k--','LineWidth',1.7);
lg{end} = 'steady state';
set(gca,'YDir','reverse'); grid on; ylim([0 30]);
xlabel('T  (\circC)'); ylabel('depth below the surface  (mm)');
title('Depth profile on the vertical x = 6 mm');
legend(lg,'Location','southeast');

subplot(1,2,2);
prof = double(R.Th_ss_y0(ix6,:));
plot(sm, Tb + sdes*prof, 'LineWidth',2); hold on;
[pk,ipk] = max(prof);
plot(sm(ipk), Tb+sdes*pk, 'ro','MarkerFaceColor','r');
text(sm(ipk)+0.6, Tb+sdes*pk, sprintf(' peak %.1f \\circC at %.2f mm', ...
     Tb+sdes*pk, sm(ipk)));
plot(0, Tb, 'bs','MarkerFaceColor','b');
text(0.4, Tb+3, ' surface pinned at 37 \circC by the blood');
grid on; xlim([0 8]); xlabel('depth below the surface  (mm)'); ylabel('T (\circC)');
title({'The cooled surface pushes the peak inward', ...
       'the catheter lies here and is cooled too, so it adds no hot spot'});
exportgraphics(f, fullfile(fd,'fig5_depth_profiles.png'),'Resolution',160); close(f);

% ---- Fig 6: thermal dose ------------------------------------------------
texp = 60;
kk = find(t <= texp);
Tser = Tb + sdes*Th_y0(:,:,kk);                 % [nx x nz x nk]
D = cem43(t(kk), Tser);                         % [nx x nz]
f = figure('Visible','off','Position',[100 100 1100 430]);
subplot(1,2,1);
imagesc(xm, sm, log10(max(D.',1e-3))); set(gca,'YDir','reverse');
clim([-2 6]); colormap(parula); cb = colorbar;
cb.Label.String = 'log_{10} CEM43  (equivalent minutes at 43 \circC)';
hold on;
[c1,h1] = contour(xm, sm, D.', [60 240], 'w-','LineWidth',2);
clabel(c1,h1,'Color','w','FontSize',9);
plot([0 Lant],[0 0],'r-','LineWidth',4);
xlabel('x (mm)'); ylabel('depth (mm)');
title(sprintf('Thermal dose after %g s at %.2f W', texp, Pdes));

subplot(1,2,2);
Pset = logspace(log10(0.2*Pdes), log10(3*Pdes), 22);
dep = zeros(size(Pset)); wid = zeros(size(Pset));
for p = 1:numel(Pset)
    sp = Pset(p)*eta/Pref;
    Dp = cem43(t(kk), Tb + sp*Th_y0(:,:,kk));
    m  = Dp >= 240;
    if any(m(:))
        [~,jj] = find(m);  dep(p) = max(sm(jj));
        [ii,~] = find(m);  wid(p) = max(xm(ii)) - min(xm(ii));
    end
end
plot(Pset, dep, 'LineWidth',2); hold on; plot(Pset, wid, 'LineWidth',2);
grid on; xlabel('generator forward power  (W)'); ylabel('mm');
legend({'lesion depth (CEM43 = 240 min)','lesion width'},'Location','northwest');
title(sprintf('Ablation-equivalent zone after %g s', texp));
exportgraphics(f, fullfile(fd,'fig6_thermal_dose.png'),'Resolution',160); close(f);

% ---- Fig 7: where the power goes ----------------------------------------
f = figure('Visible','off','Position',[100 100 820 460]);
yyaxis left;
semilogx(t(2:end), 100*R.P_loss(2:end)/R.P_in, 'LineWidth',2); hold on;
ylabel('fraction of the deposited power leaving through the boundaries  (%)');
ylim([0 105]); yline(100,'k:');
yyaxis right;
semilogx(t(2:end), R.U(2:end)/R.P_in, 'LineWidth',2);
ylabel('stored energy / deposited power  (s)');
grid on; xlabel('time  (s)');
title(sprintf(['Energy audit at the reference %.0f W delivered: ' ...
    '%.3f W deposited in the block'], Pref, R.P_in));
exportgraphics(f, fullfile(fd,'fig7_energy.png'),'Resolution',160); close(f);

% ---- Fig 8: pulsed vs CW at the same average power ----------------------
f = figure('Visible','off','Position',[100 100 900 470]);
tq = linspace(0, 60, 3000);
kern = R.Th_max;
Pavg = Pdes;
cases = {@(tt) Pavg*ones(size(tt)),  sprintf('CW at %.2f W', Pavg); ...
         @(tt) 2*Pavg*(mod(tt,2) < 1), sprintf('50 %% duty, 0.5 Hz, %.2f W peak', 2*Pavg); ...
         @(tt) 4*Pavg*(mod(tt,8) < 2), sprintf('25 %% duty, 0.125 Hz, %.2f W peak', 4*Pavg)};
co = lines(3);
for p = 1:3
    th = duhamel_fig(kern, t, cases{p,1}, tq, eta, Pref);
    plot(tq, Tb + th, 'LineWidth',1.6, 'Color', co(p,:)); hold on;
end
yline(100,'k--','boiling'); yline(43,'k:','43 \circC');
grid on; xlabel('time  (s)'); ylabel('peak temperature  (\circC)');
legend(cases(:,2),'Location','southeast','Interpreter','none');
title({sprintf('Pulsed drive at the same average power, %.2f W', Pavg), ...
       'duty cycling raises the peak, it does not lower it'});
exportgraphics(f, fullfile(fd,'fig8_pulsed.png'),'Resolution',160); close(f);

fprintf('wrote 8 figures to Heatmodel/figures/\n');
end

% =========================================================================
function th = duhamel_fig(kern, tk, Pfun, tq, eta, Pref)
tmax = max(tq);  nu = 6000;  dtu = tmax/nu;  tu = (0:nu)*dtu;
Pm = 0.5*(Pfun(tu(1:end-1)) + Pfun(tu(2:end)))*eta/Pref;
g  = interp1(tk, kern, tu, 'linear');  g(1) = 0;
y  = conv(Pm, diff(g));
y  = [0 y(1:nu)];
th = interp1(tu, y, tq, 'linear');
end

