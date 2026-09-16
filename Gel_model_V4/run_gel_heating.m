%RUN_GEL_HEATING  Microwave heating of the gel for the photographed setup
%                 (model_defaults): ONE figure, gel_heating.png.
%
%   Chain:  E_field_2d (E_rms on the plane containing the antenna, medium
%           from the insulation radius outward)
%           -> q = sig_eff E_rms^2  (axisymmetric about the wire)
%           -> 3-D transient conduction, Robin faces (or pinned top), from the
%              gel in equilibrium with the bath at T_bath0
%           -> T = T_ss + (P_set/P_ref) Theta(t)  (linear, one solve per
%              set of constants: gel_heating_solve, cached in solutions/)
%           -> catheter: lumped energy balance, cylinder in cross-flow.
%
%   Panels: the antenna-plane map at the end of the photographed run (50 W,
%   332 s), the peak temperature vs time for three powers, and the peak and
%   catheter temperature vs generator setting.  Constants from
%   model_defaults.m; change them there or in the GUI's Constants tab.
%
%   Run from this folder:  >> run_gel_heating

P_list = [25 50 100];           % [W]  generator settings shown
t_end  = 332;                   % [s]  the photographed series (IMG_7907 -> 7918)

R = gel_heating_solve();        % model_defaults, cached
Pc = R.P;  D = R.D;  c = R.cath;
T_w = R.T_bath0;                % water as the gel started: everything at 25 degC
x = R.x*1e3;  z = R.z*1e3;
t_hist = R.t(R.t <= t_end);
Tpk = zeros(numel(P_list), numel(t_hist));  Tc = zeros(size(P_list));
for i = 1:numel(P_list)
    for j = 1:numel(t_hist)
        [~, M] = gel_heating_eval(R, P_list(i), T_w, 0, t_hist(j));
        Tpk(i, j) = M.T_peak;
    end
    Tc(i) = M.T_cath;
end
[T50, M50] = gel_heating_eval(R, 50, T_w, 0, t_end);
i92 = find(t_hist >= 92, 1);

fprintf('antenna: %g mm element, axis %.2f mm above the gel, k_e = %.1f%+.1fj 1/m (%s), braid arm %g mm\n', ...
        Pc.L*1e3, Pc.h_axis*1e3, real(D.ke), imag(D.ke), D.ke_source, Pc.L_sheath*1e3);
fprintf('E-field: %.1f W delivered -> gel block %.2f W (%.0f %%), water %.2f W, film round catheter %.3f W\n', ...
        R.eta*R.P_ref, R.P_gel, 100*R.P_gel/(R.eta*R.P_ref), R.P_water, R.P_film);
fprintf('boundaries: h_top %s (%s), h_wall %.1f W/m^2K; gel starts at %g degC, room %g degC\n', ...
        num2str(D.h_top, '%.0f'), D.h_top_source, D.h_wall, R.T_bath0, R.T_room);
fprintf('grid %dx%dx%d, dh %g mm, dt %.4f s, solved in %.0f s\n', D.nx, D.ny, D.nz, Pc.dh*1e3, D.dt, R.t_solve);
fprintf('%6s %9s %11s %11s %11s\n', 'P_set', 'P_del', 'T_peak@92s', 'T_peak@end', 'T_cath');
for i = 1:numel(P_list)
    fprintf('%6.0f %9.1f %11.1f %11.1f %11.2f\n', P_list(i), R.eta*P_list(i), Tpk(i, i92), Tpk(i, end), Tc(i));
end
fprintf('peak (50 W, %g s) at x = %.1f mm, z = %.1f mm;  photographs: x ~ 3-3.6 mm, z ~ -3.3..-3.6 mm\n', ...
        t_end, x(M50.ix_peak), z(M50.iz_peak));

% ---------------------------------------------------------------- FIGURE --
f = figure('Name', 'gel heating V4', 'Color', 'w', 'Position', [40 60 1560 540]);
tl = tiledlayout(f, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf(['Microwave heating of the gel -- %g mm element, axis %.1f mm above the gel, ' ...
    'k_e = %.0f%+.0fj m^{-1}; gel and %g \\circC water at %g L/min; tank walls to %g \\circC room'], ...
    Pc.L*1e3, Pc.h_axis*1e3, real(D.ke), imag(D.ke), T_w, Pc.Q_circ*60000, R.T_room), 'FontWeight', 'bold');

ax = nexttile(tl);  ax.Toolbar.Visible = 'off';
Tpl = T50.';
imagesc(ax, x, z, Tpl);  set(ax, 'YDir', 'normal');  axis(ax, 'image');  hold(ax, 'on');
colormap(ax, hot);  clim(ax, [T_w max(Tpl(:))]);
lv = [30 35 40 45 50 60 70 80];  lv = lv(lv > T_w & lv < max(Tpl(:)));
if ~isempty(lv)
    [C, hc] = contour(ax, x, z, Tpl, lv, 'c-', 'LineWidth', 0.7);
    clabel(C, hc, 'FontSize', 8, 'Color', 'c');
end
plot(ax, [0 R.L_el]*1e3, [1 1]*R.h_ax*1e3, '-', 'Color', [.3 .3 .3], 'LineWidth', 3);
plot(ax, x(M50.ix_peak), z(M50.iz_peak), 'wo', 'MarkerSize', 8, 'LineWidth', 1.5);
ylim(ax, [z(1) 2]);  xlabel(ax, 'x (mm)');  ylabel(ax, 'z (mm)');
cb = colorbar(ax);  cb.Label.String = 'T (\circC)';
title(ax, sprintf('T on the antenna plane at t = %.0f s, P_{set} = 50 W (%g W delivered)', t_end, R.eta*50));

ax = nexttile(tl);  ax.Toolbar.Visible = 'off';  hold(ax, 'on');
for i = 1:numel(P_list)
    plot(ax, t_hist, Tpk(i, :), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('P_{set} = %g W  (%g W)', P_list(i), R.eta*P_list(i)));
end
yline(ax, 100, 'k--', 'boiling', 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
yline(ax, T_w, 'k:', sprintf('%g \\circC water', T_w), 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
xlabel(ax, 'time (s)');  ylabel(ax, 'peak gel temperature (\circC)');  grid(ax, 'on');
legend(ax, 'Location', 'southeast');  title(ax, 'Peak gel temperature vs time');

ax = nexttile(tl);  ax.Toolbar.Visible = 'off';  hold(ax, 'on');
Pp = 0:5:150;  Tp_end = zeros(size(Pp));  Tp_92 = Tp_end;  Tcp = Tp_end;
for i = 1:numel(Pp)
    [~, M] = gel_heating_eval(R, Pp(i), T_w, 0, t_end);   Tp_end(i) = M.T_peak;  Tcp(i) = M.T_cath;
    [~, M] = gel_heating_eval(R, Pp(i), T_w, 0, 92);      Tp_92(i)  = M.T_peak;
end
plot(ax, Pp, Tp_end, 'r-', 'LineWidth', 1.8, 'DisplayName', sprintf('gel peak, t = %.0f s', t_end));
plot(ax, Pp, Tp_92, 'r--', 'LineWidth', 1.2, 'DisplayName', 'gel peak, t = 92 s');
plot(ax, Pp, Tcp, 'b-', 'LineWidth', 1.8, 'DisplayName', 'catheter (film + gap conduction)');
yline(ax, 100, 'k--', 'HandleVisibility', 'off');
xlabel(ax, sprintf('generator setting P_{set} (W)   [delivered = %g %%]', 100*R.eta));  ylabel(ax, 'temperature (\circC)');
grid(ax, 'on');  legend(ax, 'Location', 'northwest');
title(ax, 'Temperature vs supplied power');

exportgraphics(f, fullfile(fileparts(mfilename('fullpath')), 'gel_heating.png'), 'Resolution', 150);
