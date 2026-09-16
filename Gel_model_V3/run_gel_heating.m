%RUN_GEL_HEATING  Microwave heating of the gel with convective boundaries,
%                 and the convective cooling of the catheter.  ONE figure.
%
%   Chain:  E_field_2d (E_rms on the plane containing the antenna)
%           -> q = sig_eff E_rms^2  (axisymmetric about the wire, so the plane
%              gives the whole 3-D source)
%           -> 3-D transient conduction in the 30 x 30 x 30 mm gel, starting
%              from the STEADY STATE of gel_heat_2d (37 degC water above,
%              tank walls to 25 degC room), same convective boundaries
%           -> T(t) = T_steady + Theta(t) * P_del/P_ref  (the heat equation is
%              linear and q is proportional to delivered power, so one solve
%              covers every power)
%           -> catheter: lumped energy balance, cylinder in cross-flow.
%   The solve itself lives in gel_heating_solve.m (run once, cached in
%   gel_heating_solution.mat); gel_heating_gui.m animates it and
%   gel_heating_optimise.m balances cooling against heating.
%
%   WHY 3-D AND NOT THE 2-D PLANE.  The source is concentrated within a few mm
%   of the wire in y as well as in z (near field ~ 1/rho).  A 2-D solve in the
%   y = 0 plane treats the source as infinitely long in y and ignores the
%   sideways spreading, which over-predicts the rise substantially
%   (../Heatmodel/Heat_model_writeup.md s7).  The figure still shows the y = 0
%   plane.
%
%   CATHETER (no CFD).  The antenna's own losses are negligible (ohmic
%   ~0.01 W), so the catheter is heated by (i) the microwave power the field
%   model deposits in the thin water film clinging to it, P_film, and (ii)
%   conduction from the hot gel surface across the water gap, and cooled by
%   the 5 L/min cross-flow, h_c from Churchill-Bernstein (Incropera Ch. 7, as
%   in ../Heatmodel/matlab/convection_h.m):
%       h_c A_c (T_c - T_w) = P_film + G_gap (T_s - T_c),   G_gap = k_w D L/g
%   Everything on the right is proportional to P, so T_c - T_w is too.
%
%   Run from this folder:  >> run_gel_heating

P_list = [50 100 150];          % [W]  generator settings shown (10, 20, 30 W delivered)
t_end  = 360;                   % [s]  exposure (the photographed run was 6 min)
T_w    = 37;                    % [degC] water, as on the bench

R = gel_heating_solve();
x = R.x*1e3;  z = R.z*1e3;  c = R.cath;
t_hist = R.t(R.t <= t_end);
Tpk = zeros(numel(P_list), numel(t_hist));  Tc = zeros(size(P_list));
for i = 1:numel(P_list)
    for j = 1:numel(t_hist)
        [~, M] = gel_heating_eval(R, P_list(i), T_w, 0, t_hist(j));
        Tpk(i, j) = M.T_peak;
    end
    Tc(i) = M.T_cath;
end
[T100, M100] = gel_heating_eval(R, 100, T_w, 0, t_end);
i60 = find(t_hist >= 60, 1);

fprintf('E-field: %.1f W delivered -> gel block %.1f W (%.0f %%), water %.1f W, film round catheter %.2f W\n', ...
        R.eta*R.P_ref, R.P_gel, 100*R.P_gel/(R.eta*R.P_ref), R.P_water, R.P_film);
fprintf('steady state (no power): gel %.1f..%.1f degC\n', R.T_ss_range);
fprintf('catheter: Re_D %.0f, Nu %.1f, h_c %.0f W/m^2K, film %.2f mm, G_gap %.3f W/K, h_c A_c %.3f W/K\n', ...
        c.Re, c.Nu, c.h_c, c.d_t*1e3, c.G_gap, c.h_c*c.A_c);
fprintf('%6s %9s %11s %11s %11s\n', 'P_set', 'P_del', 'T_peak@60s', 'T_peak@end', 'T_cath');
for i = 1:numel(P_list)
    fprintf('%6.0f %9.1f %11.1f %11.1f %11.2f\n', P_list(i), R.eta*P_list(i), Tpk(i, i60), Tpk(i, end), Tc(i));
end
fprintf('peak (100 W) at x = %.1f mm, z = %.1f mm\n', x(M100.ix_peak), z(M100.iz_peak));

% ---------------------------------------------------------------- FIGURE --
f = figure('Name','gel heating','Color','w','Position',[40 60 1560 540]);
tl = tiledlayout(f, 1, 3, 'TileSpacing','compact', 'Padding','compact');
title(tl, sprintf(['Microwave heating of the gel (27 mm monopole, 2.45 GHz, catheter %.1f mm ' ...
    'above the gel) -- %g \\circC water at 5 L/min above, tank walls to 25 \\circC room'], ...
    R.gap*1e3, T_w), 'FontWeight','bold');

ax = nexttile(tl);  ax.Toolbar.Visible = 'off';
Tpl = T100.';
imagesc(ax, x, z, Tpl);  set(ax,'YDir','normal');  axis(ax,'image');  hold(ax,'on');
colormap(ax, hot);  clim(ax, [T_w max(Tpl(:))]);
[C, hc] = contour(ax, x, z, Tpl, [40 45 50 60 70 80 90 100], 'c-', 'LineWidth', 0.7);
clabel(C, hc, 'FontSize', 8, 'Color', 'c');
plot(ax, [0 R.L_el]*1e3, [1 1]*R.h_ax*1e3, '-', 'Color', [.3 .3 .3], 'LineWidth', 3);
ylim(ax, [-30 2]);  xlabel(ax,'x (mm)');  ylabel(ax,'z (mm)');
cb = colorbar(ax);  cb.Label.String = 'T (\circC)';
title(ax, sprintf('T on the antenna plane at t = %.0f s, P_{set} = 100 W (%g W delivered)', ...
                  t_end, R.eta*100));

ax = nexttile(tl);  ax.Toolbar.Visible = 'off';  hold(ax,'on');
for i = 1:numel(P_list)
    plot(ax, t_hist, Tpk(i, :), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('P_{set} = %g W  (%g W)', P_list(i), R.eta*P_list(i)));
end
yline(ax, 100, 'k--', 'boiling', 'LabelHorizontalAlignment','left', 'HandleVisibility','off');
yline(ax, T_w, 'k:', sprintf('%g \\circC water', T_w), 'LabelHorizontalAlignment','left', 'HandleVisibility','off');
xlabel(ax,'time (s)');  ylabel(ax,'peak gel temperature (\circC)');  grid(ax,'on');
legend(ax, 'Location','southeast');  title(ax, 'Peak gel temperature vs time');

ax = nexttile(tl);  ax.Toolbar.Visible = 'off';  hold(ax,'on');
Pp = 0:5:150;  Tp_end = zeros(size(Pp));  Tp_60 = Tp_end;  Tcp = Tp_end;
for i = 1:numel(Pp)
    [~, M] = gel_heating_eval(R, Pp(i), T_w, 0, t_end);   Tp_end(i) = M.T_peak;  Tcp(i) = M.T_cath;
    [~, M] = gel_heating_eval(R, Pp(i), T_w, 0, 60);      Tp_60(i)  = M.T_peak;
end
plot(ax, Pp, Tp_end, 'r-', 'LineWidth', 1.8, 'DisplayName', sprintf('gel peak, t = %.0f s', t_end));
plot(ax, Pp, Tp_60, 'r--', 'LineWidth', 1.2, 'DisplayName', 'gel peak, t = 60 s');
plot(ax, Pp, Tcp, 'b-', 'LineWidth', 1.8, 'DisplayName', 'catheter (film + gap conduction)');
plot(ax, Pp, T_w + Pp/R.P_ref*R.P_water/(c.h_c*c.A_c), 'b:', 'LineWidth', 1.2, ...
     'DisplayName', 'catheter, bound: all water-side power');
yline(ax, 100, 'k--', 'HandleVisibility','off');
xlabel(ax,'generator setting P_{set} (W)   [delivered = 20 %]');  ylabel(ax,'temperature (\circC)');
grid(ax,'on');  legend(ax, 'Location','northwest');
title(ax, 'Temperature vs supplied power');

exportgraphics(f, 'gel_heating.png', 'Resolution', 150);
