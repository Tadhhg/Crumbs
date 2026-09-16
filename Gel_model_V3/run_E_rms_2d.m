%RUN_E_RMS_2D  ONE figure: the time-RMS electric field E_rms(x,z) on the gel
%              plane -- the quantity that sets the heating, q = sig_eff*E_rms^2.
%
%   Uses E_field_2d.m.  Run from this folder:   >> run_E_rms_2d
%   Same frame as run_E_field_2d: feed at x = 0, tip at x = 27 mm, gel from
%   x = -1.5 to 28.5 mm, z = 0 is the gel surface.

P_set = 100;                                % [W] generator setting -> 20 W delivered
dx    = 0.1e-3;                             % [m] 0.1 mm -> 301 x 301 map

% ------------------------------------------------------------------ MAP ---
x = (-1.5e-3 : dx : 28.5e-3);   z = (-30e-3 : dx : 0);
[X, Z] = meshgrid(x, z);
[~, S] = E_field_2d(X, Z, 0, P_set);
Erms = S.E_rms;                             % [V/m]

% ---------------------------------------------------------------- FIGURE --
f = figure('Name','E_rms','Color','w','Position',[80 80 900 820]);
ax = axes(f);  ax.Toolbar.Visible = 'off';
imagesc(ax, x*1e3, z*1e3, log10(Erms));
set(ax, 'YDir','normal');  axis(ax, 'image');  hold(ax, 'on');
colormap(ax, turbo);  clim(ax, log10([30 10000]));
levels = [50 70 100 150 200 300 500 700 1000 1500 2000 3000 5000];
[C, hc] = contour(ax, x*1e3, z*1e3, Erms, levels, 'k-', 'LineWidth', 0.6);
clabel(C, hc, 'FontSize', 8, 'Color', 'k', 'LabelSpacing', 400);
plot(ax, [0 27], [0.5 0.5], 'k-', 'LineWidth', 4);                        % element
plot(ax, [x(1) 0]*1e3, [0.5 0.5], '-', 'Color', [.5 .5 .5], 'LineWidth', 4); % braid
ylim(ax, [-30 2]);  xlabel(ax, 'x (mm)');  ylabel(ax, 'z (mm)');
cb = colorbar(ax, 'Ticks', log10([30 50 100 200 300 500 1000 2000 3000 5000 10000]), ...
              'TickLabels', {'30','50','100','200','300','500','1k','2k','3k','5k','10k'});
cb.Label.String = 'E_{rms}  (V/m)';
title(ax, {sprintf(['Time-RMS E-field in the gel plane, P_{set} = %g W (%g W delivered), ' ...
                    '2.45 GHz'], P_set, S.P_del), ...
           sprintf(['heating rate q = \\sigma_{eff} E_{rms}^2 with \\sigma_{eff} = %.2f S/m;  ' ...
                    'peak %.0f V/m at the surface under the feed'], S.sig_eff, max(Erms(:)))});

exportgraphics(f, 'E_rms.png', 'Resolution', 300);
