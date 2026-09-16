function [best, tab] = gel_heating_optimise(R, varargin)
%GEL_HEATING_OPTIMISE  Balance water cooling against microwave heating so
%                      that the surface layer stays below T_protect while a
%                      lesion above T_lesion forms beneath it -- and make the
%                      protected layer as thick as possible.
%
%   best = GEL_HEATING_OPTIMISE(R)          R from gel_heating_solve
%   GEL_HEATING_OPTIMISE(R)                 no output -> draws ONE figure
%   [best, tab] = ...                       tab: every cooling setting tried
%
%   Options (name/value):
%     'T_w'       water temperatures tried        default 5:2.5:R.T_bath0  [degC]
%     't_pre'     pre-cooling times tried         default [0 60 120 300]  [s]
%     't_exp'     exposure times tried            default [30 60 120 180 240 300 360]
%     'T_protect' surface layer must stay below   default 45  [degC]
%     'T_lesion'  lesion must reach               default 50  [degC]
%     'T_max'     nothing may exceed              default 95  [degC]
%     'P_max'     generator ceiling               default 150 [W]
%     'dt_check'  the surface constraint is checked every dt_check seconds
%                 of the exposure, not just at the end   default 10 [s]
%
%   THE ALGORITHM
%   The protected depth D(P) shrinks as the power rises and the lesion grows,
%   so for a given cooling setting (T_w, t_pre, t_exp) the best power is the
%   SMALLEST one that still makes a lesion.  That is found by bisection:
%     feasible(P) :  D_min = min over the exposure of the protected depth
%                    (the layer that stayed <= T_protect the whole time),
%                    lesion = max T at the end of exposure at depth >= D_min,
%                    feasible if lesion >= T_lesion and peak <= T_max.
%   Then the cooling setting with the largest D_min wins (ties: lower power,
%   shorter exposure).  Everything is instant because the fields are linear
%   combinations of the stored unit responses (gel_heating_solve).

o = struct('T_w', 5:2.5:R.T_bath0, 't_pre', [0 60 120 300], ...
           't_exp', [30 60 120 180 240 300 360], ...
           'T_protect', 45, 'T_lesion', 50, 'T_max', 95, 'P_max', 150, ...
           'dt_check', 10, 'verbose', nargout == 0);
for q1 = 1:2:numel(varargin), o.(varargin{q1}) = varargin{q1+1}; end

rows = {};
for T_w = o.T_w
  for t_pre = o.t_pre
    for t_exp = o.t_exp
        if t_pre + t_exp > R.t(end), continue; end
        tc = unique([0:o.dt_check:t_exp, t_exp]);
        % bisection on P for the smallest feasible power
        [okHi, mHi] = feasible(o.P_max);
        if ~okHi, continue; end
        lo = 0;  hi = o.P_max;  m = mHi;
        for it = 1:14
            mid = 0.5*(lo + hi);
            [ok, mm] = feasible(mid);
            if ok, hi = mid;  m = mm;  else, lo = mid; end
        end
        rows(end+1, :) = {T_w, t_pre, t_exp, hi, m.D_min, m.T_les, m.T_peak, m.T_cath}; %#ok<AGROW>
    end
  end
end
tab = cell2table(rows, 'VariableNames', ...
      {'T_w','t_pre','t_exp','P_set','D_mm','T_lesion','T_peak','T_cath'});
tab.D_mm = tab.D_mm*1e3;
if isempty(tab)
    error('gel_heating_optimise:infeasible', ...
          'no setting makes a %g degC lesion below a %g degC layer', o.T_lesion, o.T_protect);
end
tab = sortrows(tab, {'D_mm','P_set','t_exp'}, {'descend','ascend','ascend'});
best = table2struct(tab(1, :));
best.P_del = R.eta*best.P_set;
[best.T, best.M] = gel_heating_eval(R, best.P_set, best.T_w, best.t_pre, best.t_exp, ...
                                    o.T_protect, o.T_lesion);
best.opts = o;

if o.verbose
    fprintf(['optimum: water %.1f degC, pre-cool %.0f s, P_set %.1f W (%.1f W delivered), ' ...
             'exposure %.0f s\n  protected layer %.2f mm (<= %g degC throughout), lesion %.1f degC ' ...
             'below it, peak %.1f degC, catheter %.1f degC\n'], best.T_w, best.t_pre, best.P_set, ...
             best.P_del, best.t_exp, best.D_mm, o.T_protect, best.T_lesion, best.T_peak, best.T_cath);
end
if nargout == 0, draw(R, best, tab, o); end

    % ------------------------------------------------------------------
    function [ok, m] = feasible(P)
        D_min = Inf;  T_pk = -Inf;
        for k = 1:numel(tc)
            [~, M] = gel_heating_eval(R, P, T_w, t_pre, tc(k), o.T_protect, o.T_lesion);
            D_min = min(D_min, M.D);  T_pk = max(T_pk, M.T_peak);
        end
        [T, M] = gel_heating_eval(R, P, T_w, t_pre, t_exp, o.T_protect, o.T_lesion);
        depth = -R.z(:).';
        if isfinite(D_min), T_les = max(T(:, depth >= D_min - 1e-12), [], 'all');
        else,               T_les = NaN;   end
        ok = isfinite(D_min) && T_les >= o.T_lesion && T_pk <= o.T_max;
        m  = struct('D_min', D_min, 'T_les', T_les, 'T_peak', T_pk, 'T_cath', M.T_cath);
    end
end

% =========================================================================
function draw(R, best, tab, o)
f = figure('Name','surface sparing','Color','w','Position',[40 60 1500 520]);
tl = tiledlayout(f, 1, 3, 'TileSpacing','compact', 'Padding','compact');
title(tl, sprintf(['Surface sparing: keep the top layer <= %g \\circC, make a >= %g \\circC ' ...
    'lesion below it, maximise the layer thickness'], o.T_protect, o.T_lesion), 'FontWeight','bold');
x = R.x*1e3;  z = R.z*1e3;  T = best.T.';  M = best.M;

ax = nexttile(tl);  hold(ax,'on');
imagesc(ax, x, z, T);  set(ax,'YDir','normal');  axis(ax,'image');
colormap(ax, turbo);  clim(ax, [min(T(:)) max(T(:))]);
contour(ax, x, z, T, [o.T_protect o.T_protect], 'w--', 'LineWidth', 1.2);
contour(ax, x, z, T, [o.T_lesion o.T_lesion], 'w-', 'LineWidth', 1.5);
yline(ax, -M.D*1e3, 'c-', sprintf('protected %.1f mm', M.D*1e3), 'LineWidth', 1.2, ...
      'LabelHorizontalAlignment','right');
plot(ax, [0 R.L_el]*1e3, [1 1]*R.h_ax*1e3, '-', 'Color', [.3 .3 .3], 'LineWidth', 3);
ylim(ax, [-30 2]);  xlabel(ax,'x (mm)');  ylabel(ax,'z (mm)');
cb = colorbar(ax);  cb.Label.String = 'T (\circC)';
title(ax, sprintf(['water %.1f \\circC, pre-cool %.0f s, P_{set} = %.0f W (%.1f W), %.0f s\n' ...
    'dashed %g \\circC, solid %g \\circC'], best.T_w, best.t_pre, best.P_set, best.P_del, ...
    best.t_exp, o.T_protect, o.T_lesion));

ax = nexttile(tl);  hold(ax,'on');
plot(ax, best.T(M.ix_peak, :), z, 'r-', 'LineWidth', 1.8);
xline(ax, o.T_protect, 'k--', sprintf('%g', o.T_protect));
xline(ax, o.T_lesion,  'k-',  sprintf('%g', o.T_lesion));
yline(ax, -M.D*1e3, 'c-', 'protected layer', 'LineWidth', 1.2);
xlabel(ax,'T (\circC)');  ylabel(ax,'z (mm)');  grid(ax,'on');  ylim(ax, [-30 0]);
title(ax, sprintf('depth profile through the hottest column, x = %.1f mm', x(M.ix_peak)));

ax = nexttile(tl);  hold(ax,'on');
tb = tab(tab.t_pre == best.t_pre, :);
for te = unique(tb.t_exp).'
    s = tb(tb.t_exp == te, :);  s = sortrows(s, 'T_w');
    plot(ax, s.T_w, s.D_mm, '-o', 'MarkerSize', 4, 'DisplayName', sprintf('%.0f s exposure', te));
end
plot(ax, best.T_w, best.D_mm, 'kp', 'MarkerSize', 14, 'MarkerFaceColor', 'y', 'DisplayName', 'optimum');
xlabel(ax,'water temperature (\circC)');  ylabel(ax,'protected depth (mm)');  grid(ax,'on');
legend(ax, 'Location','best');
title(ax, sprintf('protected depth at the least power that still makes a lesion (pre-cool %.0f s)', best.t_pre));
exportgraphics(f, fullfile(fileparts(mfilename('fullpath')), 'gel_sparing.png'), 'Resolution', 150);
end
