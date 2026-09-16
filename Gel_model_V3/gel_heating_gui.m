function fig = gel_heating_gui(R)
%GEL_HEATING_GUI  Interactive animation of the gel heating.
%
%   >> gel_heating_gui
%
%   Sliders:  generator power, water temperature, pre-cooling time, exposure
%   time, and a time scrub.  Play animates the exposure.  "Optimise" runs
%   gel_heating_optimise (least power that still makes a >= 50 degC lesion
%   below a layer kept <= 45 degC, deepest such layer wins) and sets the
%   sliders to the answer.
%
%   Everything is instant because the fields are linear combinations of the
%   unit responses stored by gel_heating_solve (first call solves, ~2 min).

if nargin < 1, R = gel_heating_solve(); end
T_protect = 45;  T_lesion = 50;
x = R.x*1e3;  z = R.z*1e3;  depth = -R.z;

% ------------------------------------------------------------- STATE ----
S = struct('P', 100, 'T_w', 37, 't_pre', 0, 't_exp', 360, 't', 0, 'playing', false);
t_max_total = R.t(end);

% ----------------------------------------------------------- LAYOUT -----
fig = uifigure('Name', 'Gel heating', 'Position', [40 40 1440 780], 'Color', 'w');
gl  = uigridlayout(fig, [1 2], 'ColumnWidth', {'1.55x', '1x'});
axMap = uiaxes(gl);  axMap.Layout.Column = 1;
right = uigridlayout(gl, [3 1], 'RowHeight', {'1x', '1x', 390}, 'Padding', [0 0 0 0]);
right.Layout.Column = 2;
axProf = uiaxes(right);  axHist = uiaxes(right);
ctrl = uigridlayout(right, [7 3], 'ColumnWidth', {150, '1x', 90}, ...
                    'RowHeight', {36 36 36 36 36 40 '1x'}, 'RowSpacing', 4);

sl = struct();
sl.P     = addslider(1, 'Generator power P_{set} (W)', [0 150], S.P);
sl.T_w   = addslider(2, 'Water temperature (\circC)', [2 37], S.T_w);
sl.t_pre = addslider(3, 'Pre-cool before power (s)', [0 300], S.t_pre);
sl.t_exp = addslider(4, 'Exposure time (s)', [10 420], S.t_exp);
sl.t     = addslider(5, 'Time since power on (s)', [0 S.t_exp], 0);
sl.P.MajorTicks = 0:25:150;         sl.T_w.MajorTicks = [2 5:5:35 37];
sl.t_pre.MajorTicks = 0:60:300;     sl.t_pre.MinorTicks = [];
sl.t_exp.MajorTicks = 60:60:720;    sl.t_exp.MinorTicks = [];
setTimeTicks();

btnRow = uigridlayout(ctrl, [1 4], 'Padding', [0 0 0 0]);
btnRow.Layout.Row = 6;  btnRow.Layout.Column = [1 3];
bPlay = uibutton(btnRow, 'Text', 'Play', 'ButtonPushedFcn', @(~,~) togglePlay());
uibutton(btnRow, 'Text', 'Reset', 'ButtonPushedFcn', @(~,~) setTime(0));
uibutton(btnRow, 'Text', 'Optimise', 'ButtonPushedFcn', @(~,~) optimise());
uibutton(btnRow, 'Text', 'Baseline', 'ButtonPushedFcn', @(~,~) baseline());
readout = uilabel(ctrl, 'Text', '', 'FontName', 'Consolas', 'FontSize', 11, ...
                  'VerticalAlignment', 'top', 'WordWrap', 'on');
readout.Layout.Row = 7;  readout.Layout.Column = [1 3];

% ------------------------------------------------------------ MAP -------
hold(axMap, 'on');
im = imagesc(axMap, x, z, zeros(numel(z), numel(x)));
set(axMap, 'YDir', 'normal');  axis(axMap, 'image');  ylim(axMap, [-30 2]);
colormap(axMap, turbo);  clim(axMap, [20 90]);
cb = colorbar(axMap);  cb.Label.String = 'T (\circC)';
[~, c45] = contour(axMap, x, z, zeros(numel(z), numel(x)), [T_protect T_protect], 'w--', 'LineWidth', 1.2);
[~, c50] = contour(axMap, x, z, zeros(numel(z), numel(x)), [T_lesion  T_lesion],  'w-',  'LineWidth', 1.6);
hD   = plot(axMap, x([1 end]), [0 0], 'c-', 'LineWidth', 1.2);
hPk  = plot(axMap, 0, 0, 'wo', 'MarkerSize', 8, 'LineWidth', 1.5);
plot(axMap, [0 R.L_el]*1e3, [1 1]*R.h_ax*1e3, '-', 'Color', [.3 .3 .3], 'LineWidth', 4);
xlabel(axMap, 'x (mm)');  ylabel(axMap, 'z (mm)');
ttl = title(axMap, '');
warning('off', 'MATLAB:contour:ConstantData');

% -------------------------------------------------------- PROFILE -------
hold(axProf, 'on');
hProf = plot(axProf, zeros(size(z)), z, 'r-', 'LineWidth', 1.8);
xline(axProf, T_protect, 'k--', sprintf('%g', T_protect), 'LabelOrientation', 'horizontal');
xline(axProf, T_lesion,  'k-',  sprintf('%g', T_lesion),  'LabelOrientation', 'horizontal');
hProfD = yline(axProf, 0, 'c-', 'protected layer', 'LineWidth', 1.2);
xlim(axProf, [0 100]);  ylim(axProf, [-30 0]);  grid(axProf, 'on');
xlabel(axProf, 'T (\circC)');  ylabel(axProf, 'z (mm)');
ttlProf = title(axProf, '');

% -------------------------------------------------------- HISTORY -------
hold(axHist, 'on');
hTop = plot(axHist, 0, 0, 'b-', 'LineWidth', 1.6, 'DisplayName', 'max T in the top 3 mm');
hBel = plot(axHist, 0, 0, 'r-', 'LineWidth', 1.6, 'DisplayName', 'max T below 3 mm');
yline(axHist, T_protect, 'k--', 'HandleVisibility', 'off');
yline(axHist, T_lesion,  'k-',  'HandleVisibility', 'off');
hNow = xline(axHist, 0, 'k:', 'LineWidth', 1.2, 'HandleVisibility', 'off');
grid(axHist, 'on');  xlabel(axHist, 'time since power on (s)');  ylabel(axHist, 'T (\circC)');
legend(axHist, 'Location', 'southeast');
title(axHist, 'History for the current settings');

% ----------------------------------------------------------- TIMER ------
tmr = timer('ExecutionMode', 'fixedRate', 'Period', 0.05, 'BusyMode', 'drop', ...
            'TimerFcn', @(~,~) tick());
fig.CloseRequestFcn = @(~,~) closeApp();
fig.UserData = struct('getState', @getState, 'setScenario', @setScenario);

refreshHistory();
update();
if nargout == 0, clear fig; end

% =========================================================================
    function h = addslider(row, label, lim, val)
        lb = uilabel(ctrl, 'Text', label, 'Interpreter', 'tex');
        lb.Layout.Row = row;  lb.Layout.Column = 1;
        h = uislider(ctrl, 'Limits', lim, 'Value', val, ...
                     'ValueChangingFcn', @(s, e) sliderMoved(row, e.Value), ...
                     'ValueChangedFcn',  @(s, e) sliderMoved(row, e.Value));
        h.Layout.Row = row;  h.Layout.Column = 2;
        v = uilabel(ctrl, 'Text', '', 'HorizontalAlignment', 'right', 'Interpreter', 'tex');
        v.Layout.Row = row;  v.Layout.Column = 3;
        h.UserData = v;
    end

    function sliderMoved(row, val)
        switch row
            case 1, S.P = val;
            case 2, S.T_w = val;
            case 3
                S.t_pre = val;
                S.t_exp = min(S.t_exp, t_max_total - S.t_pre);
                sl.t_exp.Value = S.t_exp;  sl.t_exp.Limits = [10 t_max_total - S.t_pre];
                S.t = min(S.t, S.t_exp);  sl.t.Limits = [0 S.t_exp];  sl.t.Value = S.t;  setTimeTicks();
            case 4
                S.t_exp = val;
                S.t = min(S.t, S.t_exp);  sl.t.Limits = [0 S.t_exp];  sl.t.Value = S.t;  setTimeTicks();
            case 5
                S.t = val;
        end
        if row < 5, refreshHistory(); end
        update();
    end

    function setScenario(P, T_w, t_pre, t_exp, t)
        S.P = P;  S.T_w = T_w;  S.t_pre = t_pre;  S.t_exp = t_exp;  S.t = t;
        sl.P.Value = P;  sl.T_w.Value = T_w;  sl.t_pre.Value = t_pre;
        sl.t_exp.Limits = [10 t_max_total - t_pre];  sl.t_exp.Value = t_exp;
        sl.t.Limits = [0 t_exp];  sl.t.Value = t;  setTimeTicks();
        refreshHistory();  update();
    end

    function setTime(t)
        S.t = t;  sl.t.Value = t;  update();
    end

    function setTimeTicks()
        steps = [5 10 20 30 60 120];
        st = steps(find(S.t_exp./steps <= 12, 1));
        sl.t.MajorTicks = 0:st:S.t_exp;  sl.t.MinorTicks = [];
    end

    function baseline(), setScenario(100, 37, 0, 360, 0); end

    function s = getState(), s = S; end

    function togglePlay()
        if S.playing
            S.playing = false;  stop(tmr);  bPlay.Text = 'Play';
        else
            if S.t >= S.t_exp, S.t = 0; end
            S.playing = true;  bPlay.Text = 'Pause';  start(tmr);
        end
    end

    function tick()
        S.t = S.t + 3;                                   % 3 s of heating per tick
        if S.t >= S.t_exp
            S.t = S.t_exp;  S.playing = false;  stop(tmr);  bPlay.Text = 'Play';
        end
        sl.t.Value = S.t;
        update();
    end

    function optimise()
        dlg = uiprogressdlg(fig, 'Title', 'Optimising', 'Indeterminate', 'on', ...
                            'Message', 'least power that makes a 50 degC lesion below a 45 degC layer, for every cooling setting ...');
        try
            best = gel_heating_optimise(R, 'verbose', false);
            close(dlg);
            setScenario(best.P_set, best.T_w, best.t_pre, best.t_exp, best.t_exp);
        catch err
            close(dlg);
            uialert(fig, err.message, 'No feasible setting');
        end
    end

    function refreshHistory()
        tt = 0:2:S.t_exp;  top = zeros(size(tt));  bel = top;
        for k = 1:numel(tt)
            [~, M] = gel_heating_eval(R, S.P, S.T_w, S.t_pre, tt(k), T_protect, T_lesion);
            top(k) = M.T_top3;  bel(k) = M.T_below3;
        end
        set(hTop, 'XData', tt, 'YData', top);  set(hBel, 'XData', tt, 'YData', bel);
        % colour range: fixed during an animation, adapted to the settings
        Tend = gel_heating_eval(R, S.P, S.T_w, S.t_pre, S.t_exp);
        clim(axMap, [min(20, floor(S.T_w)) max(ceil(max(Tend(:))) + 2, 50)]);
        xlim(axHist, [0 max(S.t_exp, 10)]);
        ylim(axHist, [min([top bel 20]) max([top bel 60]) + 2]);
    end

    function update()
        [T, M] = gel_heating_eval(R, S.P, S.T_w, S.t_pre, S.t, T_protect, T_lesion);
        Tpl = T.';
        im.CData = Tpl;  c45.ZData = Tpl;  c50.ZData = Tpl;
        if isfinite(M.D), hD.YData = [-1 -1]*M.D*1e3;  hD.Visible = 'on';
        else,             hD.Visible = 'off';  end
        hPk.XData = x(M.ix_peak);  hPk.YData = z(M.iz_peak);
        ttl.String = sprintf(['t = %.0f s   P_{set} = %.0f W (%.1f W delivered),  water %.1f \\circC' ...
                              ', pre-cool %.0f s      dashed %g \\circC, solid %g \\circC'], ...
                              S.t, S.P, R.eta*S.P, S.T_w, S.t_pre, T_protect, T_lesion);
        hProf.XData = T(M.ix_peak, :);
        if isfinite(M.D), hProfD.Value = -M.D*1e3;  hProfD.Visible = 'on';
        else,             hProfD.Visible = 'off';  end
        xlim(axProf, [min(15, floor(min(T(:)))) max(60, ceil(M.T_peak) + 5)]);
        ttlProf.String = sprintf('depth profile at x = %.1f mm (hottest column)', x(M.ix_peak));
        hNow.Value = S.t;
        sl.P.UserData.Text     = sprintf('%.0f W', S.P);
        sl.T_w.UserData.Text   = sprintf('%.1f \\circC', S.T_w);
        sl.t_pre.UserData.Text = sprintf('%.0f s', S.t_pre);
        sl.t_exp.UserData.Text = sprintf('%.0f s', S.t_exp);
        sl.t.UserData.Text     = sprintf('%.0f s', S.t);
        if isfinite(M.D), Dtxt = sprintf('%.2f mm', M.D*1e3);  Ltxt = sprintf('%.1f', M.T_les);
        else,             Dtxt = 'whole gel';                  Ltxt = '--';  end
        readout.Text = sprintf([ ...
            'delivered power   %6.1f W  (%.0f %% of setting, %.0f %% of it into the gel)\n' ...
            'protected layer   %s  (T <= %g degC, min over x)\n' ...
            'max T below it    %s degC     %g degC first reached at %.2f mm\n' ...
            'peak              %6.1f degC at x = %.1f mm, z = %.1f mm\n' ...
            'top 3 mm / below  %6.1f / %.1f degC\n' ...
            'catheter          %6.1f degC  (cross-flow h = %.0f W/m^2K)'], ...
            R.eta*S.P, 100*R.eta, 100*R.P_gel/(R.eta*R.P_ref), Dtxt, T_protect, Ltxt, T_lesion, ...
            M.D_les*1e3, M.T_peak, x(M.ix_peak), z(M.iz_peak), M.T_top3, M.T_below3, ...
            M.T_cath, R.cath.h_c);
        drawnow limitrate;
    end

    function closeApp()
        try, stop(tmr); delete(tmr); catch, end
        delete(fig);
    end
end
