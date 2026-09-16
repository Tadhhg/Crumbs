function fig = gel_heating_gui(R)
%GEL_HEATING_GUI  Interactive gel heating, two tabs.
%
%   >> gel_heating_gui                 % model_defaults (the photographed setup)
%   >> gel_heating_gui(R)              % a solution from gel_heating_solve
%
%   Tab "Heating": sliders for generator power, water temperature, pre-cooling
%   time, exposure time and a time scrub; Play animates the exposure; Optimise
%   runs gel_heating_optimise; Baseline = 50 W, bath water, 332 s (the
%   photographed run).  Instant, because the fields are linear combinations
%   of the stored unit responses.
%
%   Tab "Constants & boundaries": every constant of model_defaults as a
%   number you can type (lengths in mm), the quantities derived from them
%   (King's k_e, film coefficients, grid size, time step, estimated solve
%   time), and "Apply & re-solve" -- the 3-D problem is solved again for the
%   new constants (~10 s at dh = 0.5 mm, ~5 min at 0.25 mm), cached, and the
%   Heating tab switches to it.  Empty = the model's own value (King's k_e,
%   the flat-plate h_top, the wall in series); h_top = Inf pins the surface.

if nargin < 1, R = gel_heating_solve(); end
T_protect = 45;  T_lesion = 50;
x = R.x*1e3;  z = R.z*1e3;

% ------------------------------------------------------------- STATE ----
S = struct('P', 50, 'T_w', R.T_bath0, 't_pre', 0, 't_exp', 332, 't', 0, 'playing', false);
t_max_total = R.t(end);

% ----------------------------------------------------------- WINDOW -----
fig = uifigure('Name', 'Gel heating V4', 'Position', [40 40 1460 800], 'Color', 'w');
glf  = uigridlayout(fig, [1 1], 'Padding', [0 0 0 0]);
tg   = uitabgroup(glf);
tabH = uitab(tg, 'Title', '  Heating  ');
tabC = uitab(tg, 'Title', '  Constants & boundaries  ');

% ============================================================ TAB 1 =====
gl  = uigridlayout(tabH, [1 2], 'ColumnWidth', {'1.55x', '1x'});
axMap = uiaxes(gl);  axMap.Layout.Column = 1;
right = uigridlayout(gl, [3 1], 'RowHeight', {'1x', '1x', 400}, 'Padding', [0 0 0 0]);
right.Layout.Column = 2;
axProf = uiaxes(right);  axHist = uiaxes(right);
ctrl = uigridlayout(right, [7 3], 'ColumnWidth', {150, '1x', 90}, ...
                    'RowHeight', {36 36 36 36 36 40 '1x'}, 'RowSpacing', 4);

sl = struct();
sl.P     = addslider(1, 'Generator power P_{set} (W)', [0 150], S.P);
sl.T_w   = addslider(2, 'Water temperature (\circC)', [2 45], S.T_w);
sl.t_pre = addslider(3, 'Pre-cool before power (s)', [0 300], S.t_pre);
sl.t_exp = addslider(4, 'Exposure time (s)', [10 t_max_total], S.t_exp);
sl.t     = addslider(5, 'Time since power on (s)', [0 S.t_exp], 0);
sl.P.MajorTicks = 0:25:150;         sl.T_w.MajorTicks = [2 5:5:45];
sl.t_pre.MajorTicks = 0:60:300;     sl.t_pre.MinorTicks = [];
sl.t_exp.MajorTicks = 60:60:t_max_total;  sl.t_exp.MinorTicks = [];
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
warning('off', 'MATLAB:contour:ConstantData');   % the t = 0 field is uniform
hold(axMap, 'on');
im = imagesc(axMap, x, z, zeros(numel(z), numel(x)));
set(axMap, 'YDir', 'normal');  axis(axMap, 'image');
colormap(axMap, turbo);  clim(axMap, [20 90]);
cb = colorbar(axMap);  cb.Label.String = 'T (\circC)';
[~, c45] = contour(axMap, x, z, zeros(numel(z), numel(x)), [T_protect T_protect], 'w--', 'LineWidth', 1.2);
[~, c50] = contour(axMap, x, z, zeros(numel(z), numel(x)), [T_lesion  T_lesion],  'w-',  'LineWidth', 1.6);
hD   = plot(axMap, x([1 end]), [0 0], 'c-', 'LineWidth', 1.2);
hPk  = plot(axMap, 0, 0, 'wo', 'MarkerSize', 8, 'LineWidth', 1.5);
hAnt = plot(axMap, [0 R.L_el]*1e3, [1 1]*R.h_ax*1e3, '-', 'Color', [.3 .3 .3], 'LineWidth', 4);
xlabel(axMap, 'x (mm)');  ylabel(axMap, 'z (mm)');
ttl = title(axMap, '');

% -------------------------------------------------------- PROFILE -------
hold(axProf, 'on');
hProf = plot(axProf, zeros(size(z)), z, 'r-', 'LineWidth', 1.8);
xline(axProf, T_protect, 'k--', sprintf('%g', T_protect), 'LabelOrientation', 'horizontal');
xline(axProf, T_lesion,  'k-',  sprintf('%g', T_lesion),  'LabelOrientation', 'horizontal');
hProfD = yline(axProf, 0, 'c-', 'protected layer', 'LineWidth', 1.2);
xlim(axProf, [0 100]);  grid(axProf, 'on');
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

% ============================================================ TAB 2 =====
%  field   label                                   unit     display scale   kind
spec = {
 'ANTENNA', '', '', 1, 'head';
 'f',        'frequency',                          'GHz',   1e9,   'num';
 'L',        'element length',                     'mm',    1e-3,  'num';
 'a',        'conductor radius',                   'mm',    1e-3,  'num';
 'b_ins',    'insulation outer radius (= a: bare)','mm',    1e-3,  'num';
 'epsr_ins', 'insulation permittivity',            '-',     1,     'num';
 'h_axis',   'wire axis above the gel surface',    'mm',    1e-3,  'num';
 'L_sheath', 'return-current arm on the braid',    'mm',    1e-3,  'num';
 'ke',       'current wavenumber (empty = King)',  '1/m',   1,     'cplx';
 'eta',      'delivered / set power',              '-',     1,     'num';
 'P_ref',    'reference generator setting',        'W',     1,     'num';
 'GEL', '', '', 1, 'head';
 'rho',      'density',                            'kg/m^3', 1,    'num';
 'c',        'specific heat',                      'J/(kg K)', 1,  'num';
 'k',        'thermal conductivity',               'W/(m K)', 1,   'num';
 'epsr_gel', 'permittivity, real part',            '-',     1,     'num';
 'epsi_gel', 'permittivity, loss factor',          '-',     1,     'num';
 'WATER ABOVE THE GEL', '', '', 1, 'head';
 'T_bath0',  'temperature the gel started at',     'degC',  1,     'num';
 'Q_circ',   'circulation flow',                   'L/min', 1/60000, 'num';
 'A_flow',   'flow cross-section over the gel',    'mm^2',  1e-6,  'num';
 'nu_w',     'kinematic viscosity',                'm^2/s', 1,     'num';
 'k_w',      'thermal conductivity',               'W/(m K)', 1,   'num';
 'Pr_w',     'Prandtl number',                     '-',     1,     'num';
 'epsr_w',   'permittivity (for King''s k_e)',     '-',     1,     'cplx';
 'h_top',    'h at the surface (empty = plate, Inf = pinned)', 'W/(m^2 K)', 1, 'opt';
 'TANK WALLS (sides, bottom)', '', '', 1, 'head';
 'T_room',   'room air',                           'degC',  1,     'num';
 'h_air',    'h outside the wall (convection + radiation)', 'W/(m^2 K)', 1, 'num';
 't_wall',   'wall thickness',                     'mm',    1e-3,  'num';
 'k_wall',   'wall conductivity',                  'W/(m K)', 1,   'num';
 'h_wall',   'h of the wall (empty = series, 0 = adiabatic)', 'W/(m^2 K)', 1, 'opt';
 'BLOCK AND GRID', '', '', 1, 'head';
 'x_min',    'block edge behind the feed',         'mm',    1e-3,  'num';
 'x_max',    'block edge beyond the tip',          'mm',    1e-3,  'num';
 'W_gel',    'block width across the antenna',     'mm',    1e-3,  'num';
 'H_gel',    'block depth',                        'mm',    1e-3,  'num';
 'dh',       'grid spacing (cost ~ dh^-5)',        'mm',    1e-3,  'num';
 't_end',    'length of the stored response',      's',     1,     'num';
 'dt_frame', 'stored frame interval',              's',     1,     'num';
};
nrow = size(spec, 1);
glC = uigridlayout(tabC, [1 2], 'ColumnWidth', {'1.3x', '1x'}, 'Padding', [12 12 12 12]);
left = uigridlayout(glC, [nrow 3], 'ColumnWidth', {330, 130, 90}, 'RowHeight', repmat({24}, 1, nrow), ...
                    'RowSpacing', 2, 'Scrollable', 'on');
left.Layout.Column = 1;
ef = struct();  efmeta = struct();
for i = 1:nrow
    [fld, lab, unit, scale, kind] = spec{i, :};
    if strcmp(kind, 'head')
        h = uilabel(left, 'Text', fld, 'FontWeight', 'bold', 'FontColor', [0.1 0.3 0.6]);
        h.Layout.Row = i;  h.Layout.Column = [1 3];
        continue
    end
    h = uilabel(left, 'Text', lab);  h.Layout.Row = i;  h.Layout.Column = 1;
    e = uieditfield(left, 'text', 'HorizontalAlignment', 'right', 'FontName', 'Consolas', ...
                    'ValueChangedFcn', @(~,~) refreshDerived());
    e.Layout.Row = i;  e.Layout.Column = 2;
    u = uilabel(left, 'Text', unit, 'FontColor', [0.4 0.4 0.4]);  u.Layout.Row = i;  u.Layout.Column = 3;
    ef.(fld) = e;  efmeta.(fld) = struct('scale', scale, 'kind', kind);
end

rightC = uigridlayout(glC, [4 1], 'RowHeight', {44, '1x', 30, 60}, 'Padding', [0 0 0 0]);
rightC.Layout.Column = 2;
bRow = uigridlayout(rightC, [1 3], 'Padding', [0 0 0 0]);
bRow.Layout.Row = 1;
uibutton(bRow, 'Text', 'Apply & re-solve', 'FontWeight', 'bold', 'ButtonPushedFcn', @(~,~) applySolve());
uibutton(bRow, 'Text', 'Photo defaults', 'ButtonPushedFcn', @(~,~) loadFields(model_defaults()));
uibutton(bRow, 'Text', 'Current solution', 'ButtonPushedFcn', @(~,~) loadFields(R.P));
derived = uitextarea(rightC, 'Editable', 'off', 'FontName', 'Consolas', 'FontSize', 11);
derived.Layout.Row = 2;
status = uilabel(rightC, 'Text', '', 'FontColor', [0.6 0 0], 'WordWrap', 'on');
status.Layout.Row = 3;
note = uilabel(rightC, 'WordWrap', 'on', 'FontColor', [0.35 0.35 0.35], 'VerticalAlignment', 'top', 'Text', ...
    ['Lengths in mm, flow in L/min, frequency in GHz.  Empty = the model''s own value.  Every distinct set of ' ...
     'constants is solved once and cached in solutions/.  Sources and the photograph measurements: model_defaults.m, ' ...
     '../Gel_model_V3/gel_heating_notes.md.']);
note.Layout.Row = 4;

% ----------------------------------------------------------- TIMER ------
tmr = timer('ExecutionMode', 'fixedRate', 'Period', 0.05, 'BusyMode', 'drop', ...
            'TimerFcn', @(~,~) tick());
fig.CloseRequestFcn = @(~,~) closeApp();
fig.UserData = struct('getState', @getState, 'setScenario', @setScenario, 'getR', @() R);

loadFields(R.P);
rebuild();
if nargout == 0, clear fig; end

% =========================================================================
%  TAB 1
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
                sl.t_exp.Limits = [10 max(11, t_max_total - S.t_pre)];  sl.t_exp.Value = S.t_exp;
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

    function setScenario(P_set, T_w, t_pre, t_exp, t)
        S.P = P_set;  S.T_w = T_w;  S.t_pre = t_pre;  S.t_exp = t_exp;  S.t = t;
        sl.P.Value = min(max(P_set, sl.P.Limits(1)), sl.P.Limits(2));
        sl.T_w.Value = min(max(T_w, sl.T_w.Limits(1)), sl.T_w.Limits(2));
        sl.t_pre.Value = t_pre;
        sl.t_exp.Limits = [10 max(11, t_max_total - t_pre)];  sl.t_exp.Value = t_exp;
        sl.t.Limits = [0 t_exp];  sl.t.Value = t;  setTimeTicks();
        refreshHistory();  update();
    end

    function setTime(t)
        S.t = t;  sl.t.Value = t;  update();
    end

    function setTimeTicks()
        steps = [5 10 20 30 60 120];
        st = steps(find(S.t_exp./steps <= 12, 1));
        if isempty(st), st = steps(end); end
        sl.t.MajorTicks = 0:st:S.t_exp;  sl.t.MinorTicks = [];
    end

    function baseline(), setScenario(50, R.T_bath0, 0, min(332, t_max_total), 0); end

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

    function rebuild()
        % the solution R changed (new constants): re-point every R-dependent object
        x = R.x*1e3;  z = R.z*1e3;  t_max_total = R.t(end);
        im.XData = x;  im.YData = z;  im.CData = zeros(numel(z), numel(x));
        for hc = [c45 c50]
            hc.XData = x;  hc.YData = z;  hc.ZData = zeros(numel(z), numel(x));
        end
        hD.XData = x([1 end]);
        hAnt.XData = [0 R.L_el]*1e3;  hAnt.YData = [1 1]*R.h_ax*1e3;
        xlim(axMap, [x(1) x(end)]);  ylim(axMap, [z(1) max(2, R.h_ax*1e3 + 1)]);
        hProf.YData = z;  ylim(axProf, [z(1) 0]);
        S.t_pre = min(S.t_pre, t_max_total - 10);  sl.t_pre.Value = S.t_pre;
        S.t_exp = min(S.t_exp, t_max_total - S.t_pre);
        sl.t_exp.Limits = [10 max(11, t_max_total - S.t_pre)];  sl.t_exp.Value = S.t_exp;
        sl.t_exp.MajorTicks = 60:60:t_max_total;
        S.t = min(S.t, S.t_exp);  sl.t.Limits = [0 S.t_exp];  sl.t.Value = S.t;  setTimeTicks();
        refreshHistory();  update();
    end

    function refreshHistory()
        tt = 0:2:S.t_exp;  top = zeros(size(tt));  bel = top;
        for kk = 1:numel(tt)
            [~, M] = gel_heating_eval(R, S.P, S.T_w, S.t_pre, tt(kk), T_protect, T_lesion);
            top(kk) = M.T_top3;  bel(kk) = M.T_below3;
        end
        set(hTop, 'XData', tt, 'YData', top);  set(hBel, 'XData', tt, 'YData', bel);
        Tend = gel_heating_eval(R, S.P, S.T_w, S.t_pre, S.t_exp);
        clim(axMap, [min(20, floor(min(S.T_w, R.T_bath0))) max(ceil(max(Tend(:))) + 2, 50)]);
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
                              ', pre-cool %.0f s      peak %.1f \\circC at (%.1f, %.1f) mm      dashed %g, solid %g \\circC'], ...
                              S.t, S.P, R.eta*S.P, S.T_w, S.t_pre, M.T_peak, x(M.ix_peak), z(M.iz_peak), T_protect, T_lesion);
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
            'antenna           %g mm element, axis %.2f mm above the gel, k_e = %.0f%+.0fj 1/m, arm %g mm\n' ...
            'delivered power   %6.1f W  (%.0f %% of setting, %.0f %% of it into the gel)\n' ...
            'protected layer   %s  (T <= %g degC, min over x)\n' ...
            'max T below it    %s degC     %g degC first reached at %.2f mm\n' ...
            'peak              %6.1f degC at x = %.1f mm, z = %.1f mm\n' ...
            'top 3 mm / below  %6.1f / %.1f degC\n' ...
            'catheter          %6.1f degC  (cross-flow h = %.0f W/m^2K)     grid %g mm, %d nodes'], ...
            R.L_el*1e3, R.h_ax*1e3, real(R.ke), imag(R.ke), R.P.L_sheath*1e3, ...
            R.eta*S.P, 100*R.eta, 100*R.P_gel/(R.eta*R.P_ref), Dtxt, T_protect, Ltxt, T_lesion, ...
            M.D_les*1e3, M.T_peak, x(M.ix_peak), z(M.iz_peak), M.T_top3, M.T_below3, ...
            M.T_cath, R.cath.h_c, R.dh*1e3, R.D.N);
        drawnow limitrate;
    end

% =========================================================================
%  TAB 2
    function loadFields(Pin)
        fn = fieldnames(ef);
        for i = 1:numel(fn)
            m = efmeta.(fn{i});  v = Pin.(fn{i});
            if isempty(v),        ef.(fn{i}).Value = '';
            elseif isinf(v),      ef.(fn{i}).Value = 'Inf';
            elseif strcmp(m.kind, 'cplx'), ef.(fn{i}).Value = sprintf('%.6g%+.6gj', real(v), imag(v));
            else,                 ef.(fn{i}).Value = sprintf('%.10g', v/m.scale);
            end
        end
        refreshDerived();
    end

    function Pout = readFields()
        Pout = model_defaults();                       % field order and any field not on the tab
        fn = fieldnames(ef);
        for i = 1:numel(fn)
            m = efmeta.(fn{i});  s = strtrim(ef.(fn{i}).Value);
            if isempty(s)
                if strcmp(m.kind, 'num')
                    error('gel_heating_gui:empty', '"%s" cannot be empty', fn{i});
                end
                Pout.(fn{i}) = [];
                continue
            end
            v = str2double(s);                         % accepts 172-26j, 1e-3, Inf
            if ~isscalar(v) || isnan(v)
                error('gel_heating_gui:parse', '"%s": cannot read "%s"', fn{i}, s);
            end
            if strcmp(m.kind, 'cplx'), Pout.(fn{i}) = v;
            else,                      Pout.(fn{i}) = v*m.scale;
            end
        end
        if Pout.b_ins < Pout.a,        error('gel_heating_gui:geom', 'insulation radius must be >= conductor radius'); end
        if Pout.h_axis <= Pout.b_ins,  error('gel_heating_gui:geom', 'the wire axis must be above the gel by more than its outer radius'); end
        if Pout.x_min >= 0 || Pout.x_max <= Pout.L, error('gel_heating_gui:geom', 'the block must contain the element (x_min < 0 < L < x_max)'); end
        if Pout.dh <= 0 || Pout.dt_frame <= 0 || Pout.t_end < 2*Pout.dt_frame, error('gel_heating_gui:grid', 'dh, dt_frame > 0 and t_end >= 2 dt_frame'); end
    end

    function refreshDerived()
        try
            Pn = readFields();  D = model_derived(Pn);
            status.Text = '';
            same = strcmp(D.key, R.D.key);
            if same, tag = 'these are the constants of the current solution';
            else,    tag = 'CHANGED -- press "Apply & re-solve"'; end
            derived.Value = sprintf([ ...
                'DERIVED  (%s)\n\n' ...
                'current wave   k_e = %.1f %+.1fj 1/m   [%s]\n' ...
                '               lambda_e = %.1f mm,  element = %.2f lambda_e\n' ...
                'gel            k_c = %.1f %+.1fj 1/m,  lambda = %.1f mm,  power 1/e depth %.1f mm\n' ...
                '               sigma_eff = %.2f S/m,  alpha = %.3g m^2/s\n' ...
                'top face       h_top = %s W/(m^2 K)  [%s]\n' ...
                '               U = %.3f m/s, Re_L = %.0f, plate value %.0f, Bi = %s\n' ...
                'walls          h_wall = %.2f W/(m^2 K)  [%s],  Bi = %.2f\n' ...
                'catheter gap   %.2f mm of water between insulation and gel\n\n' ...
                'grid           %d x %d x %d = %d nodes,  dh = %g mm\n' ...
                'time step      dt = %.4f s,  %d steps to %g s\n' ...
                'estimated solve time  ~%.0f s%s\n'], ...
                tag, real(D.ke), imag(D.ke), D.ke_source, D.lambda_e*1e3, D.L_over_lambda, ...
                real(D.kc), imag(D.kc), D.lambda_gel*1e3, D.d_1e_power*1e3, D.sig_eff, D.alpha, ...
                num2str(D.h_top, '%.0f'), D.h_top_source, D.U_w, D.Re_L, D.h_top_corr, num2str(D.Bi_top, '%.1f'), ...
                D.h_wall, D.h_wall_source, D.Bi_wall, D.gap*1e3, ...
                D.nx, D.ny, D.nz, D.N, Pn.dh*1e3, D.dt, D.nsteps, Pn.t_end, D.t_est, ...
                tern(D.t_est > 120, '  -- consider a coarser dh while exploring', ''));
        catch err
            status.Text = err.message;
        end
    end

    function applySolve()
        try
            Pn = readFields();  D = model_derived(Pn);
        catch err
            uialert(fig, err.message, 'Cannot read the constants');  return
        end
        if strcmp(D.key, R.D.key)
            status.Text = 'nothing changed -- the current solution already uses these constants';  return
        end
        if S.playing, togglePlay(); end
        dlg = uiprogressdlg(fig, 'Title', 'Solving', 'Indeterminate', 'on', 'Message', ...
              sprintf('%d nodes, %d time steps, about %.0f s ...', D.N, D.nsteps, D.t_est));
        try
            Rn = gel_heating_solve(Pn);
            close(dlg);
        catch err
            close(dlg);  uialert(fig, err.message, 'Solve failed');  return
        end
        R = Rn;
        S.T_w = R.T_bath0;  sl.T_w.Value = min(max(S.T_w, sl.T_w.Limits(1)), sl.T_w.Limits(2));
        rebuild();
        loadFields(R.P);
        status.Text = sprintf('solved in %.0f s: %.2f W of %.1f W into the gel; cached as %s', ...
                              R.t_solve, R.P_gel, R.eta*R.P_ref, R.file);
        tg.SelectedTab = tabH;
    end

    function closeApp()
        try, stop(tmr); delete(tmr); catch, end
        delete(fig);
    end
end

function s = tern(c, a, b)
if c, s = a; else, s = b; end
end
