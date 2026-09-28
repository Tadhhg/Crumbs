function varargout = v5_gui(preset)
%V5_GUI  Coupled dielectric + catheter heating of the thermochromic gel.
%
%   >> v5_gui              % setup table (60 W, 27 mm, 2 mm, 37 degC water)
%   >> v5_gui('photo')     % the photographed run (50 W, 12 mm, 1.1 mm, 25 degC)
%
%   Tab "Model"
%     left    starting conditions and experimental settings (sliders + typed
%             values), the pigment band, the view, and the resolution with
%             the solve runtime.  A circled-arrow mark (re-solve) means the
%             stored responses are recomputed; everything else is instant.
%     centre  the antenna plane coloured like the photographs (or in degC,
%             rise, dielectric / diffusion parts, SAR); click to set a probe;
%             time slider and Play; history of catheter, peak, surface, probe.
%     right   the photograph nearest the current time (50 W run), the depth
%             profile through the hottest column split into the dielectric
%             and diffusion parts, and a readout of the derived quantities.
%   Tab "Advanced constants"
%     every other constant of v5_defaults as a typed value, marked instant or
%     re-solve, with the derived quantities and the runtime estimate.
%
%   Model: ../Catheter_heating/Convection model.docx; E-field:
%   ../Insulated_monopole_Erms; solver v5_solve (exact-in-time modal), v5_eval.

if nargin < 1, preset = 'table'; end
here = fileparts(mfilename('fullpath'));
addpath(here);                            % callbacks need v5_*.m after run() restores the folder

P  = v5_defaults(preset);                 % current constants (instant + solve)
R  = [];  O = [];                         % solution and its evaluation
cache = containers.Map();  cacheOrder = {};
S = struct('n', 1, 'playing', false, 'probe', [NaN NaN], 'view', '', 'cmap', turbo(256), 'lim', [0 1], 'nEval', 0, 'envTitle', '', 'envLim', []);
rate = struct('cE', 2.7e-8, 'cF', 1.6e-8);  % s per (x-station * rho * segment), s per (node * frame); refit after each solve
photos = readPhotos();

% ============================================================= WINDOW ====
fig = uifigure('Name', 'V5 coupled heating: thermochromic gel', 'Position', [20 40 1560 880], 'Color', 'w');
g0  = uigridlayout(fig, [1 1], 'Padding', [0 0 0 0]);
tg  = uitabgroup(g0);
tabM = uitab(tg, 'Title', '  Model  ');
tabA = uitab(tg, 'Title', '  Advanced constants  ');

gM = uigridlayout(tabM, [1 3], 'ColumnWidth', {372, '1x', 400}, 'Padding', [6 6 6 6], 'ColumnSpacing', 8);

% ---------------------------------------------------------- CONTROLS -----
cp = uigridlayout(gM, [34 3], 'ColumnWidth', {150, '1x', 62}, 'RowSpacing', 2, 'Padding', [2 2 2 2], ...
                  'RowHeight', repmat({24}, 1, 34), 'Scrollable', 'on');
cp.Layout.Column = 1;
row = 0;
hdr('Preset');
row = row + 1;
ddPreset = uidropdown(cp, 'Items', {'Setup table (60 W, 27 mm)', 'Photographed run (50 W, 12 mm)'}, ...
                      'ItemsData', {'table', 'photo'}, 'Value', preset);
place(ddPreset, row, [1 2]);
b = uibutton(cp, 'Text', 'Load', 'ButtonPushedFcn', @(~,~) loadPreset(ddPreset.Value));  place(b, row, 3);

hdr('Starting conditions (setup table)');
C = struct();
C.P_set   = ctl('Input power (W)',           'P_set',     [0 150],  1,    'instant');
C.L       = ctl('Antenna length (mm) ↻',     'L',         [4 40],   1e-3, 'solve');
C.d       = ctl('Antenna–gel dist. (mm) ↻', 'd',         [0.5 8],  1e-3, 'solve');
hdr('Experiment');
C.P_dropoff = ctl('Power drop-off (lost)',   'P_dropoff', [0 0.98], 1,    'instant');
C.T_b     = ctl('Water temp. (°C)',          'T_b',       [10 50],  1,    'instant');
C.T_0     = ctl('Initial gel temp. (°C)',    'T_0',       [10 50],  1,    'instant');
C.flow    = ctl('Circulation (L/min)',       'flow',      [0.2 10], 1,    'instant');
hdr('Power envelope (duty cycle)');
C.t_on    = ctl('Treatment time (s)',        't_on',      [0 P.t_end], 1, 'instant');
C.duty    = ctl('Duty cycle (%)',            'duty',      [1 100],  0.01, 'instant');
C.T_pulse = ctl('Pulse period (s)',          'T_pulse',   [0.5 120], 1,   'instant');
C.t_burst = ctl('Burst length (s)',          't_burst',   [2 600],  1,    'instant');
C.t_rest  = ctl('Rest between bursts (s)',   't_rest',    [0 300],  1,    'instant');
hdr('Pigment (colouring)');
C.T_tlc   = ctl('Colour start (°C)',         'T_tlc',     [20 70],  1,    'display');
C.bw_tlc  = ctl('Band width (K)',            'bw_tlc',    [1 40],   1,    'display');
hdr('View');
ddView = uidropdown(cp, 'Items', {'Pigment (as photographed)', 'Pigment (vivid)', 'Temperature (°C)', 'Rise above water (K)', ...
                    'Dielectric part θ_d (K)', 'Diffusion part θ_s (K)', 'SAR (W/kg)'}, ...
                    'ValueChangedFcn', @(s,~) setView(s.Value));
row = row + 1;  place(ddView, row, [1 3]);
row = row + 1;
cbIso = uicheckbox(cp, 'Text', 'isotherm at colour start', 'Value', true, 'ValueChangedFcn', @(~,~) drawFrame());
place(cbIso, row, [1 3]);
hdr('Resolution and runtime');
ddDim = uidropdown(cp, 'Items', {'3-D block (accurate)', '2-D slice (fast, overestimates)'}, 'ItemsData', {3, 2}, ...
                   'Value', P.dim, 'ValueChangedFcn', @(~,~) resolutionChanged());
row = row + 1;  place(ddDim, row, [1 3]);
row = row + 1;
lb = uilabel(cp, 'Text', 'Grid spacing Δh');  place(lb, row, 1);
dhList = [1 0.75 0.5 0.35 0.25 0.2 0.15]*1e-3;
ddDh = uidropdown(cp, 'Items', compose('%.2f mm', dhList*1e3), 'ItemsData', num2cell(dhList), ...
                  'Value', P.dh, 'ValueChangedFcn', @(~,~) resolutionChanged());
place(ddDh, row, 2);
bSolve = uibutton(cp, 'Text', 'Solve', 'FontWeight', 'bold', 'ButtonPushedFcn', @(~,~) solveNow(true));
place(bSolve, row, 3);
row = row + 1;
lblRun = uilabel(cp, 'Text', '', 'FontName', 'Consolas', 'FontSize', 10.5, 'WordWrap', 'on', 'VerticalAlignment', 'top');
place(lblRun, row, [1 3]);  cp.RowHeight{row} = 92;
row = row + 1;
bExp = uibutton(cp, 'Text', 'Export frame (PNG + settings)', 'ButtonPushedFcn', @(~,~) exportFrame());
place(bExp, row, [1 3]);
cp.RowHeight = cp.RowHeight(1:row);

% ------------------------------------------------------------ CENTRE -----
gc = uigridlayout(gM, [4 1], 'RowHeight', {22, '1.7x', 34, '1x'}, 'Padding', [0 0 0 0], 'RowSpacing', 4);
gc.Layout.Column = 2;
ttlMap = uilabel(gc, 'Text', '', 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Interpreter', 'tex');
ttlMap.Layout.Row = 1;         % labels, not axes titles: a changing title re-lays the axes out
axMap = uiaxes(gc);  axMap.Layout.Row = 2;
tr = uigridlayout(gc, [1 5], 'ColumnWidth', {70, 60, '1x', 70, 110}, 'Padding', [0 0 0 0]);
tr.Layout.Row = 3;
bPlay = uibutton(tr, 'Text', '▶ Play', 'ButtonPushedFcn', @(~,~) togglePlay());
uibutton(tr, 'Text', '⟲ 0 s', 'ButtonPushedFcn', @(~,~) setFrame(1));
slT = uislider(tr, 'Limits', [0 P.t_end], 'Value', 0, 'ValueChangingFcn', @(~,e) setTime(e.Value), ...
               'ValueChangedFcn', @(s,~) setTime(s.Value));
lblT = uilabel(tr, 'Text', 't = 0 s', 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
ddSpeed = uidropdown(tr, 'Items', {'10 s / s', '30 s / s', '60 s / s', '120 s / s'}, 'ItemsData', {10, 30, 60, 120}, 'Value', 30);
axHist = uiaxes(gc);  axHist.Layout.Row = 4;

% ------------------------------------------------------------- RIGHT -----
gr = uigridlayout(gM, [6 1], 'RowHeight', {34, '1x', 0, 20, '1x', 250}, 'Padding', [0 0 0 0], 'RowSpacing', 2);
gr.Layout.Column = 3;
gph = uigridlayout(gr, [1 2], 'ColumnWidth', {'1x', 128}, 'Padding', [0 0 0 0]);  gph.Layout.Row = 1;
ttlPh = uilabel(gph, 'Text', '', 'HorizontalAlignment', 'center', 'WordWrap', 'on', 'FontSize', 11);
cbPhoto = uicheckbox(gph, 'Text', 'show photographs', 'Value', false, 'ValueChangedFcn', @(~,~) togglePhoto());
axEnv = uiaxes(gr);  axEnv.Layout.Row = 2;     % power envelope (default) ...
axPh = uiaxes(gr);   axPh.Layout.Row = 3;      % ... or the photographs: the hidden one's row is
axPh.Visible = 'off';                          % collapsed, so it cannot sit over the other and eat the mouse
ttlProf = uilabel(gr, 'Text', '', 'HorizontalAlignment', 'center', 'FontSize', 11);
ttlProf.Layout.Row = 4;
axProf = uiaxes(gr);  axProf.Layout.Row = 5;
txtOut = uitextarea(gr, 'Editable', 'off', 'FontName', 'Consolas', 'FontSize', 10.5);
txtOut.Layout.Row = 6;

% --------------------------------------------------------- MAP OBJECTS ----
warning('off', 'MATLAB:contour:ConstantData');
hold(axMap, 'on');  axMap.Toolbar.Visible = 'off';  disableDefaultInteractivity(axMap);
hWater = patch(axMap, 'XData', [0 1 1 0], 'YData', [0 0 1 1], 'FaceColor', [0.44 0.67 0.75], 'EdgeColor', 'none');   % photographed water
hImg = image(axMap, 'XData', [0 1], 'YData', [0 -1], 'CData', zeros(2, 2, 3), 'ButtonDownFcn', @(~,e) mapClick(e));
set(axMap, 'YDir', 'normal');
[~, hIso] = contour(axMap, [0 1], [0 -1], zeros(2), [1 1], '--', 'Color', [1 1 1], 'LineWidth', 1.1);
hCable = patch(axMap, 'XData', [0 1 1 0], 'YData', [0 0 1 1], 'FaceColor', [0.86 0.58 0.48], 'EdgeColor', [0.6 0.35 0.3]);
hFeed  = patch(axMap, 'XData', [0 1 1 0], 'YData', [0 0 1 1], 'FaceColor', [0.88 0.88 0.82], 'EdgeColor', [0.5 0.5 0.5]);
hIns   = patch(axMap, 'XData', [0 1 1 0], 'YData', [0 0 1 1], 'FaceColor', [0.86 0.93 0.97], 'FaceAlpha', 0.7, 'EdgeColor', [0.55 0.65 0.7]);
hCond  = patch(axMap, 'XData', [0 1 1 0], 'YData', [0 0 1 1], 'FaceColor', [0.97 0.97 0.99], 'EdgeColor', [0.6 0.6 0.65]);
hPk    = plot(axMap, NaN, NaN, 'o', 'Color', [0.1 0.1 0.1], 'MarkerSize', 9, 'LineWidth', 1.4);
hProbe = plot(axMap, NaN, NaN, '+', 'Color', [0.8 0 0.4], 'MarkerSize', 12, 'LineWidth', 2);
axis(axMap, 'image');  box(axMap, 'on');  axMap.Layer = 'top';
axMap.Color = [0.80 0.83 0.85];              % outside the gel block: tank
set([hWater hIso hCable hFeed hIns hCond hPk hProbe], 'HitTest', 'off', 'PickableParts', 'none');
xlabel(axMap, 'x along the antenna (mm), feed at 0');  ylabel(axMap, 'z (mm), gel surface at 0');
cbMap = colorbar(axMap);

% ------------------------------------------------------ OTHER AXES --------
hold(axHist, 'on');  grid(axHist, 'on');  axHist.Toolbar.Visible = 'off';
hBand  = patch(axHist, 'XData', [0 1 1 0], 'YData', [0 0 1 1], 'FaceColor', [0.55 0.85 0.75], 'FaceAlpha', 0.18, ...
               'EdgeColor', 'none', 'DisplayName', 'pigment band');
hCath  = plot(axHist, NaN, NaN, '-', 'Color', [0.85 0.45 0.1], 'LineWidth', 1.6, 'DisplayName', 'catheter tip');
hPeakL = plot(axHist, NaN, NaN, '-', 'Color', [0.8 0.1 0.1], 'LineWidth', 1.8, 'DisplayName', 'peak in the gel plane');
hSurf  = plot(axHist, NaN, NaN, '-', 'Color', [0.1 0.4 0.8], 'LineWidth', 1.4, 'DisplayName', 'max at the gel surface');
hPrH   = plot(axHist, NaN, NaN, '-', 'Color', [0.8 0 0.4], 'LineWidth', 1.4, 'DisplayName', 'probe');
hNow   = xline(axHist, 0, 'k:', 'LineWidth', 1.2, 'HandleVisibility', 'off');
hOff   = xline(axHist, P.t_on, '--', 'treatment end', 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off', ...
               'LabelVerticalAlignment', 'bottom');
xlabel(axHist, 'time since power on (s)');  ylabel(axHist, 'T (°C)');
legend(axHist, 'Location', 'northoutside', 'NumColumns', 5, 'FontSize', 9);

hold(axPh, 'on');  axPh.Toolbar.Visible = 'off';  disableDefaultInteractivity(axPh);
hPh = image(axPh, 'XData', [-8 26], 'YData', [4 -16], 'CData', zeros(2, 2, 3));
set(axPh, 'YDir', 'normal');  axis(axPh, 'image');  xlim(axPh, [-8 26]);  ylim(axPh, [-16 4]);
xlabel(axPh, 'x (mm)');  ylabel(axPh, 'z (mm)');
hPh.Visible = 'off';

hold(axEnv, 'on');  grid(axEnv, 'on');         % pan (drag), zoom (scroll) and the toolbar stay on
hEnvF = area(axEnv, [0 1], [0 0], 'FaceColor', [0.85 0.35 0.2], 'FaceAlpha', 0.35, 'EdgeColor', [0.7 0.2 0.1], ...
             'DisplayName', 'at the antenna');
hEnvA = plot(axEnv, NaN, NaN, '-', 'Color', [0.2 0.2 0.2], 'LineWidth', 1.2, 'DisplayName', 'frame average (used)');
hEnvM = yline(axEnv, 0, ':', '', 'Color', [0.1 0.4 0.8], 'LineWidth', 1.4, 'DisplayName', 'mean over treatment', ...
              'LabelHorizontalAlignment', 'left');
hEnvN = xline(axEnv, 0, 'k:', 'LineWidth', 1.2, 'HandleVisibility', 'off');
xlabel(axEnv, 'time since power on (s)');  ylabel(axEnv, 'P at the antenna (W)');
legend(axEnv, 'Location', 'north', 'NumColumns', 3, 'FontSize', 8);

hold(axProf, 'on');  grid(axProf, 'on');  axProf.Toolbar.Visible = 'off';
hPrT = plot(axProf, NaN, NaN, 'k-', 'LineWidth', 1.8, 'DisplayName', 'T');
hPrD = plot(axProf, NaN, NaN, '-', 'Color', [0.8 0.1 0.1], 'LineWidth', 1.3, 'DisplayName', 'T_b + θ_d (dielectric)');
hPrS = plot(axProf, NaN, NaN, '-', 'Color', [0.1 0.4 0.8], 'LineWidth', 1.3, 'DisplayName', 'T_b + θ_s (catheter diffusion)');
hPrX = xline(axProf, 40, '--', 'colour start', 'Color', [0.2 0.6 0.5], 'HandleVisibility', 'off');
xlabel(axProf, 'T (°C)');  ylabel(axProf, 'z (mm)');
legend(axProf, 'Location', 'southeast', 'FontSize', 8);

% ======================================================= ADVANCED TAB ====
%  field       label                                          unit        scale   group
spec = {
 'POWER',     '', '', 1, '';
 'P_refl',    'reflected power',                              'W',        1,     'instant';
 'gamma',     'share of the lost power heating the tip, γ',   '-',        1,     'instant';
 'ANTENNA AND CABLE', '', '', 1, '';
 'f',         'frequency',                                    'GHz',      1e9,   'solve';
 'b',         'insulated element outer radius (r0)',          'mm',       1e-3,  'solve';
 'a',         'conductor radius (ASSUMED b/3.35)',            'mm',       1e-3,  'solve';
 'eps_ins',   'insulation permittivity (PTFE)',               '-',        1,     'solve';
 'l_f',       'exposed feed block behind the feed',           'mm',       1e-3,  'solve';
 'r_cable',   'coax outer radius (drawing only)',             'mm',       1e-3,  'display';
 'GEL (1 % Phytagel)', '', '', 1, '';
 'eps_g1',    'permittivity, real ε''',                       '-',        1,     'solve';
 'eps_g2',    'permittivity, loss ε'''' (σ = ωε0ε'''')',      '-',        1,     'solve';
 'rho',       'density ρ',                                    'kg/m^3',   1,     'solve';
 'c',         'specific heat c',                              'J/(kg K)', 1,     'solve';
 'k',         'thermal conductivity k',                       'W/(m K)',  1,     'solve';
 'WATER', '', '', 1, '';
 'eps_w1',    'permittivity, real',                           '-',        1,     'solve';
 'eps_w2',    'permittivity, loss',                           '-',        1,     'solve';
 'BOUNDARIES', '', '', 1, '';
 'h_s',       'water-to-gel surface h_s (Inf = pinned)',      'W/(m^2 K)',1,     'solve';
 'h_w',       'tank walls to room h_w (0 = adiabatic)',       'W/(m^2 K)',1,     'solve';
 'T_room',    'room air',                                     'degC',     1,     'instant';
 'SURFACE TRANSFER  T_diff = T_b + α(T_cath − T_b)', '', '', 1, '';
 'beta',      'β in α = min(1, β/(d√f))',                     'm (L/min)^0.5', 1, 'instant';
 'L_c',       'heated catheter length (NaN = element)',       'mm',       1e-3,  'solve';
 'CATHETER TIP (lumped)', '', '', 1, '';
 'm_cat',     'heated tip mass m',                            'g',        1e-3,  'instant';
 'c_cat',     'tip specific heat',                            'J/(kg K)', 1,     'instant';
 'k_cat',     'tip-to-water loss k_c',                        'W/K',      1,     'instant';
 'GEL BLOCK', '', '', 1, '';
 'gel_L',     'length along the antenna',                     'mm',       1e-3,  'solve';
 'gel_H',     'height (depth)',                               'mm',       1e-3,  'solve';
 'gel_W',     'width across (ASSUMED)',                       'mm',       1e-3,  'solve';
 'x_start',   'gel edge behind the feed (feed at 0)',         'mm',       1e-3,  'solve';
 'NUMERICS', '', '', 1, '';
 't_end',     'length of the stored run',                     's',        1,     'solve';
 'dt_frame',  'frame interval (exact in time)',               's',        1,     'solve';
};
nsp = size(spec, 1);
gA = uigridlayout(tabA, [1 2], 'ColumnWidth', {'1.25x', '1x'}, 'Padding', [12 12 12 12]);
la = uigridlayout(gA, [nsp 4], 'ColumnWidth', {320, 120, 95, 70}, 'RowHeight', repmat({23}, 1, nsp), ...
                  'RowSpacing', 2, 'Scrollable', 'on');
la.Layout.Column = 1;
EF = struct();
for ii = 1:nsp
    [fld, lab, unit, ~, grp] = spec{ii, :};
    if isempty(grp)
        h = uilabel(la, 'Text', fld, 'FontWeight', 'bold', 'FontColor', [0.1 0.3 0.6]);
        h.Layout.Row = ii;  h.Layout.Column = [1 4];
        continue
    end
    h = uilabel(la, 'Text', lab);  h.Layout.Row = ii;  h.Layout.Column = 1;
    e = uieditfield(la, 'text', 'HorizontalAlignment', 'right', 'FontName', 'Consolas', ...
                    'ValueChangedFcn', @(~,~) refreshDerived());
    e.Layout.Row = ii;  e.Layout.Column = 2;
    u = uilabel(la, 'Text', unit, 'FontColor', [0.4 0.4 0.4]);  u.Layout.Row = ii;  u.Layout.Column = 3;
    switch grp
        case 'solve',   tagtxt = 're-solve';  col = [0.7 0.3 0];
        case 'instant', tagtxt = 'instant';   col = [0 0.5 0.2];
        otherwise,      tagtxt = 'display';   col = [0.4 0.4 0.4];
    end
    tgl = uilabel(la, 'Text', tagtxt, 'FontColor', col, 'FontSize', 11);  tgl.Layout.Row = ii;  tgl.Layout.Column = 4;
    EF.(fld) = struct('h', e, 'scale', spec{ii, 4}, 'group', grp);
end
ra = uigridlayout(gA, [4 1], 'RowHeight', {40, '1x', 36, 90}, 'Padding', [0 0 0 0]);
ra.Layout.Column = 2;
bra = uigridlayout(ra, [1 3], 'Padding', [0 0 0 0]);  bra.Layout.Row = 1;
uibutton(bra, 'Text', 'Apply', 'FontWeight', 'bold', 'ButtonPushedFcn', @(~,~) applyAdvanced());
uibutton(bra, 'Text', 'Reset tab to current', 'ButtonPushedFcn', @(~,~) loadAdvanced());
uibutton(bra, 'Text', 'Reset to preset', 'ButtonPushedFcn', @(~,~) loadPreset(ddPreset.Value));
txtDer = uitextarea(ra, 'Editable', 'off', 'FontName', 'Consolas', 'FontSize', 10.5);  txtDer.Layout.Row = 2;
lblAdv = uilabel(ra, 'Text', '', 'FontColor', [0.6 0 0], 'WordWrap', 'on');  lblAdv.Layout.Row = 3;
nt = uilabel(ra, 'WordWrap', 'on', 'VerticalAlignment', 'top', 'FontColor', [0.35 0.35 0.35], 'Text', ...
    ['"instant" constants only rescale the stored responses; "re-solve" constants recompute them ' ...
     '(the Model tab shows the time).  Values marked ASSUMED or guessed are first guesses from ' ...
     'Convection model.docx s6 and v5_defaults.m, not measurements.  Lengths in mm, mass in g.']);
nt.Layout.Row = 4;

% ============================================================== START ====
tmr = timer('ExecutionMode', 'fixedRate', 'Period', 0.05, 'BusyMode', 'drop', 'TimerFcn', @(~,~) tick());
fig.CloseRequestFcn = @(~,~) closeApp();
fig.UserData = struct('getP', @getP, 'getR', @getR, 'getO', @getO, 'setFrame', @setFrame, ...
                      'setParam', @setParam, 'solve', @() solveNow(true), 'setView', @setView);
syncControls();  loadAdvanced();
solveNow(true);
setFrame(nearestFrame(min(332, P.t_on)));
if nargout > 0, varargout{1} = fig; end   % fig must stay alive: the callbacks share it

% =========================================================================
%  CONTROL HELPERS
    function hdr(txt)
        row = row + 1;                   % row = the last row used
        h = uilabel(cp, 'Text', txt, 'FontWeight', 'bold', 'FontColor', [0.1 0.3 0.6]);
        place(h, row, [1 3]);
    end

    function place(h, r, c)
        h.Layout.Row = r;  h.Layout.Column = c;
    end

    function s = ctl(label, fld, lim, scale, grp)
        row = row + 1;
        lbl = uilabel(cp, 'Text', label, 'FontSize', 11.5);
        place(lbl, row, 1);
        sl = uislider(cp, 'Limits', lim, 'MajorTicks', [], 'MinorTicks', [], 'Value', min(max(P.(fld)/scale, lim(1)), lim(2)));
        place(sl, row, 2);
        ed = uieditfield(cp, 'numeric', 'Value', P.(fld)/scale, 'ValueDisplayFormat', '%.4g', ...
                         'HorizontalAlignment', 'right');
        place(ed, row, 3);
        s = struct('sl', sl, 'ed', ed, 'fld', fld, 'scale', scale, 'grp', grp);
        if strcmp(grp, 'display')
            sl.ValueChangingFcn = @(~,e) valueIn(fld, e.Value, false);
        end
        sl.ValueChangedFcn = @(~,e) valueIn(fld, e.Value, true);
        ed.ValueChangedFcn = @(~,e) valueIn(fld, e.Value, true);
    end

    function valueIn(fld, v, final)
        c = C.(fld);
        P.(fld) = v*c.scale;
        c.ed.Value = v;
        if v >= c.sl.Limits(1) && v <= c.sl.Limits(2), c.sl.Value = v; end
        switch c.grp
            case 'display', if final, refreshHistory(); end, drawFrame();
            case 'instant', if final, evaluate(); end
            case 'solve',   if final, resolutionChanged(); end
        end
        loadAdvanced();
    end

    function v = getP(), v = P; end      % scripted access (nested: current values)
    function v = getR(), v = R; end
    function v = getO(), v = O; end

    function setParam(fld, v)            % scripted access: SI value
        P.(fld) = v;  syncControls();  loadAdvanced();
        if any(strcmp(fld, solveFields())), solveNow(false); else, evaluate(); end
    end

    function syncControls()
        fn = fieldnames(C);
        for i = 1:numel(fn)
            c = C.(fn{i});  v = P.(c.fld)/c.scale;
            c.ed.Value = v;
            c.sl.Value = min(max(v, c.sl.Limits(1)), c.sl.Limits(2));
        end
        ddDim.Value = P.dim;
        [~, j] = min(abs(dhList - P.dh));  ddDh.Value = dhList(j);  P.dh = dhList(j);
        C.t_on.sl.Limits = [0 P.t_end];
    end

    function f = solveFields()
        [~, f] = v5_key(P);
    end

    function loadPreset(name)
        keep = struct('dim', P.dim, 'dh', P.dh);
        P = v5_defaults(name);  P.dim = keep.dim;  P.dh = keep.dh;
        ddPreset.Value = name;
        syncControls();  loadAdvanced();
        solveNow(false);
    end

% =========================================================================
%  SOLVE / EVALUATE
    function resolutionChanged()
        P.dim = ddDim.Value;  P.dh = ddDh.Value;
        est = estimate(P);
        if isKey(cache, v5_key(P)) || est < 4
            solveNow(false);
        else
            lblRun.Text = sprintf('Settings changed: press Solve.\nEstimated %s for %s nodes.\n%s', ...
                                  fmtTime(est), commas(nodes(P)), runText());
            lblRun.FontColor = [0.7 0.2 0];
        end
    end

    function solveNow(~)
        if S.playing, togglePlay(); end
        P.dim = ddDim.Value;  P.dh = ddDh.Value;
        if P.d <= P.b
            uialert(fig, sprintf('The antenna axis must be further from the gel than the element radius (%.2f mm).', P.b*1e3), ...
                    'Geometry');  return
        end
        key = v5_key(P);
        if isKey(cache, key)
            R = cache(key);
        else
            est = estimate(P);
            dlg = [];
            if est > 0.8
                dlg = uiprogressdlg(fig, 'Title', 'Solving', 'Message', sprintf('about %s ...', fmtTime(est)));
            end
            try
                R = v5_solve(P, @(f, m) progress(dlg, f, m, est));
            catch err
                if ~isempty(dlg), close(dlg); end
                uialert(fig, err.message, 'Solve failed');  return
            end
            if ~isempty(dlg), close(dlg); end
            fitRates(R);
            cache(key) = R;  cacheOrder{end+1} = key;
            if numel(cacheOrder) > 6, remove(cache, cacheOrder{1}); cacheOrder(1) = []; end
        end
        lblRun.FontColor = [0 0 0];
        slT.Limits = [0 R.t(end)];
        C.t_on.sl.Limits = [0 R.t(end)];
        drawStatic();
        evaluate();
        refreshDerived();
    end

    function progress(dlg, f, m, est)
        if isempty(dlg), return, end
        dlg.Value = min(max(f, 0), 1);
        dlg.Message = sprintf('%s\n(estimated %s)', m, fmtTime(est));
        drawnow limitrate;
    end

    function evaluate()
        if isempty(R), return, end
        t0 = tic;
        O = v5_eval(R, P);
        O.t_eval = toc(t0);
        % colour limits and history, fixed for the whole run
        T = reshape(O.T, [], numel(R.t));
        O.peak = max(T, [], 1);
        O.surf = max(reshape(O.T(1,:,:), [], numel(R.t)), [], 1);
        O.Tlo = min([P.T_b P.T_0 P.T_room]);
        O.Thi = max(max(O.peak), O.Tlo + 1);
        S.nEval = S.nEval + 1;  O.stamp = S.nEval;   % new evaluation: colour limits are refreshed
        xlim(axProf, [floor(O.Tlo) - 1, ceil(O.Thi) + 1]);
        refreshHistory();
        refreshEnvelope();
        drawFrame();
    end

% =========================================================================
%  DRAWING
    function drawStatic()
        x = R.x*1e3;  z = R.z*1e3;
        hImg.XData = x([1 end]);  hImg.YData = z([1 end]);
        hIso.XData = x;  hIso.YData = z;  hIso.ZData = zeros(numel(z), numel(x));
        top = max(P.d + P.r_cable, P.d + P.b)*1e3 + 1.5;
        xa = [min(x(1), -P.l_f*1e3 - 2) max(x(end), P.L*1e3 + 1)];
        set(hWater, 'XData', xa([1 2 2 1]), 'YData', [0 0 top top]);
        dz = P.d*1e3;
        rc = P.r_cable*1e3;  rb = P.b*1e3;  rcond = max(P.a*1e3, 0.08);  lf = P.l_f*1e3;  Lm = P.L*1e3;
        set(hCable, 'XData', [xa(1) -lf -lf xa(1)], 'YData', dz + [-rc -rc rc rc]);
        set(hFeed,  'XData', [-lf 0 0 -lf], 'YData', dz + 1.05*[-rb -rb rb rb]);
        set(hIns,   'XData', [0 Lm Lm 0], 'YData', dz + [-rb -rb rb rb]);
        set(hCond,  'XData', [0 Lm Lm 0], 'YData', dz + [-rcond -rcond rcond rcond]);
        xlim(axMap, xa);  ylim(axMap, [z(end) top]);
        hOff.Value = P.t_on;
        xlim(axHist, [0 R.t(end)]);
        ylim(axProf, [z(end) 0]);
    end

    function refreshHistory()
        t = R.t;
        set(hCath,  'XData', t, 'YData', O.T_cath);
        set(hPeakL, 'XData', t, 'YData', O.peak);
        set(hSurf,  'XData', t, 'YData', O.surf);
        if all(isfinite(S.probe))
            [ix, iz] = probeIdx();
            set(hPrH, 'XData', t, 'YData', squeeze(O.T(iz, ix, :)), 'DisplayName', ...
                sprintf('probe (%.1f, %.1f) mm', R.x(ix)*1e3, R.z(iz)*1e3));
        else
            set(hPrH, 'XData', NaN, 'YData', NaN, 'DisplayName', 'probe (click the map)');
        end
        yl = [floor(min([O.Tlo, O.T_cath])) - 1, ceil(max([O.Thi, O.T_cath, P.T_tlc + P.bw_tlc])) + 2];
        ylim(axHist, yl);  xlim(axHist, [0 R.t(end)]);
        hOff.Value = P.t_on;
        set(hBand, 'XData', [0 R.t(end) R.t(end) 0], 'YData', P.T_tlc + [0 0 P.bw_tlc P.bw_tlc]);
    end

    function drawFrame()
        if isempty(O), return, end
        t0 = tic;
        S.n = min(S.n, numel(R.t));
        n = S.n;  x = R.x*1e3;  z = R.z*1e3;
        T = double(O.T(:,:,n));
        v = ddView.Value;
        key = sprintf('%s|%g|%g|%g', v, P.T_tlc, P.bw_tlc, O.stamp);
        if ~strcmp(key, S.view)          % colour bar and limits: only when they change
            S.view = key;  setColours(v);
        end
        switch v
            case {'Pigment (as photographed)', 'Pigment (vivid)'}
                hImg.CData = v5_tlc(T, P.T_tlc, P.bw_tlc, tern(strcmp(v, 'Pigment (vivid)'), 'vivid', 'photo'));
            otherwise
                switch v
                    case 'Temperature (°C)',        F = T;
                    case 'Rise above water (K)',    F = T - P.T_b;
                    case 'Dielectric part θ_d (K)', F = double(O.theta_d(:,:,n));
                    case 'Diffusion part θ_s (K)',  F = double(O.theta_s(:,:,n));
                    case 'SAR (W/kg)',              F = R.Q*O.P_del/P.rho;
                end
                lim = S.lim;
                idx = round(1 + 255*min(max((F - lim(1))/(lim(2) - lim(1)), 0), 1));
                hImg.CData = reshape(S.cmap(idx, :), [size(F) 3]);
        end
        hIso.ZData = T;  hIso.LevelList = P.T_tlc;  hIso.Visible = matlab.lang.OnOffSwitchState(cbIso.Value);
        [Tpk, i] = max(T(:));  [izp, ixp] = ind2sub(size(T), i);
        set(hPk, 'XData', x(ixp), 'YData', z(izp));
        if all(isfinite(S.probe)), [ix, iz] = probeIdx(); set(hProbe, 'XData', x(ix), 'YData', z(iz));
        else, set(hProbe, 'XData', NaN, 'YData', NaN); end
        ttlMap.Text = sprintf('t = %.0f s   %s   P_{in} = %.0f W → %.1f W at the antenna   peak %.1f °C at (%.1f, %.1f) mm', ...
                                R.t(n), tern(O.env.p_now(n) > 0, 'power ON', 'power OFF'), P.P_set, O.P_del, Tpk, x(ixp), z(izp));
        % time readouts
        lblT.Text = sprintf('t = %.0f s', R.t(n));
        if ~S.playing || mod(n, 3) == 0, slT.Value = R.t(n); end
        hNow.Value = R.t(n);
        % profile through the hottest column (or the probe column)
        if all(isfinite(S.probe)), [ixc, ~] = probeIdx(); else, ixc = ixp; end
        hPrT.XData = T(:, ixc);  hPrT.YData = z;
        hPrD.XData = P.T_b + double(O.theta_d(:, ixc, n));  hPrD.YData = z;
        hPrS.XData = P.T_b + double(O.theta_s(:, ixc, n));  hPrS.YData = z;
        hPrX.Value = P.T_tlc;
        ttlProf.Text = sprintf('depth profile at x = %.1f mm (%s)', x(ixc), tern(all(isfinite(S.probe)), 'probe', 'hottest column'));
        hEnvN.Value = R.t(n);
        if cbPhoto.Value, drawPhoto(R.t(n)); end
        if ~S.playing || mod(n, 4) == 0 || n == numel(R.t)
            readout(n, T, Tpk, ixp, izp);
            lblRun.Text = runText();
        end
        O.t_draw = toc(t0);
        drawnow limitrate;
    end

    function setColours(v)
        switch v
            case {'Pigment (as photographed)', 'Pigment (vivid)'}
                [~, cmap, Tt] = v5_tlc([], P.T_tlc, P.bw_tlc, tern(strcmp(v, 'Pigment (vivid)'), 'vivid', 'photo'));
                lim = Tt([1 end]);  lab = 'pigment colour vs T (°C)';
            otherwise
                cmap = turbo(256);
                switch v
                    case 'Temperature (°C)',        lim = [O.Tlo O.Thi];  lab = 'T (°C)';
                    case 'Rise above water (K)',    lim = [min(0, O.Tlo - P.T_b) O.Thi - P.T_b];  lab = 'T − T_b (K)';
                    case 'Dielectric part θ_d (K)', lim = [0 double(max(O.theta_d(:)))];  lab = 'θ_d (K)';
                    case 'Diffusion part θ_s (K)',  lim = [min(0, double(min(O.theta_s(:)))) max(double(max(O.theta_s(:))), 1e-3)];  lab = 'θ_s (K)';
                    case 'SAR (W/kg)'
                        Fs = sort(R.Q(:)*O.P_del/P.rho);  lim = [0 Fs(ceil(0.997*numel(Fs)))];
                        lab = 'SAR = σE_{rms}^2/ρ (W/kg), while on';
                end
        end
        if lim(2) <= lim(1), lim(2) = lim(1) + 1; end
        S.cmap = cmap;  S.lim = lim;
        colormap(axMap, cmap);  clim(axMap, lim);  cbMap.Label.String = lab;
    end

    function refreshEnvelope()
        E = O.env;
        set(hEnvF, 'XData', E.ts, 'YData', O.P_del*E.ps);
        xa = [R.t(1:end-1); R.t(2:end)];  ya = O.P_del*[E.pbar; E.pbar];
        set(hEnvA, 'XData', xa(:).', 'YData', ya(:).');
        hEnvM.Value = O.P_del*E.mean_on;
        hEnvM.Label = sprintf('mean %.1f W', O.P_del*E.mean_on);
        lim = [0 R.t(end) 0 max(1, 1.45*O.P_del)];                          % headroom for the legend
        if isempty(S.envLim) || isequal([axEnv.XLim axEnv.YLim], S.envLim)  % not zoomed/panned by hand
            axis(axEnv, lim);  S.envLim = lim;
            resetplotview(axEnv, 'SaveCurrentView');                          % "restore view" returns here
        end
        s = sprintf('Power envelope: %.1f W when on, duty %.0f %%', O.P_del, 100*P.duty);
        if P.t_rest > 0, s = [s sprintf(', %g s bursts / %g s rests', P.t_burst, P.t_rest)]; end
        s = [s sprintf(';  on %.0f %% of the %g s treatment', 100*E.mean_on, P.t_on)];
        if E.fine, s = [s sprintf('  (pulses shorter than 2 frames are averaged per %g s frame)', R.P.dt_frame)]; end
        S.envTitle = s;
        if ~cbPhoto.Value, ttlPh.Text = s; end
    end

    function togglePhoto()
        on = cbPhoto.Value;
        axPh.Visible = matlab.lang.OnOffSwitchState(on);   hPh.Visible = axPh.Visible;
        axEnv.Visible = matlab.lang.OnOffSwitchState(~on);
        set([hEnvF hEnvA hEnvM hEnvN], 'Visible', axEnv.Visible);
        axEnv.Legend.Visible = axEnv.Visible;
        gr.RowHeight(2:3) = tern(on, {0, '1x'}, {'1x', 0});
        photos.shown = 0;
        if on, drawPhoto(R.t(S.n)); else, ttlPh.Text = S.envTitle; end
    end

    function drawPhoto(t)
        if isempty(photos.t), ttlPh.Text = 'photos/ not found'; return, end
        [dt, k] = min(abs(photos.t - t));
        if k ~= photos.shown
            hPh.CData = photos.img{k};  photos.shown = k;
        end
        s = sprintf('Photographed 50 W run, t = %d s', photos.t(k));
        if dt > 20, s = [s sprintf(' (nearest to %.0f s)', t)]; end
        if ~strcmp(ddPreset.Value, 'photo') || abs(P.L - 12e-3) > 1e-4
            s = [s ' (12 mm antenna, 25 °C bath: compare colours, not sizes)'];
        end
        ttlPh.Text = s;
    end

    function readout(n, T, Tpk, ixp, izp)
        x = R.x*1e3;  z = R.z*1e3;
        above = T >= P.T_tlc;
        cellA = (x(2) - x(1))*(z(1) - z(2));
        first = find(above(:, ixp), 1);  % uncoloured layer above the hottest column
        if ~any(above(:)),     topTxt = 'nothing coloured yet';
        elseif isempty(first), topTxt = 'hottest column not coloured';
        else,                  topTxt = sprintf('clear top layer %.2f mm', -z(first));
        end
        if R.dim == 3
            gelTxt = sprintf('%.2f W absorbed in the gel block (%.0f %% of it)', R.P_gel*O.P_del, 100*R.P_gel);
        else
            gelTxt = '2-D slice: no absorbed-power total';
        end
        txtOut.Value = sprintf([ ...
            'POWER\n' ...
            '  at the antenna  %.1f W = (%.0f − %.0f W) × (1 − %.2f) when on\n' ...
            '  envelope        %.0f %% on, mean %.1f W, %d pulse(s)\n' ...
            '  %s\n' ...
            'CATHETER (lumped)\n' ...
            '  P_c = γ P_fwd P_drop = %.2f W  ->  ΔT_c = %.1f K, τ_c = %.0f s\n' ...
            '  tip now %.1f °C;  α = min(1, β/(d√f)) = %.3f\n' ...
            '  water over the catheter  +%.2f K above T_b\n' ...
            'GEL PLANE at t = %.0f s\n' ...
            '  peak %.1f °C at x = %.1f, z = %.1f mm\n' ...
            '     dielectric θ_d %.1f K + diffusion θ_s %.2f K\n' ...
            '  surface max %.1f °C;  %s\n' ...
            '  coloured (T ≥ %.1f °C) area %.0f mm²\n' ...
            'FIELD\n' ...
            '  k_L = %.0f %+.0fj 1/m,  λ_L/4 = %.1f mm,  σ = %.2f S/m\n' ...
            '  Z_in = %.1f %+.1fj Ω,  peak SAR %.0f W/kg (on)'], ...
            O.P_del, P.P_set, P.P_refl, P.P_dropoff, 100*O.env.mean_on, O.P_del*O.env.mean_on, O.env.n_pulses, gelTxt, ...
            O.P_c, O.dTc, O.tau_c, O.T_cath(n), O.alpha, O.Theta(n), ...
            R.t(n), Tpk, x(ixp), z(izp), O.theta_d(izp, ixp, n), O.theta_s(izp, ixp, n), ...
            max(T(1, :)), topTxt, P.T_tlc, nnz(above)*cellA, ...
            real(R.info.k_L), imag(R.info.k_L), R.info.lambda_L/4*1e3, R.info.sigma, ...
            real(R.info.Z_in), imag(R.info.Z_in), max(R.Q(:))*O.P_del/P.rho);
    end

    function s = runText()
        if isempty(R), s = ''; return, end
        tm = R.timing;
        s = sprintf(['last solve  %s  (%s, Δh %.2f mm)\n' ...
                     '  E-field %s · modes %s · frames %s\n' ...
                     '  %s nodes, %d frames of %g s\n' ...
                     'evaluate %.0f ms · draw %.0f ms'], ...
                     fmtTime(tm.total), tern(R.dim == 3, '3-D', '2-D'), R.P.dh*1e3, ...
                     fmtTime(tm.E), fmtTime(tm.modes), fmtTime(tm.frames), commas(R.N), numel(R.t), R.P.dt_frame, ...
                     1e3*getf(O, 't_eval'), 1e3*getf(O, 't_draw'));
    end

% =========================================================================
%  TIME, PLAY, PROBE, VIEW
    function setView(v), ddView.Value = v; drawFrame(); end

    function setTime(t), setFrame(nearestFrame(t)); end

    function n = nearestFrame(t)
        if isempty(R), n = 1; return, end
        [~, n] = min(abs(R.t - t));
    end

    function setFrame(n)
        S.n = max(1, min(n, numel(R.t)));  drawFrame();
    end

    function togglePlay()
        if S.playing
            S.playing = false;  stop(tmr);  bPlay.Text = '▶ Play';
        else
            if S.n >= numel(R.t), S.n = 1; end
            S.playing = true;  bPlay.Text = '❚❚ Pause';  start(tmr);
        end
    end

    function tick()
        step = max(1, round(ddSpeed.Value*tmr.Period/R.P.dt_frame));
        S.n = S.n + step;
        if S.n >= numel(R.t)
            S.n = numel(R.t);  S.playing = false;  stop(tmr);  bPlay.Text = '▶ Play';
        end
        drawFrame();
    end

    function mapClick(e)
        p = e.IntersectionPoint;
        S.probe = [p(1) p(2)]*1e-3;
        refreshHistory();  drawFrame();
    end

    function [ix, iz] = probeIdx()
        [~, ix] = min(abs(R.x - S.probe(1)));  [~, iz] = min(abs(R.z - S.probe(2)));
    end

    function exportFrame()
        d = fullfile(here, 'exports');  if ~exist(d, 'dir'), mkdir(d); end
        stem = fullfile(d, sprintf('v5_%s_t%03.0fs_%s', P.preset, R.t(S.n), char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'))));
        exportgraphics(axMap, [stem '.png'], 'Resolution', 200);
        Pexp = P;  t = R.t(S.n);  Tplane = O.T(:,:,S.n);  x = R.x;  z = R.z;
        save([stem '.mat'], 'Pexp', 't', 'Tplane', 'x', 'z');
        uialert(fig, sprintf('Saved %s.png and .mat', stem), 'Exported', 'Icon', 'success');
    end

% =========================================================================
%  ADVANCED TAB
    function loadAdvanced()
        fn = fieldnames(EF);
        for i = 1:numel(fn)
            v = P.(fn{i})/EF.(fn{i}).scale;
            EF.(fn{i}).h.Value = num2str(v, 8);
        end
        refreshDerived();
    end

    function Pn = readAdvanced()
        Pn = P;  fn = fieldnames(EF);
        for i = 1:numel(fn)
            s = strtrim(EF.(fn{i}).h.Value);
            v = str2double(s);
            if strcmpi(s, 'nan'), v = NaN; elseif isnan(v)
                error('"%s": cannot read "%s"', fn{i}, s);
            end
            v = v*EF.(fn{i}).scale;
            if abs(v - P.(fn{i})) <= 1e-9*abs(P.(fn{i})), v = P.(fn{i}); end   % display rounding only
            Pn.(fn{i}) = v;
        end
        if Pn.a >= Pn.b, error('conductor radius a must be smaller than the element radius b'); end
        if Pn.d <= Pn.b, error('antenna axis to gel (%.2f mm) must exceed b', Pn.d*1e3); end
        if Pn.h_s < 0 || Pn.h_w < 0, error('h_s, h_w must be >= 0'); end
        if Pn.dt_frame <= 0 || Pn.t_end < 4*Pn.dt_frame, error('need dt_frame > 0 and t_end >= 4 dt_frame'); end
        if any([Pn.gel_L Pn.gel_H Pn.gel_W Pn.rho Pn.c Pn.k Pn.k_cat Pn.m_cat Pn.c_cat] <= 0), error('sizes and properties must be > 0'); end
    end

    function refreshDerived()
        try
            Pn = readAdvanced();
        catch err
            lblAdv.Text = err.message;  return
        end
        lblAdv.Text = '';
        [~, info] = v5_erms(Pn, 0, Pn.d);
        P_fwd = max(Pn.P_set - Pn.P_refl, 0);  P_c = Pn.gamma*P_fwd*Pn.P_dropoff;
        same = ~isempty(R) && strcmp(v5_key(Pn), R.key);
        txtDer.Value = sprintf([ ...
            'DERIVED   (%s)\n\n' ...
            'power      at the antenna %.2f W,  heating the tip P_c = %.2f W\n' ...
            'catheter   ΔT_c = P_c/k_c = %.1f K,  τ_c = m c/k_c = %.1f s\n' ...
            'surface    α = min(1, β/(d√f)) = %.3f  (d = %.2f mm, f = %.2g L/min)\n' ...
            'gel        σ = ωε0ε'''' = %.2f S/m,  a = k/(ρc) = %.3g m²/s\n' ...
            '           √(a t) = %.1f mm at 92 s, %.1f mm at 360 s\n' ...
            'antenna    k_L = %.1f %+.1fj 1/m,  λ_L/4 = %.1f mm,  L = %.2f λ_L\n' ...
            '           Z_in = %.1f %+.1fj Ω,  I0 = %.3f A per √W\n\n' ...
            'grid       %s nodes (%s),  %d frames\n' ...
            'estimated solve time  %s\n'], ...
            tern(same, 'matches the current solution', 'CHANGED: Apply re-solves'), ...
            P_fwd*(1 - Pn.P_dropoff), P_c, P_c/Pn.k_cat, Pn.m_cat*Pn.c_cat/Pn.k_cat, ...
            min(1, Pn.beta/(Pn.d*sqrt(Pn.flow))), Pn.d*1e3, Pn.flow, ...
            info.sigma, Pn.k/(Pn.rho*Pn.c), 1e3*sqrt(Pn.k/(Pn.rho*Pn.c)*92), 1e3*sqrt(Pn.k/(Pn.rho*Pn.c)*360), ...
            real(info.k_L), imag(info.k_L), info.lambda_L/4*1e3, Pn.L/info.lambda_L, ...
            real(info.Z_in), imag(info.Z_in), info.I0, ...
            commas(nodes(Pn)), tern(Pn.dim == 3, '3-D', '2-D'), floor(Pn.t_end/Pn.dt_frame) + 1, fmtTime(estimate(Pn)));
    end

    function applyAdvanced()
        try
            Pn = readAdvanced();
        catch err
            uialert(fig, err.message, 'Cannot read the constants');  return
        end
        resolve = ~strcmp(v5_key(Pn), v5_key(P));
        P = Pn;  syncControls();
        if resolve, solveNow(true); else, evaluate(); end
        refreshDerived();
    end

% =========================================================================
%  RUNTIME ESTIMATE
    function N = nodes(Pq)
        nx = round(Pq.gel_L/Pq.dh) + 1;  nz = round(Pq.gel_H/Pq.dh) + 1;
        ny = tern(Pq.dim == 3, round(Pq.gel_W/2/Pq.dh) + 1, 1);
        N = nx*ny*nz;
    end

    function t = estimate(Pq)
        nx = round(Pq.gel_L/Pq.dh) + 1;
        nrho = tern(Pq.dim == 3, 300, round(Pq.gel_H/Pq.dh) + 1);
        dxs = min(0.05e-3, Pq.L/200);
        nseg = round(Pq.L/dxs) + round(Pq.l_f/dxs);
        nt_ = floor(Pq.t_end/Pq.dt_frame) + 1;
        t = rate.cE*nx*nrho*nseg + rate.cF*nodes(Pq)*nt_ + 0.05;
    end

    function fitRates(Rs)
        tm = Rs.timing;
        if tm.E > 0.05, rate.cE = tm.E/(Rs.nx*tm.nrho*tm.nseg); end
        if tm.frames > 0.05, rate.cF = (tm.frames + tm.modes)/(Rs.N*numel(Rs.t)); end
    end

% =========================================================================
    function ph = readPhotos()
        ph = struct('t', [], 'img', {{}}, 'shown', 0);
        f = fullfile(here, 'photos', 'photos.csv');
        if ~exist(f, 'file'), return, end
        tb = readtable(f, 'TextType', 'string');
        for i = 1:height(tb)
            ph.t(i) = tb.t_s(i);
            ph.img{i} = imread(fullfile(here, 'photos', tb.file(i)));
        end
    end

    function closeApp()
        try
            stop(tmr);  delete(tmr);
        catch
        end
        delete(fig);
    end
end

% =========================================================================
function s = tern(c, a, b)
if c, s = a; else, s = b; end
end

function v = getf(s, f)
if isstruct(s) && isfield(s, f), v = s.(f); else, v = NaN; end
end

function s = fmtTime(t)
if t < 1
    s = sprintf('%.0f ms', 1e3*t);
elseif t < 90
    s = sprintf('%.1f s', t);
else
    s = sprintf('%.1f min', t/60);
end
end

function s = commas(n)
s = regexprep(sprintf('%d', round(n)), '(\d)(?=(\d{3})+$)', '$1,');
end
