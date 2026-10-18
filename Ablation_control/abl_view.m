function fig = abl_view(H, S, P, info)
% ABL_VIEW  The slider window used by abl_plot and abl_pid: the layout of the V5
% GUI, with play / reset / time slider / speed, and a choice of colouring
% (temperature, or the V5 pigment colours that match the photographs).
%
%   fig = abl_view(H, S, P, info)
%
% Inputs:   H      history from abl_plot or abl_pid, one entry per frame:
%                  H.t, H.P_gen, H.P_ant, H.peak, H.surf, H.Tc, H.depth, and H.T (nx x nz x frames)
%           S      final simulation state (grid, final T_max)
%           P      parameters (abl.params)
%           info   struct with
%                    .envTitle     text above the power plot
%                    .name         window name (optional)
%                    .depthLimit   optional: depth of the spared layer [m], drawn as a dotted line
% Output:   fig    the window
%
% The window:
%   left, top      temperature on the antenna plane at the chosen time, with the
%                  ablated zone (everywhere that has reached $T_{abl}$ so far) outlined
%   left, middle   play / reset / time slider / playback speed / colouring
%   left, bottom   temperatures and ablation depth over time
%   right, top     the power reaching the antenna over time
%   right, bottom  temperature against depth through the hottest column

if ~isfield(info, 'name'), info.name = 'Ablation model'; end
if ~isfield(info, 'depthLimit'), info.depthLimit = []; end
t  = H.t;  nt = numel(t);
P_max  = max([H.P_ant, 1]);                                 % highest power at the antenna, for the plot range  [W]
mean_P = sum(H.P_ant(2:end).*diff(t))/max(t(end), eps);     % average power at the antenna over the run  [W]

% --- 1. fixed things for the whole window ---
x   = S.x*1e3;  z = S.z*1e3;                                % grid in mm
Tlo = min(P.T_b, P.T_0);                                    % colour scale: fixed over the whole run
Thi = max([H.peak, H.Tc, P.T_abl + 2, Tlo + 1]);
[~, ic] = max(S.Tmax(:));  [ixc, ~] = ind2sub(size(S.T), ic);   % hottest column over the whole run (for the profile)
ui.n = nt;  ui.playing = false;  ui.timer = [];             % window state: frame shown, playing or not
% Pigment views: the V5 thermochromic colouring, which matches the photographs,
% but with the colour starting at $T_{abl}$, so only gel above the cutoff lights up.
bw_pig = 23;                                                % [K] colour start -> blue; V5's value fitted to the photographs (v5_defaults P.bw_tlc)

% --- 2. the window: same layout as the V5 GUI ---
fig = uifigure('Name', info.name, 'Position', [60 60 1500 860], 'CloseRequestFcn', @(~,~) closeWindow());
gM  = uigridlayout(fig, [1 2], 'ColumnWidth', {'1.5x', '1x'});

% centre column: map title, map, time controls, history
gc = uigridlayout(gM, [4 1], 'RowHeight', {62, '1.7x', 34, '1x'}, 'Padding', [0 0 0 0], 'RowSpacing', 4);
ttlMap  = uilabel(gc, 'Text', '', 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Interpreter', 'tex', 'WordWrap', 'on');
axMap   = uiaxes(gc);
tr      = uigridlayout(gc, [1 6], 'ColumnWidth', {80, 60, '1x', 90, 110, 200}, 'Padding', [0 0 0 0]);
bPlay   = uibutton(tr, 'Text', '▶ Play', 'ButtonPushedFcn', @(~,~) togglePlay());
uibutton(tr, 'Text', '⟲ 0 s', 'ButtonPushedFcn', @(~,~) showFrame(1));
slT     = uislider(tr, 'Limits', [0 max(t(end), eps)], 'Value', t(end), ...
                   'ValueChangingFcn', @(~,e) setTime(e.Value), 'ValueChangedFcn', @(s,~) setTime(s.Value));
lblT    = uilabel(tr, 'Text', '', 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
ddSpeed = uidropdown(tr, 'Items', {'10 s / s', '30 s / s', '60 s / s', '120 s / s'}, 'ItemsData', {10, 30, 60, 120}, 'Value', 30);
ddView  = uidropdown(tr, 'Items', {'Temperature (°C)', 'Pigment (as photographed)', 'Pigment (vivid)'}, ...
                     'Value', 'Temperature (°C)', 'ValueChangedFcn', @(~,~) setView());
axHist  = uiaxes(gc);

% right column: power envelope and depth profile
gr = uigridlayout(gM, [4 1], 'RowHeight', {22, '1x', 22, '1.3x'}, 'Padding', [0 0 0 0], 'RowSpacing', 2);
uilabel(gr, 'Text', info.envTitle, ...
        'HorizontalAlignment', 'center', 'FontSize', 11);
axEnv  = uiaxes(gr);
uilabel(gr, 'Text', sprintf('depth profile at x = %.1f mm (hottest column)', x(ixc)), 'HorizontalAlignment', 'center', 'FontSize', 11);
axProf = uiaxes(gr);

% --- 2a. the temperature map ---
warning('off', 'MATLAB:contour:ConstantData');              % no ablated zone yet: nothing to outline
hold(axMap, 'on');  axMap.Toolbar.Visible = 'off';  disableDefaultInteractivity(axMap);
d_mm = P.d*1e3;  r_cable = 1.5;                             % antenna axis height and coax radius (drawing only)  [mm]
rb = P.b*1e3;  rcond = max(P.a*1e3, 0.08);  lf = P.l_f*1e3;  Lm = P.L*1e3;
top = d_mm + max(r_cable, rb) + 1.5;                        % top of the drawing  [mm]
xa  = [min(x(1), -lf - 2), max(x(end), Lm + 1)];            % x range of the drawing  [mm]
patch(axMap, xa([1 2 2 1]), [0 0 top top], [0.44 0.67 0.75], 'EdgeColor', 'none');          % water above the gel
if max(abs(diff(x, 2))) < 1e-6 && max(abs(diff(z, 2))) < 1e-6                              % evenly spaced grid
    hImg = imagesc(axMap, x, z, double(H.T(:,:,nt)).');                                     % gel temperature
else                                                                                        % graded grid: an image would stretch it, so draw a surface
    hImg = surface(axMap, x, z, zeros(numel(z), numel(x)), double(H.T(:,:,nt)).', 'EdgeColor', 'none', 'FaceColor', 'interp');
    axMap.SortMethod = 'childorder';                                                        % keep the lines drawn after it on top
end
[~, hIso] = contour(axMap, x, z, S.Tmax.', [P.T_abl P.T_abl], '--', 'Color', 'w', 'LineWidth', 1.4);   % ablated zone edge
patch(axMap, [xa(1) -lf -lf xa(1)], d_mm + [-r_cable -r_cable r_cable r_cable], [0.86 0.58 0.48], 'EdgeColor', [0.6 0.35 0.3]);  % coax cable
patch(axMap, [-lf 0 0 -lf], d_mm + 1.05*[-rb -rb rb rb], [0.88 0.88 0.82], 'EdgeColor', [0.5 0.5 0.5]);                          % feed block
patch(axMap, [0 Lm Lm 0], d_mm + [-rb -rb rb rb], [0.86 0.93 0.97], 'FaceAlpha', 0.7, 'EdgeColor', [0.55 0.65 0.7]);             % insulation
patch(axMap, [0 Lm Lm 0], d_mm + [-rcond -rcond rcond rcond], [0.97 0.97 0.99], 'EdgeColor', [0.6 0.6 0.65]);                    % conductor
if ~isempty(info.depthLimit)                                 % spared layer: no ablation allowed above this line
    plot(axMap, xa, -1e3*info.depthLimit*[1 1], ':', 'Color', [1 0.3 0.6], 'LineWidth', 1.6);
end
hPk = plot(axMap, NaN, NaN, 'o', 'Color', [0.1 0.1 0.1], 'MarkerSize', 9, 'LineWidth', 1.4);                                     % peak
set(axMap, 'YDir', 'normal', 'Color', [0.80 0.83 0.85], 'Layer', 'top');
axis(axMap, 'image');  box(axMap, 'on');  xlim(axMap, xa);  ylim(axMap, [z(end) top]);
cb = colorbar(axMap);                                       % colours and labels are set by setView
xlabel(axMap, 'x along the antenna (mm), feed at 0');  ylabel(axMap, 'z (mm), gel surface at 0');

% --- 2b. history ---
hold(axHist, 'on');  grid(axHist, 'on');  axHist.Toolbar.Visible = 'off';
yyaxis(axHist, 'left');
plot(axHist, t, H.Tc,   '-', 'Color', [0.85 0.45 0.1], 'LineWidth', 1.6, 'DisplayName', 'catheter tip');
plot(axHist, t, H.peak, '-', 'Color', [0.8 0.1 0.1],  'LineWidth', 1.8, 'DisplayName', 'peak in the gel plane');
plot(axHist, t, H.surf, '-', 'Color', [0.1 0.4 0.8],  'LineWidth', 1.4, 'DisplayName', 'max at the gel surface');
yline(axHist, P.T_abl, '--', sprintf('T_{abl} = %g °C', P.T_abl), 'LineWidth', 1.2, ...
      'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
hNow = xline(axHist, t(end), ':', 'LineWidth', 1.4, 'HandleVisibility', 'off');            % the moment shown
ylim(axHist, [floor(Tlo) - 1, ceil(Thi) + 2]);  ylabel(axHist, 'T (°C)');
yyaxis(axHist, 'right');
plot(axHist, t, 1e3*H.depth, ':', 'Color', [0.1 0.55 0.3], 'LineWidth', 1.8, 'DisplayName', 'ablation depth');
ylabel(axHist, 'ablation depth (mm)');  axHist.YColor = [0.1 0.55 0.3];
ylim(axHist, [0 max(1, 1.15*1e3*max(H.depth))]);
xlim(axHist, [0 max(t(end), eps)]);  xlabel(axHist, 'time since power on (s)');
legend(axHist, 'Location', 'northoutside', 'NumColumns', 4, 'FontSize', 9);

% --- 2c. power envelope ---
hold(axEnv, 'on');  grid(axEnv, 'on');
ts = reshape([t(1:end-1); t(2:end)], 1, []);                % stair-step: each step's power over its interval
ps = reshape([H.P_ant(2:end); H.P_ant(2:end)], 1, []);
area(axEnv, ts, ps, 'FaceColor', [0.85 0.35 0.2], 'FaceAlpha', 0.35, 'EdgeColor', [0.7 0.2 0.1], 'DisplayName', 'at the antenna');
yline(axEnv, mean_P, ':', sprintf('mean %.1f W', mean_P), 'Color', [0.1 0.4 0.8], 'LineWidth', 1.4, ...
      'LabelHorizontalAlignment', 'left', 'DisplayName', 'mean over treatment');
hEnvN = xline(axEnv, t(end), ':', 'LineWidth', 1.4, 'HandleVisibility', 'off');            % the moment shown
axis(axEnv, [0 max(t(end), eps) 0 1.45*P_max]);
xlabel(axEnv, 'time since power on (s)');  ylabel(axEnv, 'P at the antenna (W)');
legend(axEnv, 'Location', 'north', 'NumColumns', 2, 'FontSize', 8);

% --- 2d. depth profile ---
hold(axProf, 'on');  grid(axProf, 'on');  axProf.Toolbar.Visible = 'off';
hPrT = plot(axProf, NaN, NaN, '-', 'SeriesIndex', 'none', 'LineWidth', 1.8, 'DisplayName', 'T now');
hPrM = plot(axProf, NaN, NaN, '--', 'Color', [0.8 0.1 0.1], 'LineWidth', 1.3, 'DisplayName', 'T_{max} (hottest so far)');
if ~isempty(info.depthLimit)
    yline(axProf, -1e3*info.depthLimit, ':', 'ABLATION\_DEPTH', 'Color', [1 0.3 0.6], 'LineWidth', 1.6, ...
          'HandleVisibility', 'off', 'LabelHorizontalAlignment', 'left');
end
xline(axProf, P.T_abl, '--', 'T_{abl}', 'Color', [0.2 0.6 0.5], 'LineWidth', 1.2, 'HandleVisibility', 'off');
hDep = yline(axProf, 0, ':', '', 'Color', [0.1 0.55 0.3], 'LineWidth', 1.3, 'HandleVisibility', 'off', 'Visible', 'off');
xlim(axProf, [floor(Tlo) - 1, ceil(Thi) + 1]);  ylim(axProf, [z(end) 0]);
xlabel(axProf, 'T (°C)');  ylabel(axProf, 'z (mm)');
legend(axProf, 'Location', 'southeast', 'FontSize', 8);

setView();                                                  % colour the map and open on the end of the treatment

% =========================================================================
% WINDOW CALLBACKS (nested functions: they can see the variables above)
% =========================================================================

    function showFrame(n)
        % draw frame n: map, ablated zone, peak, title, "now" lines and the profile
        ui.n = n;
        Tn = double(H.T(:,:,n));                            % temperature at this moment
        Tm = double(max(H.T(:,:,1:n), [], 3));              % hottest so far, which is what abl.depth uses
        if strcmp(ddView.Value, 'Temperature (°C)')
            hImg.CData = Tn.';                              % temperatures, coloured through the colour map
        else
            hImg.CData = pigment(Tn.', P.T_abl, bw_pig, pigStyle());   % ready-made photo colours (RGB)
        end
        hIso.ZData = Tm.';
        [Tpk, ip] = max(Tn(:));  [ixp, izp] = ind2sub(size(Tn), ip);
        set(hPk, 'XData', x(ixp), 'YData', z(izp));
        k = max(n, 2);                                      % the step that ends at this frame
        if H.P_ant(k) > 0, state = 'power ON'; else, state = 'power OFF'; end
        ttlMap.Text = sprintf(['t = %.0f s   %s   P_{gen} = %.0f W \\rightarrow %.1f W at the antenna   ' ...
                               'peak %.1f °C at (%.1f, %.1f) mm\n%s dashed: ablated zone ' ...
                               '(T_{max} \\geq T_{abl} = %g °C), depth %.2f mm'], ...
                               t(n), state, H.P_gen(k), H.P_ant(k), Tpk, x(ixp), z(izp), ui.isoName, P.T_abl, 1e3*H.depth(n));
        hNow.Value = t(n);  hEnvN.Value = t(n);
        set(hPrT, 'XData', Tn(ixc,:), 'YData', z);
        set(hPrM, 'XData', Tm(ixc,:), 'YData', z);
        if H.depth(n) > 0
            set(hDep, 'Value', -1e3*H.depth(n), 'Label', sprintf('depth %.2f mm', 1e3*H.depth(n)), 'Visible', 'on');
        else
            hDep.Visible = 'off';
        end
        lblT.Text = sprintf('t = %.0f s', t(n));
        slT.Value = t(n);
        drawnow limitrate;
    end

    function setView()
        % switch the map between plain temperatures and the V5 pigment colours
        if strcmp(ddView.Value, 'Temperature (°C)')
            colormap(axMap, turbo(256));  clim(axMap, [Tlo Thi]);
            cb.Label.String = 'T (°C)';
            hIso.LineColor = 'w';  ui.isoName = 'white';            % white shows best on turbo
        else
            [~, cmap, Tt] = pigment([], P.T_abl, bw_pig, pigStyle());  % the colour bar for the pigment
            colormap(axMap, cmap);  clim(axMap, Tt([1 end]));
            cb.Label.String = sprintf('pigment colour vs T (°C), starts at T_{abl} = %g °C', P.T_abl);
            hIso.LineColor = [0.15 0.15 0.15];  ui.isoName = 'dark';  % the pigment colours are pale
        end
        showFrame(ui.n);
    end

    function s = pigStyle()
        % which pigment colours: the photographed ones or the saturated ones
        if strcmp(ddView.Value, 'Pigment (vivid)'), s = 'vivid'; else, s = 'photo'; end
    end

    function setTime(tv)
        % slider moved: show the frame nearest to time tv
        [~, n] = min(abs(t - tv));
        showFrame(n);
    end

    function togglePlay()
        % start or stop playback; it runs at the chosen speed (simulated seconds per real second)
        if ui.playing
            stopPlay();  return
        end
        if ui.n >= nt, showFrame(1); end                    % at the end: start again from 0 s
        ui.playing = true;  bPlay.Text = '❚❚ Pause';
        ui.t0 = t(ui.n);  ui.clock = tic;
        ui.timer = timer('ExecutionMode', 'fixedSpacing', 'Period', 0.05, 'BusyMode', 'drop', ...
                         'TimerFcn', @(~,~) tick());
        start(ui.timer);
    end

    function tick()
        % one playback frame: jump to the time reached at the chosen speed
        if ~isvalid(fig), stopPlay(); return, end
        tv = ui.t0 + ddSpeed.Value*toc(ui.clock);
        if tv >= t(end)
            stopPlay();  showFrame(nt);  return
        end
        setTime(tv);
    end

    function stopPlay()
        ui.playing = false;
        if ~isempty(ui.timer) && isvalid(ui.timer), stop(ui.timer); delete(ui.timer); end
        ui.timer = [];
        if isvalid(fig), bPlay.Text = '▶ Play'; end
    end

    function closeWindow()
        stopPlay();
        delete(fig);
    end
end


% =========================================================================
% PIGMENT: colour temperatures the way the thermochromic gel looks
% =========================================================================
function [rgb, cmap, Tticks] = pigment(T, T_start, bw, style)
% A copy of V5_convection/Coupled_heating_GUI/v5_tlc.m, so this folder works on
% its own. It colours a temperature field like the thermochromic gel in the
% photographs (50 W CW, 6 min).
%
%   rgb = pigment(T, T_start, bw, style)          T any 2-D array [degC] -> m x n x 3 colours
%   [~, cmap, Tticks] = pigment([], T_start, bw, style)   256-colour bar and its temperatures
%
% Inputs:   T_start   temperature where the colour starts  [degC]  (here $T_{abl}$)
%           bw        band width: colour start -> blue  [K]
%           style     'photo' = colours sampled from the photographs, 'vivid' = same hues, saturated
%
% Below T_start the gel keeps its own hazy blue-grey. Through the band the
% colour runs pink -> yellow -> green -> teal (the photographed core), then
% blue at the top of the band, and fades back to clear above it. The 'photo'
% colours are 9 x 9 pixel averages taken from photos/photo_332s.jpg.
s_k = [-0.06  0.00  0.12  0.35  0.55  0.78  1.00  1.12  1.35];   % position in the band: (T - T_start)/bw
switch style
    case 'photo'
        gel = [0.645 0.730 0.785];
        c_k = [gel;
               0.70 0.70 0.76;             % onset: first pink specks
               0.82 0.72 0.78;             % pink ring
               0.83 0.80 0.76;             % yellow
               0.70 0.86 0.78;             % green
               0.56 0.855 0.82;            % teal core
               0.50 0.72 0.86;             % blue, top of the band (not photographed)
               0.56 0.70 0.82;             % fading
               gel];                       % clear above the band
    otherwise
        gel = [0.62 0.71 0.78];
        c_k = [gel;
               0.74 0.66 0.74;
               0.88 0.62 0.70;
               0.90 0.82 0.62;
               0.66 0.86 0.64;
               0.46 0.82 0.74;
               0.30 0.58 0.86;
               0.40 0.60 0.80;
               gel];
end
colour = @(s) interp1(s_k, c_k, min(max(s(:), s_k(1)), s_k(end)));   % band position -> RGB
if isempty(T)
    Tticks = linspace(T_start - 0.15*bw, T_start + 1.45*bw, 256);
    cmap   = colour((Tticks - T_start)/bw);
    rgb    = [];
    return
end
rgb = reshape(colour((T(:) - T_start)/bw), [size(T) 3]);
end
