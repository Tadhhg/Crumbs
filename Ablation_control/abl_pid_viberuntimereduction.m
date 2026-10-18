function [H, S] = abl_pid_viberuntimereduction(ABLATION_DEPTH, T_END, FLOW)
% ABL_PID_VIBERUNTIMEREDUCTION  Same as abl_pid (closed-loop power control with
% two PIDs while the flow changes), but on a GRADED grid: the points are close
% together only near the antenna and the gel surface, where the gel heats up and
% the ablation boundary sits, and spread out in the far parts of the block.
% That gives fine resolution where it matters for less compute than an evenly
% fine grid. The grid is set by the "graded grid" settings below.
%
% Why it works: the gel's heat modes are products of 1-D shapes along x, y and
% z, and each 1-D line can be unevenly spaced (abl_forvibes.modes, abl_forvibes.graded). The step
% is still exact in time, so the small cells need no smaller time steps.
% Limits: the fine zone is a slab (fine in z means fine at every x and y), and
% the spacing must change gradually (GROWTH 1.1 to 1.2).
%
%   abl_pid_viberuntimereduction                  default run: spare 1.5 mm, 360 s,
%                                                 flow 5 -> 2 L/min at 120 s -> 10 L/min at 240 s
%   abl_pid_viberuntimereduction(ABLATION_DEPTH)  choose the spared surface layer, in mm (1.5) or m (1.5e-3)
%   abl_pid_viberuntimereduction(ABLATION_DEPTH, T_END, FLOW)
%                              T_END = run length [s]; FLOW = flow schedule [t1 f1; t2 f2; ...]
%                              ([s] and [L/min]), the flow jumps to f_i at time t_i.
%                              Leave any input [] to keep its default
%   [H, S] = abl_pid_viberuntimereduction(...)    also return the history H and the final state S
%
% Example:  abl_pid_viberuntimereduction(1, 300)    spare the top 1 mm, 5 minutes, default flow schedule
%
% With LIVE = true (below) a window shows the run as it happens, paced to the
% clock (SPEED = 1 is real time). Drag its flow slider to change the flow
% yourself mid-run and watch the PIDs respond; Stop ends the run early.
%
% Every DT seconds the controller looks at the temperature map on the antenna
% plane (a "perfect map", as if measured, e.g. by MR thermometry) and sets the
% generator power:
%
%   PID 1  safety, has priority. Keeps every point shallower than ABLATION_DEPTH
%          (the spared surface layer) below the ablation temperature. It regulates
%          the hottest temperature the layer is heading for to $T_{abl}$ - MARGIN, and gives
%          the highest power that keeps it there. "Heading for" = the hottest the layer
%          gets if the last power runs one more step and then stops for T_LOOK s, so
%          heat already on its way is counted before it arrives. A look-ahead check also simulates "this power for one
%          step, then switched off for T_LOOK seconds" (heat keeps spreading upward
%          after the power stops); if that would bring the layer within MARGIN/2 of
%          $T_{abl}$, the power is halved until it is safe (down to 0 W if needed).
%   PID 2  performance. Keeps the ablated area (everywhere that has reached $T_{abl}$)
%          growing at R_TARGET mm^2/s, between P_MIN and MAX_WATTS.
%
%   Applied power = min(PID 1, PID 2): PID 2 is followed unless PID 1 asks for less.
%   The loop that is not in charge keeps its integral term at the applied power
%   ("external reset"), so it can take over smoothly and cannot wind up.
%
% All model parameters come from abl_forvibes.params (abl_forvibes.m); the controller settings are below.

if nargin < 1 || isempty(ABLATION_DEPTH), ABLATION_DEPTH = 1.5e-3; end
if nargin < 2 || isempty(T_END), T_END = 360; end
if nargin < 3 || isempty(FLOW),  FLOW = [0 5; 120 2; 240 10]; end
if ABLATION_DEPTH >= 0.05, ABLATION_DEPTH = 1e-3*ABLATION_DEPTH; end   % given in mm (no gel layer is 5 cm thick)

% --- settings ---
% ABLATION_DEPTH (input 1, default 1.5 mm): spared surface layer, no ablation shallower than this  [m]
MAX_WATTS = 60;              % %*% highest generator power PID 2 may ask for  [W]
P_MIN     = 10;              % %*% lowest power PID 2 may ask for, so the antenna stays on (PID 1 can still go lower)  [W]
P_START   = 40;              % %*% power PID 2 starts at, instead of ramping up from P_MIN; PID 1 still cuts it if the layer gets too hot  [W]
DT        = 0.5;             % controller update interval  [s]
MODE      = 'pid';           % 'pid' = PID 2 grows the zone; 'hold' = PID 2 replaced by a fixed P_HOLD (PID 1 still guards)
P_HOLD    = 30;              %#ok<NASGU> generator power in 'hold' mode  [W]
MARGIN    = 2;               % PID 1 keeps the spared layer this far below $T_{abl}$  [K]
T_LOOK    = 10;              % look-ahead with the power off  [s]
R_TARGET  = 0.5;             % target growth rate of the ablated area  [mm^2/s]
TAU_R     = 5;               % smoothing time of the measured growth rate  [s]
% graded grid (fine near the antenna and the surface, coarse far away)
H_FINE    = 0.25e-3;         % spacing in the fine zone, across (y) and in depth (z)  [m]
FINE_Z    = 6e-3;            % fine zone reaches this far below the gel surface  [m]
FINE_Y    = 5e-3;            % fine zone reaches this far either side of the antenna plane  [m]
H_X       = 0.5e-3;          % spacing along the antenna (x); the heating changes slowly along it  [m]
H_COARSE  = 2e-3;            % largest spacing, far from the antenna  [m]
GROWTH    = 1.15;            % each gap may be this much bigger than the one before it  [-]
LIVE      = true;            % show the run live in a window, paced to the clock
SPEED     = 1;               % live playback speed: 1 = real time, 5 = five times faster (also adjustable in the window)
% PID gains (u = I + Kp e - Kd dy/dt, I += Ki e DT)
Kp1 = 8;    Ki1 = 0.4;   Kd1 = 0;      % PID 1: e = headroom in K  ->  W
Kp2 = 10;   Ki2 = 1;     Kd2 = 0;      % PID 2: e = growth-rate shortfall in mm^2/s  ->  W

% --- 1. set the model up ---
P = abl_forvibes.params;
xe = P.x0 + [0, P.gel(1)];                                  % ends of the block along x  [m]
P.grid.x = abl_forvibes.graded(xe(1), xe(2), xe(1), xe(2), H_X, H_COARSE, GROWTH);
P.grid.y = abl_forvibes.graded(0, P.gel(2), 0, FINE_Y, H_FINE, H_COARSE, GROWTH);
P.grid.z = -abl_forvibes.graded(0, P.gel(3), 0, FINE_Z, H_FINE, H_COARSE, GROWTH);
P.flow = flowAt(FLOW, 0);
S = abl_forvibes.init(P);
n_modes = numel(P.grid.x)*numel(P.grid.y)*numel(P.grid.z);  % number of heat modes = grid points in the 3-D block
layer = -S.z < ABLATION_DEPTH - 1e-12;                      % z rows inside the spared layer (depth < ABLATION_DEPTH)
T_safe = P.T_abl - MARGIN;                                  % PID 1 set point for the hottest point in the layer  [degC]

% --- 2. history: one entry per controller tick ---
t  = 0:DT:T_END;  nt = numel(t);
H.t = t;
z0 = zeros(1, nt);
H.P_gen = z0;  H.P_ant = z0;                                % generator power, power at the antenna (step ending at t)  [W]
H.P1 = z0;  H.P2 = z0;                                      % what PID 1 and PID 2 asked for  [W]
H.owner = z0;                                               % who set the power: 1 = PID 1, 2 = PID 2, 3 = look-ahead safety back-off
H.flow = z0;  H.y1 = z0;  H.y1pred = z0;                    % flow [L/min]; hottest point in the layer now / in the look-ahead [degC]
H.A = z0;  H.r = z0;  H.t_tick = z0;                        % ablated area [mm^2]; its smoothed growth rate [mm^2/s]; compute time [s]
H.peak = z0;  H.surf = z0;  H.Tc = z0;  H.depth = z0;
H.T = zeros([size(S.T), nt], 'single');
H.flow(1) = P.flow;  H.y1(1) = max(S.T(:, layer), [], 'all');  H.y1pred(1) = H.y1(1);
H.peak(1) = max(S.T(:));  H.surf(1) = max(S.T(:,1));  H.Tc(1) = S.Tc;  H.T(:,:,1) = S.T;

% --- 3. controller state ---
I1 = MAX_WATTS;  I2 = min(max(P_START, P_MIN), MAX_WATTS);  % integral terms; PID 2 starts at P_START  [W]
y1_prev = H.y1(1);  r = 0;  r_prev = 0;  A_prev = 0;
n_done = nt;                                                % last tick simulated (less if stopped early)
if LIVE
    L = liveWindow(S, P, FLOW, ABLATION_DEPTH, T_safe, MAX_WATTS, T_END, SPEED);
    wall = tic;  wall_target = 0;                           % pacing: when each tick is due on the clock  [s]
end

% --- 4. the closed loop ---
for n = 2:nt
    tick = tic;

    % 4a. flow: from the schedule, or from the live window's slider once you have moved it;
    %     switch the model to the new flow if it changed
    if LIVE && isvalid(L.fig) && L.cbManual.Value
        f = round(L.slFlow.Value, 1);
    else
        f = flowAt(FLOW, t(n-1));
    end
    if f ~= S.P.flow, S = abl_forvibes.reflow(S, f); end

    % 4b. measure (perfect map)
    A  = 1e6*sum(S.cell(S.Tmax >= P.T_abl));                    % ablated area on the antenna plane  [mm^2]
    r  = r + DT/TAU_R*((A - A_prev)/DT - r);                % smoothed growth rate  [mm^2/s]

    % 4c. PID 1: headroom below the safe temperature, measured on where the layer is heading
    y1m = lookAhead(S, H.P_gen(n-1), DT, T_LOOK, layer, Inf);    % hottest the layer gets with the last power, then off  [degC]
    e1 = T_safe - y1m;
    I1 = I1 + Ki1*e1*DT;
    u1 = I1 + Kp1*e1 - Kd1*(y1m - y1_prev)/DT;
    P1 = min(max(u1, 0), MAX_WATTS);

    % 4d. PID 2: growth-rate shortfall (or a fixed power in 'hold' mode)
    e2 = R_TARGET - r;
    if strcmp(MODE, 'hold')
        u2 = P_HOLD;                                        %#ok<UNRCH>
    else
        I2 = I2 + Ki2*e2*DT;
        u2 = I2 + Kp2*e2 - Kd2*(r - r_prev)/DT;
    end
    P2 = min(max(u2, P_MIN), MAX_WATTS);

    % 4e. PID 1 has priority: take the lower request
    if P1 < P2, Pg = P1; owner = 1; else, Pg = P2; owner = 2; end

    % 4f. look-ahead safety check: Pg for one step, then off for T_LOOK seconds;
    %     if unsafe, halve the power until it is safe (0 W as the last resort)
    [y1pred, ok] = lookAhead(S, Pg, DT, T_LOOK, layer, P.T_abl - MARGIN/2);
    while ~ok && Pg > 0
        owner = 3;
        if Pg > 1, Pg = Pg/2; else, Pg = 0; end
        [y1pred, ok] = lookAhead(S, Pg, DT, T_LOOK, layer, P.T_abl - MARGIN/2);
    end

    % 4g. anti-windup: the loop not in charge tracks the applied power
    if owner == 2, I1 = Pg; else, I2 = Pg; end
    if owner == 1 && u1 ~= P1, I1 = Pg - Kp1*e1; end         % PID 1 in charge but clamped
    if owner == 3, I1 = Pg - Kp1*e1; end                     % after a safety back-off, PID 1 continues from the safe power
    if owner == 2 && ~strcmp(MODE, 'hold') && u2 ~= P2, I2 = Pg - Kp2*e2; end   % PID 2 in charge but clamped

    % 4h. apply the power for one step (the "patient")
    S = abl_forvibes.step(S, Pg, DT);

    % 4i. record
    H.P_gen(n) = Pg;  H.P_ant(n) = (1 - S.P.drop)*Pg;  H.P1(n) = P1;  H.P2(n) = P2;  H.owner(n) = owner;
    H.flow(n) = S.P.flow;  H.y1(n) = max(S.T(:, layer), [], 'all');  H.y1pred(n) = y1pred;
    H.A(n) = 1e6*sum(S.cell(S.Tmax >= P.T_abl));  H.r(n) = r;
    H.peak(n) = max(S.T(:));  H.surf(n) = max(S.T(:,1));  H.Tc(n) = S.Tc;  H.depth(n) = abl_forvibes.depth(S);
    H.T(:,:,n) = S.T;
    y1_prev = y1m;  r_prev = r;  A_prev = A;
    H.t_tick(n) = toc(tick);

    % 4j. live window: redraw, then wait until this tick is due on the clock
    if LIVE
        if ~isvalid(L.fig) || L.fig.UserData.stop, n_done = n; break, end
        liveUpdate(L, H, n, S, f);
        wall_target = wall_target + DT/L.ddSpeed.Value;
        pause(max(wall_target - toc(wall), 0.001));        % pause also lets the slider and buttons respond
    end
end

% stopped early: keep only the ticks that were simulated
if n_done < nt
    for name = fieldnames(H).'
        v = H.(name{1});
        if isequal(size(v), [1 nt]), H.(name{1}) = v(1:n_done); end
    end
    H.T = H.T(:,:,1:n_done);  T_END = H.t(end);
end
if LIVE && isvalid(L.fig), L.fig.UserData.stop = true; L.bStop.Text = 'Finished'; L.bStop.Enable = 'off'; end

% --- 5. report ---
worst = max(S.Tmax(:, layer), [], 'all');
fprintf('\nabl_pid_viberuntimereduction (%s mode), %g s, ABLATION_DEPTH %.2f mm, MAX_WATTS %g W\n', MODE, T_END, 1e3*ABLATION_DEPTH, MAX_WATTS);
fprintf('  graded grid: %.2f mm fine zone (%g mm deep, %g mm wide), %.2f mm along x, up to %g mm far away; %d modes\n', ...
        1e3*H_FINE, 1e3*FINE_Z, 1e3*FINE_Y, 1e3*H_X, 1e3*H_COARSE, n_modes);
fprintf('  spared layer: hottest ever %.2f degC (T_abl %g) -> %s\n', worst, P.T_abl, tern(worst < P.T_abl, 'SPARED', 'VIOLATED'));
fprintf('  ablated area %.0f mm^2, depth %.2f mm, mean power %.1f W, safety back-offs %d\n', H.A(end), 1e3*H.depth(end), ...
        mean(H.P_gen(2:end)), nnz(H.owner == 3));
fprintf('  in charge: PID 1 %.0f %%, PID 2 %.0f %% of the time; slowest tick %.0f ms (limit %.0f ms)\n', ...
        100*mean(H.owner(2:end) == 1), 100*mean(H.owner(2:end) == 2), 1e3*max(H.t_tick), 1e3*DT);

% --- 6. show it ---
info.name = 'Ablation model: closed-loop PID control';
info.envTitle = sprintf('Generator power set by the PIDs (MAX_WATTS %g W), shown at the antenna', MAX_WATTS);
info.depthLimit = ABLATION_DEPTH;
abl_view(H, S, P, info);
controlFigure(H, P, T_safe, R_TARGET, MAX_WATTS);
end


% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function f = flowAt(FLOW, tt)
% flow rate in force at time tt from the step schedule [t1 f1; t2 f2; ...]
k = find(FLOW(:,1) <= tt + 1e-9, 1, 'last');
if isempty(k), k = 1; end
f = FLOW(k, 2);
end

function [ymax, ok] = lookAhead(S, Pg, DT, T_LOOK, layer, T_lim)
% simulate Pg for one step, then off for T_LOOK s (1 s steps, exact), on a copy
% of the state; ymax = hottest point in the spared layer seen, ok = stays below T_lim
S = abl_forvibes.step(S, Pg, DT);
ymax = max(S.T(:, layer), [], 'all');
for k = 1:ceil(T_LOOK)
    S = abl_forvibes.step(S, 0, 1);
    ymax = max(ymax, max(S.T(:, layer), [], 'all'));
end
ok = ymax < T_lim;
end

function s = tern(c, a, b)
if c, s = a; else, s = b; end
end

function controlFigure(H, P, T_safe, R_TARGET, MAX_WATTS)
% what the controller did: flow, the two PID requests and the applied power,
% the spared layer's temperature, and the ablated area and its growth rate
t = H.t;
f = figure('Name', 'Ablation model: controller', 'Position', [100 100 1100 850]);
tl = tiledlayout(f, 4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(tl);  stairs(ax, t, H.flow, 'LineWidth', 1.6);  grid(ax, 'on');
ylabel(ax, 'flow (L/min)');  title(ax, 'Circulating flow (the disturbance)');
ylim(ax, [0 1.2*max(H.flow)]);

ax = nexttile(tl);  hold(ax, 'on');  grid(ax, 'on');
in1 = double(H.owner == 1 | H.owner == 3)*1.1*MAX_WATTS;
area(ax, t, in1, 'FaceColor', [1 0.3 0.6], 'FaceAlpha', 0.12, 'EdgeColor', 'none', 'DisplayName', 'PID 1 in charge');
stairs(ax, t, min(H.P1, 1.05*MAX_WATTS), ':', 'Color', [0.85 0.2 0.5], 'LineWidth', 1.3, 'DisplayName', 'PID 1 request');
stairs(ax, t, H.P2, ':', 'Color', [0.1 0.55 0.3], 'LineWidth', 1.3, 'DisplayName', 'PID 2 request');
stairs(ax, t, H.P_gen, '-', 'SeriesIndex', 'none', 'LineWidth', 1.8, 'DisplayName', 'applied (generator)');
yline(ax, MAX_WATTS, '--', 'MAX\_WATTS', 'HandleVisibility', 'off');
ylim(ax, [0 1.1*MAX_WATTS]);  ylabel(ax, 'P_{gen} (W)');  title(ax, 'Power: min(PID 1, PID 2)');
legend(ax, 'Location', 'eastoutside');

ax = nexttile(tl);  hold(ax, 'on');  grid(ax, 'on');
plot(ax, t, H.y1, '-', 'Color', [0.85 0.2 0.5], 'LineWidth', 1.6, 'DisplayName', 'hottest in spared layer');
plot(ax, t, H.y1pred, ':', 'Color', [0.85 0.2 0.5], 'LineWidth', 1.2, 'DisplayName', 'look-ahead (power off)');
yline(ax, T_safe, '--', 'PID 1 set point', 'HandleVisibility', 'off');
yline(ax, P.T_abl, '-', 'T_{abl}', 'Color', [0.8 0.1 0.1], 'HandleVisibility', 'off');
ylabel(ax, 'T (°C)');  title(ax, 'Spared surface layer (PID 1)');  legend(ax, 'Location', 'eastoutside');

ax = nexttile(tl);  hold(ax, 'on');  grid(ax, 'on');
yyaxis(ax, 'left');   plot(ax, t, H.A, '-', 'LineWidth', 1.6, 'DisplayName', 'ablated area');  ylabel(ax, 'area (mm^2)');
yyaxis(ax, 'right');  plot(ax, t, H.r, '-', 'LineWidth', 1.2, 'DisplayName', 'growth rate');
yline(ax, R_TARGET, '--', 'R\_TARGET', 'HandleVisibility', 'off');  ylabel(ax, 'growth (mm^2/s)');
xlabel(ax, 'time (s)');  title(ax, 'Ablated region (PID 2)');  legend(ax, 'Location', 'eastoutside');
linkaxes(findall(f, 'Type', 'axes'), 'x');  xlim(ax, [0 t(end)]);
end

function L = liveWindow(S, P, FLOW, ABLATION_DEPTH, T_safe, MAX_WATTS, T_END, SPEED)
% the live window: map on the left; power and spared layer on the right;
% controls on top (time, stop, flow slider, playback speed)
x = S.x*1e3;  z = S.z*1e3;
L.fig = uifigure('Name', 'Ablation model: live PID control', 'Position', [60 60 1400 800]);
L.fig.UserData = struct('stop', false);
L.fig.CloseRequestFcn = @(f,~) set(f, 'UserData', struct('stop', true), 'Visible', 'off');   % closing = stop
g = uigridlayout(L.fig, [3 2], 'RowHeight', {34, '1x', '1x'}, 'ColumnWidth', {'1.25x', '1x'});

% controls
gc = uigridlayout(g, [1 7], 'ColumnWidth', {110, 90, 120, '1x', 90, 60, 110}, 'Padding', [0 0 0 0]);
gc.Layout.Column = [1 2];
L.lblT  = uilabel(gc, 'Text', 't = 0 s', 'FontWeight', 'bold');
L.bStop = uibutton(gc, 'Text', '■ Stop', 'ButtonPushedFcn', @(~,~) set(L.fig, 'UserData', struct('stop', true)));
L.cbManual = uicheckbox(gc, 'Text', 'manual flow', 'Value', false, ...
                        'Tooltip', 'Ticked: the flow follows the slider. Unticked: it follows the FLOW schedule.');
L.slFlow = uislider(gc, 'Limits', [0.5 15], 'Value', flowAt(FLOW, 0), 'MajorTicks', [0.5 2.5 5 7.5 10 12.5 15], 'MinorTicks', [], ...
                    'ValueChangingFcn', @(~,~) set(L.cbManual, 'Value', true));   % touching the slider takes over the flow
L.lblF  = uilabel(gc, 'Text', '', 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
uilabel(gc, 'Text', 'speed', 'HorizontalAlignment', 'right');
L.ddSpeed = uidropdown(gc, 'Items', {'real time', '2x', '5x', '20x'}, 'ItemsData', {1, 2, 5, 20}, ...
                       'Value', 1);
[~, k] = min(abs([1 2 5 20] - SPEED));  L.ddSpeed.Value = L.ddSpeed.ItemsData{k};   % nearest choice to SPEED

% map
L.axMap = uiaxes(g);  L.axMap.Layout.Row = [2 3];  L.axMap.Layout.Column = 1;
ax = L.axMap;  hold(ax, 'on');  ax.Toolbar.Visible = 'off';  disableDefaultInteractivity(ax);
d_mm = P.d*1e3;  rb = P.b*1e3;  Lm = P.L*1e3;  top = d_mm + rb + 1.5;
patch(ax, x([1 end end 1]), [0 0 top top], [0.44 0.67 0.75], 'EdgeColor', 'none');      % water above the gel
L.hImg = surface(ax, x, z, zeros(numel(z), numel(x)), S.T.', 'EdgeColor', 'none', 'FaceColor', 'interp');   % uneven grid, so a surface not an image
ax.SortMethod = 'childorder';                               % draw in the order added, so the lines stay on top
warning('off', 'MATLAB:contour:ConstantData');
[~, L.hIso] = contour(ax, x, z, S.Tmax.', [P.T_abl P.T_abl], '--', 'Color', 'w', 'LineWidth', 1.4);
plot(ax, x([1 end]), -1e3*ABLATION_DEPTH*[1 1], ':', 'Color', [1 0.3 0.6], 'LineWidth', 1.6);
patch(ax, [0 Lm Lm 0], d_mm + [-rb -rb rb rb], [0.86 0.93 0.97], 'EdgeColor', [0.55 0.65 0.7]);   % antenna
set(ax, 'YDir', 'normal', 'Layer', 'top');  axis(ax, 'image');  box(ax, 'on');
xlim(ax, x([1 end]));  ylim(ax, [z(end) top]);
colormap(ax, turbo(256));  clim(ax, [min(P.T_b, P.T_0), P.T_abl + 25]);
cb = colorbar(ax);  cb.Label.String = 'T (°C)';
xlabel(ax, 'x along the antenna (mm), feed at 0');  ylabel(ax, 'z (mm), gel surface at 0');
title(ax, 'white dashed: ablated zone   pink dotted: ABLATION\_DEPTH', 'FontWeight', 'normal');

% power
L.axP = uiaxes(g);  L.axP.Layout.Row = 2;  L.axP.Layout.Column = 2;
ax = L.axP;  hold(ax, 'on');  grid(ax, 'on');  ax.Toolbar.Visible = 'off';
L.hP1 = plot(ax, NaN, NaN, ':', 'Color', [0.85 0.2 0.5], 'LineWidth', 1.3, 'DisplayName', 'PID 1 request');
L.hP2 = plot(ax, NaN, NaN, ':', 'Color', [0.1 0.55 0.3], 'LineWidth', 1.3, 'DisplayName', 'PID 2 request');
L.hPg = plot(ax, NaN, NaN, '-', 'SeriesIndex', 'none', 'LineWidth', 1.8, 'DisplayName', 'applied');
yline(ax, MAX_WATTS, '--', 'MAX\_WATTS', 'HandleVisibility', 'off');
xlim(ax, [0 T_END]);  ylim(ax, [0 1.15*MAX_WATTS]);
xlabel(ax, 'time (s)');  ylabel(ax, 'P_{gen} (W)');  legend(ax, 'Location', 'north', 'NumColumns', 3);
title(ax, '', 'FontWeight', 'normal');

% spared layer and area
L.axL = uiaxes(g);  L.axL.Layout.Row = 3;  L.axL.Layout.Column = 2;
ax = L.axL;  hold(ax, 'on');  grid(ax, 'on');  ax.Toolbar.Visible = 'off';
yyaxis(ax, 'left');
L.hY = plot(ax, NaN, NaN, '-', 'Color', [0.85 0.2 0.5], 'LineWidth', 1.6, 'DisplayName', 'hottest in spared layer');
yline(ax, T_safe, '--', 'set point', 'HandleVisibility', 'off');
yline(ax, P.T_abl, '-', 'T_{abl}', 'Color', [0.8 0.1 0.1], 'HandleVisibility', 'off');
ylim(ax, [min(P.T_b, P.T_0) - 1, P.T_abl + 3]);  ylabel(ax, 'T (°C)');
yyaxis(ax, 'right');
L.hA = plot(ax, NaN, NaN, '-', 'Color', [0.1 0.55 0.3], 'LineWidth', 1.4, 'DisplayName', 'ablated area');
ylabel(ax, 'area (mm^2)');  ax.YColor = [0.1 0.55 0.3];  ylim(ax, [0 Inf]);
xlim(ax, [0 T_END]);  xlabel(ax, 'time (s)');  legend(ax, 'Location', 'northwest');
end

function liveUpdate(L, H, n, S, f)
% redraw the live window for tick n
L.hImg.CData = S.T.';
L.hIso.ZData = S.Tmax.';
k = 1:n;  t = H.t(k);
set(L.hPg, 'XData', t, 'YData', H.P_gen(k));
set(L.hP1, 'XData', t, 'YData', min(H.P1(k), L.axP.YLim(2)));
set(L.hP2, 'XData', t, 'YData', H.P2(k));
set(L.hY,  'XData', t, 'YData', H.y1(k));
set(L.hA,  'XData', t, 'YData', H.A(k));
who = {'PID 1 (protecting the surface layer)', 'PID 2 (growing the zone)', 'safety back-off'};
L.axP.Title.String = sprintf('%.1f W   in charge: %s', H.P_gen(n), who{H.owner(n)});
L.lblT.Text = sprintf('t = %.1f s', H.t(n));
if ~L.cbManual.Value, L.slFlow.Value = min(max(f, L.slFlow.Limits(1)), L.slFlow.Limits(2)); end
L.lblF.Text = sprintf('%.1f L/min', f);
L.axMap.Title.String = sprintf('peak %.1f °C   ablated %.0f mm^2, depth %.2f mm   (white dashed: zone, pink dotted: ABLATION\\_DEPTH)', ...
                               max(S.T(:)), H.A(n), 1e3*H.depth(n));
drawnow limitrate;
end
