function [H, S] = abl_plot(P_gen, duty, t_treat, T_pulse)
% ABL_PLOT  Run the reduced ablation model (abl.m) for a pulsed power and show
% the result in the layout of the V5 GUI, with a time slider to scrub through it.
%
%   abl_plot(P_gen, duty, t_treat)            pulse period 10 s (the V5 default)
%   abl_plot(P_gen, duty, t_treat, T_pulse)   choose the pulse period
%   [H, S] = abl_plot(...)                    also return the history H and final state S
%
% Inputs:   P_gen     $P_{gen}$  generator power while ON  [W]
%           duty      duty cycle: share of each pulse period that the power is ON.
%                     Give it as a fraction (0.5) or a percentage (50); 1 = continuous (CW)
%           t_treat   treatment time  [s]
%           T_pulse   $T_{pulse}$  pulse period: ON for duty*T_pulse, then OFF  [s]
%
% Example:  abl_plot(60, 50, 300)     60 W, 50 % duty (5 s on / 5 s off), 5 minutes
%
% The window (drawn by abl_view.m):
%   left, top      temperature on the antenna plane at the chosen time, with the
%                  ablated zone (everywhere that has reached $T_{abl}$ so far) outlined in white
%   left, middle   play / reset / time slider / playback speed
%   left, bottom   temperatures and ablation depth over time
%   right, top     the power reaching the antenna over time
%   right, bottom  temperature against depth through the hottest column
% All parameters come from abl.params, so edit abl.m to change the model.

if nargin < 4 || isempty(T_pulse), T_pulse = 10; end
if duty > 1, duty = duty/100; end           % a percentage was given
duty = min(max(duty, 0), 1);

% --- 1. set the model up ---
P = abl.params;
S = abl.init(P);

% --- 2. the power envelope: when the power is ON and OFF ---
% Each pulse period is ON for duty*T_pulse, then OFF for the rest. abl.step is
% exact while the power is held constant, so the simulation steps exactly at
% every switch. Long ON/OFF stretches are cut into pieces of at most dt_max so
% the history curves are smooth and the slider has a frame at least every dt_max.
dt_max = 1;                                                 % longest step between recorded frames  [s]
t_on   = duty*T_pulse;                                      % ON time in each period  [s]
edges  = unique([0:T_pulse:t_treat, (0:T_pulse:t_treat) + t_on, t_treat]);   % every ON/OFF switch  [s]
edges  = edges(edges >= 0 & edges <= t_treat);
t = 0;
for k = 1:numel(edges) - 1                                  % cut each stretch into steps of <= dt_max
    m = ceil((edges(k+1) - edges(k))/dt_max - 1e-9);
    t = [t, edges(k) + (1:m)*(edges(k+1) - edges(k))/m];   %#ok<AGROW>
end
on = @(tt) mod(tt, T_pulse) < t_on - 1e-9;                  % is the power ON at time tt?

% --- 3. run the simulation and record every frame ---
% The temperature on the antenna plane is stored at every step (as single
% precision, to halve the memory) so the slider can show any moment.
nt = numel(t);
H.t     = t;                                                % time  [s]
H.P_gen = zeros(1, nt);                                     % generator power during the step ending at t  [W]
H.P_ant = zeros(1, nt);                                     % power reaching the antenna during the step ending at t  [W]
H.peak  = zeros(1, nt);                                     % hottest point on the antenna plane  [degC]
H.surf  = zeros(1, nt);                                     % hottest point on the gel surface  [degC]
H.Tc    = zeros(1, nt);                                     % catheter tip  [degC]
H.depth = zeros(1, nt);                                     % ablation depth  [m]
H.T     = zeros([size(S.T), nt], 'single');                 % temperature on the antenna plane, every frame  [degC]
H.peak(1) = max(S.T(:));  H.surf(1) = max(S.T(:,1));  H.Tc(1) = S.Tc;  H.T(:,:,1) = S.T;
for k = 2:nt
    Pg = P_gen*on((t(k-1) + t(k))/2);                       % power during this step (constant within it)
    S  = abl.step(S, Pg, t(k) - t(k-1));
    H.P_gen(k) = Pg;
    H.P_ant(k) = (1 - P.drop)*Pg;
    H.peak(k)  = max(S.T(:));
    H.surf(k)  = max(S.T(:,1));                             % column 1 of S.T is the surface, z = 0
    H.Tc(k)    = S.Tc;
    H.depth(k) = abl.depth(S);
    H.T(:,:,k) = S.T;
end
P_del = (1 - P.drop)*P_gen;                                 % $P_{del}$  power at the antenna while ON  [W]

% --- 4. show it ---
info.envTitle = sprintf('Power envelope: %.1f W when on, duty %.0f %% (%g s period)', P_del, 100*duty, T_pulse);
abl_view(H, S, P, info);
end
