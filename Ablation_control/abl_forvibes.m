classdef abl_forvibes
% This is a minimal and debloated version of the heating simulation ran in the v5 prototype.
% abl_forvibes = abl.m plus graded grids (P.grid, abl_forvibes.graded): fine spacing only
% near the antenna and the surface, for finer resolution at less compute. Used by
% abl_pid_viberuntimereduction.m. Without P.grid it behaves exactly like abl.m.
% Using this reduced model we have faster simulation times, which is useful for the PID design and testing. 
%
% These are the variables we need to keep track.
%     power     $P_{del} = (1 - \delta) P_{gen}$,   $P_{lost} = \delta P_{gen}$
%     source    $Q(\mathbf{r}) = \sigma E_{rms}^2(\mathbf{r}) P_{del}$      ($E_{rms}$: King's insulated monopole at 1 W)
%     gel       $\rho c \, \partial T / \partial t = k \nabla^2 T + Q$          (3-D, symmetric about $y = 0$)
%     surface   $-k \, \partial T / \partial z = h_s (T - T_w)$,   $T_w = T_b + \alpha (T_{cat} - T_b)$  on $0 \le x \le L_c$
%               $h_s = h_{ref} (f_w / f_{ref})^n$        ($n = 1/2$: laminar flat plate, $Nu \propto Re^{1/2}$)
%     walls     $-k \, \partial T / \partial n = h_w (T - T_b)$          (surroundings held at body temperature)
%     catheter  $m_{cat} c_{cat} \, dT_{cat}/dt = \gamma P_{lost} - k_{cat} (T_{cat} - T_b)$
%     transfer  $\alpha = \min(1, \beta / (d \sqrt{f_w}))$
%     ablation  $T_{max}(\mathbf{r}) = \max_t T(\mathbf{r}, t)$,   ablated where $T_{max} \ge T_{abl}$
%
%   P = abl_forvibes.params;  S = abl_forvibes.init(P);         % once per geometry (~1 s)
%   for n = 1:N, S = abl_forvibes.step(S, P_gen(n), dt); end    % exact for P_gen held over dt
%   abl_forvibes.depth(S)                               % ablation depth (m)
%
%   Controller states (struct S):
%     S.T    $T(x, z)$ on the antenna plane $y = 0$  [degC], nx x nz
%     S.Tc   $T_{cat}$  catheter tip temperature  [degC]
%     S.Tmax $T_{max}(x, z)$  peak temperature so far  [degC]
%     S.t    $t$  elapsed time  [s]
%   Coordinates: $x$ along the element (feed at 0), $y$ across, $z$ up (gel surface at 0, gel at $z < 0$).

methods (Static)

function P = params()
% table preset of v5_defaults; SI units, temperatures in degC
% %&% denotes a changing variable acting as input to alforithm (want to make these as few as possible)
% %*% denotes a target variable that is set before procedure, e.g. ablation cutoff temp at 54 but set to 50 for safety 
% %$% denotes a changeABLE variable (whats effected by the controller)
% %#% denotes a weird variable that probably needs an exact measurement to be derived
% %^% denotes my cool experimental variables
 
% --- antenna ---

P.f       = 2.45e9;          % $f$ microwave frequency  [Hz]
P.L       = 27e-3;           % $L$ element length, feed to tip  [m]
P.a       = 0.27e-3;         % $a$ inner conductor radius  [m]
P.b       = 0.9e-3;          % $b$ outer radius of the insulation  [m]
P.eps_ins = 2.1;             % insulation (PTFE) relative permittivity  [-]
P.l_f     = 1.5e-3;          % %#% exposed feed block behind the feed (return charge)  [m]
P.d       = 2e-3;            % %&% element axis to gel surface, this changes as the cathetar moves [m]
% --- media ---

P.eps_g   = 43.98 - 13.5j;   % %#% gel complex permittivity  [-]
P.eps_w   = 73 - 8.8j;       % %#% water complex permittivity  [-]
P.rho     = 1000;            % %#% gel density  [kg/m^3]
P.c       = 4180;            % %#% gel specific heat  [J/(kg K)]
P.k       = 0.6;             % %#% gel thermal conductivity  [W/(m K)]
% --- temperatures and boundaries ---

P.T_b     = 37;              % $T_b$ circulating water (blood) temperature  [degC]
P.T_0     = 37;              % $T_0$ gel temperature at $t = 0$  [degC]
P.h_s     = 1400;            % $h_{ref}$ water-to-gel surface heat transfer coefficient at $f_{ref}$ (laminar plate, V4)  [W/(m^2 K)]
P.flow_ref = 5;              % $f_{ref}$ flow $h_{ref}$ was derived at  [L/min]
P.n_h     = 0.5;             % $n\$ turbulance factor 0 = laminar, 1 = turbulent, 0.5 = flat plate  [-]
P.h_w     = 8;               % $h_w$ [MINOR THIS IS MOSTLY FOR THEORETICAL APPLICATION} walls to surroundings at $T_b$ heat transfer coefficient  [W/(m^2 K)]
% --- surface transfer from the catheter ---

P.L_c     = 27e-3;           % $L_c$ length (heated catheter footprint on the surface)  [m]
P.beta    = 9e-4;            % %^% transfer constant in $\alpha$  [m (L/min)^{1/2}]
P.flow    = 5;               % %&%  circulation flow rate  [L/min]
% --- power split and catheter ---

P.drop    = 0.7;             % %#% fraction of energy lost before the antenna  [-]
P.gamma   = 0.043;           % %#% this is the portion of P.drop that just heats the catheter tip (the antenna), NOT the heating due to MW[-]
P.m_cat   = 2e-3;            % %#% tip mass [kg]
P.c_cat   = 830;             % %^% Specific heat capacity of the catheter tip  (experimental average) [J/(kg K)]
P.k_cat   = 0.1;             % %#% Tip-to-water heat loss conductance  [W/K]

% --- ablation ---
P.T_abl   = 50;              % %*% Ablation cutoff temperature  [degC]
% --- numerics ---
P.gel     = [36 18 39]*1e-3; % $[L_x, W_y/2, H_z]$   gel block: length, half-width, depth  [m]
P.x0      = -4.5e-3;         % $x_0$           block edge behind the feed (feed at $x = 0$)  [m]
P.dh      = 1e-3;          % $\Delta h$      grid spacing  [m]
end


function S = init(P)
% Builds everything that does not change while the controller is running:
% the grid of points in the gel, the "heat modes" of the gel block, and how
% strongly each heat source (antenna, warm water) drives each mode.
%
% Inputs:   P   parameter struct from abl_forvibes.params
% Outputs:  S   the simulation state, handed to abl_forvibes.step from then on
%
% The big idea: modes
%   Any temperature pattern in the block can be written as a sum of fixed
%   shapes, the "modes" $\phi_m(\mathbf{r})$, much like a musical chord is a sum of pure
%   notes. The useful thing about these shapes is that each one simply fades
%   away at its own rate $\lambda_m$ and never turns into a different shape. So instead
%   of tracking the temperature at every point, we track one number per mode,
%   the amount $c_m$ of that shape that is present:
%       $\theta(\mathbf{r}, t) = T - T_b = \sum_m c_m(t) \, \phi_m(\mathbf{r})$
%   Each amount obeys its own simple equation
%       $dc_m/dt = -\lambda_m c_m + g_m$
%   where $g_m$ is how hard the heat sources are pushing that mode. This
%   equation can be solved exactly, so abl_forvibes.step needs no small time steps.

% --- 1. the grid ---
% Normally the points are evenly spaced P.dh apart. A graded grid (P.grid, made
% by abl_forvibes.graded) can be given instead: small spacing where the gel heats up and
% the temperature changes fast, large spacing in the far parts of the block.
% Then dx, dy, dz are the gaps between neighbouring points, not one number.
if isfield(P, 'grid') && ~isempty(P.grid)
    x  = P.grid.x(:).';  y = P.grid.y(:).';  z = P.grid.z(:).';   % positions given by the graded grid  [m]
    dx = diff(x);  dy = diff(y);  dz = -diff(z);    % gaps between neighbouring points  [m]
    n  = [numel(x), numel(y), numel(z)];
else
    dh = P.dh;                       % $\Delta h$   grid spacing 0.5 slow, 0.1 slower[m]
    n  = round(P.gel/dh) + 1;        % $[n_x, n_y, n_z]$   number of grid points along x, y, z
    x  = P.x0 + (0:n(1)-1)*dh;       % $x$   positions along the antenna, feed at x = 0  [m]
    y  = (0:n(2)-1)*dh;              % $y$   positions across, y = 0 is the antenna plane  [m]
    z  = -(0:n(3)-1)*dh;             % $z$   depths, z = 0 is the gel surface, negative into the gel  [m]
    dx = dh;  dy = dh;  dz = dh;     % one spacing for every gap
end
rhoc = P.rho*P.c;                    % $\rho c$   heat needed to warm 1 m^3 of gel by 1 K  [J/(m^3 K)]

% --- 2. surface cooling at the current flow ---
h_s = P.h_s*(P.flow/P.flow_ref)^P.n_h;   % $h_s = h_{ref} (f_w / f_{ref})^n$  [W/(m^2 K)]

% --- 3. the heat modes of the block ---
% Because the block is a rectangular box, every 3-D mode is just a 1-D shape
% along x times a 1-D shape along y times a 1-D shape along z. So we only
% have to find the 1-D shapes (abl_forvibes.modes), which is fast.
% Outputs per direction: V = the 1-D shapes (one per column), lam = how fast
% each fades [1/m^2], w = the length of gel each grid point stands for [m].
[Vx, lam_x, wx] = abl_forvibes.modes(n(1), dx, P.h_w/P.k, P.h_w/P.k);
[Vy, lam_y, wy] = abl_forvibes.modes(n(2), dy, 0,         P.h_w/P.k);
[Vz, lam_z, wz] = abl_forvibes.modes(n(3), dz, h_s/P.k,   P.h_w/P.k);

% builds the 3-D array a(i) * b(j) * c(k)
outer3 = @(a, b, c) a(:).*b(:).'.*reshape(c, 1, 1, []);

% --- 4. the microwave heating in the gel ---
% The antenna is a straight wire, so its field only depends on how far along
% the wire you are (x) and how far you are from the wire's axis (r). We work the
% field out on a list of 300 distances r, spaced closer together near the
% wire where the field changes fastest, and then look up every grid point in
% that table. The look-up is done on log(E), because E falls off very steeply
% with distance and log(E) is much smoother to interpolate.
r_min  = P.d;                                       % closest gel point: the surface right under the wire  [m]
r_max  = 1.001*hypot(y(end), P.d - z(end));         % furthest gel point: the bottom far corner  [m]
r      = exp(linspace(log(r_min), log(r_max), 300));    % 300 distances from the axis, log-spaced  [m]
logE   = griddedInterpolant({log(r), x}, log(abl_forvibes.erms(P, x, r)));   % look-up table of log(E_rms) at 1 W
[X, Y, Z] = ndgrid(x, y, z);                        % x, y, z of every grid point
R      = hypot(Y, P.d - Z);                         % distance of every grid point from the antenna axis  [m]
Q      = abl_forvibes.sigma(P)*exp(2*logE(log(R), X));       % $Q = \sigma E_{rms}^2$  heating per W delivered  [W/m^3 per W]

% --- 5. the catheter footprint on the surface ---
% The warm water from the catheter only lies over 0 <= x <= L_c. Each grid
% point along x stands for a strip of surface reaching halfway to the points
% on either side (cut short at the ends of the block). fp is how much of each
% strip is under the catheter: the whole strip if fully under it, 0 if not,
% something in between at the edges.
gx = dx.*ones(1, n(1) - 1);                         % every gap along x (the same number when evenly spaced)  [m]
half_lo  = [gx(1), gx]/2;  half_hi = [gx, gx(end)]/2;   % half the gap to the left and right neighbour  [m]
strip_lo = max(x - half_lo, x(1));                  % start of each grid point's strip  [m]
strip_hi = min(x + half_hi, x(end));                % end of each grid point's strip  [m]
fp = max(0, min(strip_hi, P.L_c) - max(strip_lo, 0));   % length of each strip under the catheter  [m]

% --- 6. how hard each heat source pushes each mode ---
% Every source is converted from "heat at each point" into "push on each mode"
% (abl_forvibes.proj does this conversion). Dividing by $\rho c$ turns watts into kelvin per second.
% Both are worked out for a unit source, and abl_forvibes.step scales them each step:
%   $g_m = P_{del} \, g_{Q,m} + \Theta \, g_{U,m}$,   $\Theta = \alpha (T_{cat} - T_b)$
S.gQ = abl_forvibes.proj(Q.*outer3(wx, wy, wz), Vx, Vy, Vz)/rhoc;     % $g_Q$  push from the microwaves, per W delivered
S.gU = outer3(Vx.'*fp(:), Vy.'*wy, h_s*Vz(1,:))/rhoc;          % $g_U$  push from surface water 1 K warmer than $T_b$ under the catheter

% --- 7. starting amounts and fading rates ---
% At t = 0 the whole gel sits at $T_0$, so $\theta(0) = T_0 - T_b$ everywhere.
% The fading rate of a 3-D mode is the sum of its three 1-D rates times the
% gel's thermal diffusivity $\kappa = k / \rho c$ (how fast heat spreads through it).
S.c   = (P.T_0 - P.T_b)*outer3(Vx.'*wx, Vy.'*wy, Vz.'*wz);           % $c_m(0)$  starting amount of each mode  [K m^{3/2}]
S.Lam = P.k/rhoc*(lam_x + lam_y.' + reshape(lam_z, 1, 1, []));       % $\lambda_m = \kappa (\lambda_x + \lambda_y + \lambda_z)$  fading rates  [1/s]

% --- 8. catheter constants ---
S.alpha = min(1, P.beta/(P.d*sqrt(P.flow)));    % $\alpha$  share of the tip's temperature rise passed to the surface water  [-]
S.tau   = P.m_cat*P.c_cat/P.k_cat;              % $\tau_{cat} = m_{cat} c_{cat} / k_{cat}$  how long the tip takes to heat up  [s]

% --- 9. store what abl_forvibes.step needs and set the starting state ---
% Only the y = 0 row of the y-shapes is kept, since the output is the antenna plane.
S.Vx = Vx;  S.Vz = Vz;  S.vy = Vy(1,:);         % mode shapes along x, z, and their values at y = 0
S.P = P;  S.x = x;  S.z = z;  S.h_s = h_s;      % parameters, grid, surface coefficient actually used
% Kept for abl_forvibes.reflow: a flow change only alters the z-shapes, so the x and y
% parts of the fading rates and of the warm-water push can be reused as they are.
S.wz   = wz;                                    % length of gel each z point stands for  [m]
S.dz   = dz;                                    % gaps between the z points (one number if evenly spaced)  [m]
S.cell = wx*wz.';                               % area of antenna plane each (x, z) point stands for, used to measure the ablated area  [m^2]
S.Lxy  = P.k/rhoc*(lam_x + lam_y.');            % $\kappa (\lambda_x + \lambda_y)$  x and y part of the fading rates  [1/s]
S.gUxy = (Vx.'*fp(:)).*(Vy.'*wy).'/rhoc;        % x and y part of $g_U$ (multiplied by $h_s$ times the z-shape at the surface)
S.t    = 0;                                    % $t$  elapsed time  [s]
S.Tc   = P.T_b;                                 % $T_{cat}$  the tip starts at water temperature  [degC]
S.T    = abl_forvibes.plane(S);                          % $T(x, z)$  starting temperature on the antenna plane  [degC]
S.Tmax = S.T;                                   % $T_{max}(x, z)$  hottest temperature reached so far  [degC]
end


% =========================================================================
% STEP: move the simulation forward by one time step
% =========================================================================
function S = step(S, Pg, dt)
% Advances the simulation by dt seconds with the generator held at Pg watts.
% This is what the controller calls every time it picks a new power.
% Both the catheter tip and the gel are solved exactly over the step (no
% approximation in time), as long as Pg stays constant during the step.
%
% Inputs:   S    state from abl_forvibes.init or the previous abl_forvibes.step
%           Pg   $P_{gen}$  generator power for this step  [W]
%           dt   $\Delta t$  step length  [s]
% Outputs:  S    updated state (S.T, S.Tc, S.Tmax, S.t)
P = S.P;

% --- 1. catheter tip temperature ---
% The tip is heated by its share of the lost power and loses heat to the water:
%   $m_{cat} c_{cat} \, dT_{cat}/dt = \gamma \delta P_{gen} - k_{cat} (T_{cat} - T_b)$
% Over one step at constant power, the tip moves from where it was toward its
% steady temperature $T_b + \gamma \delta P_{gen} / k_{cat}$, closing the gap by the factor $1 - e^{-\Delta t/\tau_{cat}}$:
%   $T_{cat}^+ = T_b + (T_{cat} - T_b) e^{-\Delta t/\tau_{cat}} + \frac{\gamma \delta P_{gen}}{k_{cat}} (1 - e^{-\Delta t/\tau_{cat}})$
fade   = exp(-dt/S.tau);                                % $e^{-\Delta t/\tau_{cat}}$  how much of the old temperature rise survives the step
P_tip  = P.gamma*P.drop*Pg;                             % $\gamma \delta P_{gen}$  power heating the tip  [W]
Tc_new = P.T_b + (S.Tc - P.T_b)*fade + P_tip/P.k_cat*(1 - fade);   % $T_{cat}^+$  tip temperature at the end of the step  [degC]

% --- 2. warm water over the gel surface ---
% The water passing the tip warms up by a share $\alpha$ of the tip's rise. The
% tip changes during the step, so the average of its start and end values is used.
Theta = S.alpha*((S.Tc + Tc_new)/2 - P.T_b);           % $\Theta = \alpha (T_{cat} - T_b)$  surface water rise, step average  [K]

% --- 3. gel temperature (the modes) ---
% Each mode amount obeys $dc_m/dt = -\lambda_m c_m + g_m$ with a constant push $g_m$ over
% the step. Its exact solution: the old amount fades, and the push adds a new part:
%   $c_m^+ = e^{-\lambda_m \Delta t} c_m + \frac{1 - e^{-\lambda_m \Delta t}}{\lambda_m} g_m$
% (expm1(u) = e^u - 1, computed accurately even for tiny u.)
P_del = (1 - P.drop)*Pg;                                % $P_{del} = (1 - \delta) P_{gen}$  power reaching the antenna  [W]
push  = P_del*S.gQ + Theta*S.gU;                        % $g_m$  push on each mode from the microwaves and the warm water
decay = exp(-S.Lam*dt);                                 % $e^{-\lambda_m \Delta t}$  share of each mode that survives the step
gain  = -expm1(-S.Lam*dt)./S.Lam;                       % $(1 - e^{-\lambda_m \Delta t}) / \lambda_m$  how much of the push builds up  [s]
S.c   = decay.*S.c + gain.*push;                        % $c_m^+$  new mode amounts

% --- 4. read out the results ---
S.T    = abl_forvibes.plane(S);                                  % $T(x, z)$  new temperature on the antenna plane  [degC]
S.Tmax = max(S.Tmax, S.T);                              % $T_{max}$  keep the hottest value each point has reached
S.Tc   = Tc_new;
S.t    = S.t + dt;
end


% =========================================================================
% REFLOW: change the circulating flow in the middle of a run
% =========================================================================
function S = reflow(S, flow)
% Switches the simulation to a new flow rate without losing the heat that is
% already in the gel. Running abl_forvibes.init again would start from scratch, so this
% is what a controller calls when the flow changes during the procedure.
%
% Inputs:   S      current state
%           flow   $f_w$  new circulation flow rate  [L/min]
% Outputs:  S      the same moment in time, now cooled by the new flow
%
% What the flow changes:
%   - $h_s$, how hard the water cools the gel surface ($h_s = h_{ref} (f_w / f_{ref})^n$)
%   - $\alpha$, the share of the tip's temperature rise passed to the surface water
% Only the top surface feels $h_s$, so only the z-shapes change; the x and y
% shapes, and the microwave heating pattern, stay exactly as they were.
P = S.P;  P.flow = flow;
nz = numel(S.z);

% --- 1. new z-shapes for the new surface cooling ---
h_s = P.h_s*(P.flow/P.flow_ref)^P.n_h;                  % $h_s$ at the new flow  [W/(m^2 K)]
[Vz, lam_z] = abl_forvibes.modes(nz, S.dz, h_s/P.k, P.h_w/P.k);   % same call as in abl_forvibes.init, new top condition

% --- 2. carry the heat across to the new shapes ---
% The old z-shapes and the new ones both describe the same temperatures, just
% with different building blocks. M converts amounts of old shapes into amounts
% of new shapes (a weighted overlap of every new shape with every old one):
%   $M = V_{z,new}' \, \mathrm{diag}(w_z) \, V_{z,old}$
% Applying M along z to the mode amounts keeps the temperature field identical;
% the microwave push $g_Q$ is converted the same way.
M = Vz.'*(S.wz.*S.Vz);
alongz = @(A) reshape(reshape(A, [], nz)*M.', size(A));  % apply M to the z index of a 3-D array
S.c  = alongz(S.c);                                     % $c_m$  same heat, new building blocks
S.gQ = alongz(S.gQ);                                    % $g_Q$  same microwave heating, new building blocks

% --- 3. rebuild what depends on the new shapes and flow ---
S.Lam   = S.Lxy + reshape(P.k/(P.rho*P.c)*lam_z, 1, 1, []);        % $\lambda_m$  new fading rates  [1/s]
S.gU    = S.gUxy.*reshape(h_s*Vz(1,:), 1, 1, []);                   % $g_U$  warm-water push through the new surface condition
S.alpha = min(1, P.beta/(P.d*sqrt(P.flow)));                       % $\alpha$  at the new flow  [-]
S.Vz = Vz;  S.h_s = h_s;  S.P = P;
end


% =========================================================================
% DEPTH: how deep the ablation has reached
% =========================================================================
function d = depth(S)
% Finds the deepest point on the antenna plane that has ever reached the
% ablation temperature $T_{abl}$.
%
% Inputs:   S   simulation state
% Outputs:  d   ablation depth below the gel surface  [m]
%
% How: for every position x, find the deepest grid point that has been at or
% above $T_{abl}$. The true edge lies somewhere between that point and the one
% below it, so we draw a straight line between their temperatures and find
% where it crosses $T_{abl}$. This gives a depth finer than the grid spacing.
excess = S.Tmax - S.P.T_abl;        % how far each point's peak got above $T_{abl}$; >= 0 means ablated  [K]
nz = numel(S.z);                    % number of grid points in depth
d  = 0;                             % nothing ablated yet
for i = find(any(excess >= 0, 2)).'                 % every x position with at least one ablated point
    j = find(excess(i,:) >= 0, 1, 'last');          % the deepest ablated point at this x
    if j == nz                                      % ablation has reached the bottom of the block:
        d = -S.z(nz);  return                       % report the block depth (make P.gel deeper)
    end
    frac = excess(i,j)/(excess(i,j) - excess(i,j+1));   % how far down from point j to point j+1 the edge sits, 0 to 1
    d = max(d, -S.z(j) + (S.z(j) - S.z(j+1))*frac); % keep the deepest edge over all x
end
end


% =========================================================================
% SIGMA: electrical conductivity of the gel
% =========================================================================
function s = sigma(P)
% How strongly the gel turns the microwave field into heat. It comes from
% the imaginary (lossy) part of the gel's permittivity:
%   $\sigma = \omega \varepsilon_0 \varepsilon_g''$,   $\omega = 2 \pi f$
% P.eps_g is stored as $\varepsilon_g' - j\varepsilon_g''$, so imag(P.eps_g) = $-\varepsilon_g''$, hence the minus sign.
eps0 = 8.8541878128e-12;            % $\varepsilon_0$  permittivity of free space  [F/m]
s = -2*pi*P.f*eps0*imag(P.eps_g);   % $\sigma$  [S/m]
end


% =========================================================================
% ERMS: the antenna's electric field in the gel
% =========================================================================
function E = erms(P, x, r)
% The electric field strength $E_{rms}$ around the antenna when 1 W is delivered
% to it. The heating is $\sigma E_{rms}^2$, so this sets where the gel gets hot.
% Same physics as V5_convection/Insulated_monopole_Erms (see the derivation PDF there).
%
% Inputs:   P   parameters
%           x   positions along the antenna, feed at x = 0  [m]
%           r   distances from the antenna's axis  [m]
% Outputs:  E   $E_{rms}$  [V/m], one row per r and one column per x
%
% Steps:
%   1. how the current wave travels along an insulated wire in water (King's formula)
%   2. how much current at the feed delivers 1 W
%   3. the current and electric charge at each point along the antenna
%   4. add up the field from every small piece of the antenna

% --- physical constants ---
eps0 = 8.8541878128e-12;            % $\varepsilon_0$  permittivity of free space  [F/m]
mu0  = 4e-7*pi;                     % $\mu_0$  permeability of free space  [H/m]
c0   = 299792458;                   % $c_0$  speed of light  [m/s]

% --- 1. wave number along the insulated wire ---
% A wave number k tells how fast a wave's phase turns per metre: wavelength = 2 pi / real(k).
% Inside the water-cooled, insulated wire the current wave travels slower than
% in free space and fades a little; King, Trembly & Strohbehn (1983) give
%   $k_L = k_0 \sqrt{\varepsilon_{ins} \left(1 + \frac{H_0^{(2)}(k_w b)}{k_w b \, H_1^{(2)}(k_w b) \ln(b/a)}\right)}$
% with $H_0^{(2)}, H_1^{(2)}$ Hankel functions (Bessel-type functions for outgoing waves).
w    = 2*pi*P.f;                    % $\omega$  angular frequency  [rad/s]
k0   = w/c0;                        % $k_0$  wave number in free space  [1/m]
lg   = log(P.b/P.a);                % $\ln(b/a)$  appears in every coax formula  [-]
kb   = k0*sqrt(P.eps_w)*P.b;        % $k_w b$  water wave number times insulation radius  [-]
kg   = k0*sqrt(P.eps_g);            % $k_g$  wave number in the gel  [1/m]
kL   = k0*sqrt(P.eps_ins*(1 + besselh(0,2,kb)/(kb*besselh(1,2,kb)*lg)));   % $k_L$  wave number on the antenna  [1/m]

% --- 2. feed current for 1 W delivered ---
% The antenna behaves like a coax line left open at the tip. Its input
% impedance (the ratio of voltage to current at the feed) is
%   $Z_{in} = -j \frac{k_L}{\omega C_{ins}} \cot(k_L L)$,   $C_{ins} = \frac{2 \pi \varepsilon_0 \varepsilon_{ins}}{\ln(b/a)}$
% Power delivered = $\frac{1}{2} |I_0|^2 \mathrm{Re}(Z_{in})$, so 1 W needs $I_0 = \sqrt{2 / \mathrm{Re}(Z_{in})}$.
Z_in = -1j*kL*lg/(2*pi*eps0*P.eps_ins*w)*cot(kL*P.L);  % $Z_{in}$  input impedance  [ohm]
I0   = sqrt(2/real(Z_in));                              % $I_0$  feed current for 1 W  [A]

% --- 3. current and charge along the antenna ---
% The antenna is chopped into short pieces, each treated as a tiny source.
% Element (0 < x < L): a standing wave that is largest at the feed and zero at the tip,
%   $I(x) = I_0 \sin(k_L (L - x)) / \sin(k_L L)$
% Feed block (-l_f < x < 0): the return current falls in a straight line to zero,
%   $I(x) = I_0 (1 + x / l_f)$
% Charge per metre follows from how fast the current changes along the wire,
%   $q = \frac{j}{\omega} \frac{dI}{dx}$
% I and q below are already multiplied by each piece's length.
ds = min(0.05e-3, P.L/200);                 % target piece length  [m]
nb = max(2, round(P.l_f/ds));               % number of pieces on the feed block
ne = max(20, round(P.L/ds));                % number of pieces on the element
xb = ((1:nb) - 0.5)*P.l_f/nb - P.l_f;       % centres of the feed block pieces  [m]
xe = ((1:ne) - 0.5)*P.L/ne;                 % centres of the element pieces  [m]
xs = [xb xe];                               % centres of all pieces, back of the block to the tip  [m]
I_block = I0*(1 + xb/P.l_f)*P.l_f/nb;                           % current x length, feed block  [A m]
I_elem  = I0*sin(kL*(P.L - xe))/sin(kL*P.L)*P.L/ne;             % current x length, element  [A m]
q_block = (1j/w)*(I0/nb)*ones(1, nb);                           % charge per piece, feed block  [C]
q_elem  = (1j/w)*(-I0*kL*cos(kL*(P.L - xe))/sin(kL*P.L)*P.L/ne);    % charge per piece, element  [C]
I = [I_block, I_elem];
q = [q_block, q_elem];

% --- 4. add up the field from every piece ---
% Each piece makes a field in two ways:
%   - its current, through the magnetic vector potential:  $-j \omega A$
%   - its charge, through the electric potential:          $-\nabla \phi$
% Both spread out like $G = e^{-j k_g D} / D$ with distance D from the piece (an
% outgoing wave that fades in the lossy gel). The charges sit in the water just
% above the gel, so their field entering the gel is scaled by $2 \varepsilon_g / (\varepsilon_w + \varepsilon_g)$.
cA = -1j*w*mu0/(4*pi);                          % factor for the current's field
cQ = -2/(4*pi*eps0*(P.eps_w + P.eps_g));        % factor for the charge's field (includes the water-to-gel scaling)
[X, R] = meshgrid(x, r);                        % every (x, r) point where E is wanted
Ex = 0;  Er = 0;                                % field along the wire and away from it  [V/m]
for m = 1:numel(xs)
    DX = X - xs(m);                             % distance along the wire from piece m  [m]
    D  = hypot(DX, R);                          % straight-line distance from piece m  [m]
    G  = exp(-1j*kg*D)./D;                      % $G$  how the wave spreads out from the piece
    H  = -(1 + 1j*kg*D).*G./D.^2;               % $(dG/dD) / D$  used for the charge's field
    Ex = Ex + cA*I(m)*G + cQ*q(m)*H.*DX;
    Er = Er + cQ*q(m)*H.*R;
end
E = sqrt((abs(Ex).^2 + abs(Er).^2)/2);          % $E_{rms}$  (the 1/2 converts peak to rms)  [V/m]
end


% =========================================================================
% MODES: the heat shapes along one direction
% =========================================================================
function [V, lam, w] = modes(n, ds, h0, h1)
% Finds the 1-D heat shapes for a line of n grid points spaced ds apart. The
% 3-D modes in abl_forvibes.init are products of these.
%
% Inputs:   n        number of grid points
%           ds       grid spacing  [m]; or a list of the n-1 gaps between
%                    neighbouring points, for an unevenly spaced (graded) line
%           h0, h1   surface heat loss at the first and last point, as $h/k$  [1/m]
%                    (0 means no heat crosses that end)
% Outputs:  V        the shapes, one per column
%           lam      how fast each shape fades (multiplied by $\kappa$ in abl_forvibes.init)  [1/m^2]
%           w        the length of line each point stands for  [m]
%
% How: heat flows between neighbouring points in proportion to their
% temperature difference (finite volumes). Writing that out for every point gives a matrix K,
% and the shapes that keep their form while fading are the eigenvectors of K
% (weighted by w); the eigenvalues are the fading rates. Scaling makes the shapes
% orthonormal ($V' \mathrm{diag}(w) V = I$), so the amount of each shape in any pattern is a simple weighted sum.
if isscalar(ds)                                 % evenly spaced
    w = ds*ones(n, 1);  w([1 n]) = ds/2;        % each point stands for one spacing; the two end points for half
    K = toeplitz([2 -1 zeros(1, n-2)])/ds;      % heat flow between each point and its two neighbours
    K(1,1) = 1/ds + h0;                         % first point: one neighbour, plus heat lost through its surface
    K(n,n) = 1/ds + h1;                         % last point: same
else                                            % graded: every gap can differ
    % Each point stands for half the gap on its left plus half the gap on its
    % right. The heat flowing across a gap is the temperature difference
    % divided by that gap, so a small gap conducts strongly and a large one weakly.
    ds = ds(:);
    w  = ([ds; 0] + [0; ds])/2;                 % length of line each point stands for  [m]
    g  = 1./ds;                                 % how well each gap conducts  [1/m]
    K  = diag([g; 0] + [0; g]) - diag(g, 1) - diag(g, -1);   % heat flow between each point and its neighbours
    K(1,1) = K(1,1) + h0;                       % first point: plus heat lost through its surface
    K(n,n) = K(n,n) + h1;                       % last point: same
end
[Q, l] = eig(K./sqrt(w*w.'));                   % solve for the shapes (symmetric form, so eig is accurate)
lam = max(diag(l), 0);                          % fading rates (rounding can give tiny negatives; clip to 0)
V   = Q./sqrt(w);                               % undo the symmetric scaling to get the actual shapes
end


% =========================================================================
% GRADED: grid points that are close together only where they need to be
% =========================================================================
function p = graded(lo, hi, a, b, h_fine, h_coarse, growth)
% Makes the positions of grid points along one direction, from lo to hi:
% evenly spaced h_fine apart inside the fine zone a to b, then spreading out
% on both sides, each gap growth times the one before, until the gaps reach
% h_coarse. Hand the result to abl_forvibes.init as P.grid.x, .y or (negated) .z.
%
% Inputs:   lo, hi     ends of the line  [m]
%           a, b       the fine zone, lo <= a <= b <= hi  [m]
%           h_fine     spacing inside the fine zone  [m]
%           h_coarse   largest spacing far from it  [m]
%           growth     how much each gap may grow on the one before (1.1 to 1.2
%                      keeps the accuracy; sudden jumps in spacing lose it)  [-]
% Outputs:  p          the positions, lo to hi, increasing  [m]
%
% Example:  abl_forvibes.graded(0, 39e-3, 0, 6e-3, 0.25e-3, 2e-3, 1.15)
%           0.25 mm spacing in the top 6 mm, growing to 2 mm by the bottom
nf  = max(1, round((b - a)/h_fine));
mid = linspace(a, b, nf + 1);                   % the fine zone, evenly spaced
% Walk outward from each end of the fine zone, letting the gap grow each step.
% The left side is walked as distances from a, then turned back into positions.
side = {[], []};                                % points left of a, right of b
dist = [a - lo, hi - b];                        % how far each side reaches  [m]
for s = 1:2
    q = 0;  d = (b - a)/nf;                     % distance walked so far, current gap  [m]
    while q(end) < dist(s) - 1e-12
        d = min(d*growth, h_coarse);
        q(end+1) = min(q(end) + d, dist(s));    %#ok<AGROW>
    end
    if numel(q) > 2 && q(end) - q(end-1) < 0.5*(q(end-1) - q(end-2))
        q(end-1) = [];                          %#ok<AGROW> last gap too short: merge it into the one before
    end
    side{s} = q(2:end);
end
p = [a - flip(side{1}), mid, b + side{2}];
end


% =========================================================================
% PROJ: convert a 3-D heat pattern into amounts of each mode
% =========================================================================
function G = proj(F, Vx, Vy, Vz)
% Takes a value at every (x, y, z) grid point and returns the matching amount
% for every (x, y, z) mode. Because the modes are products of 1-D shapes, this
% can be done one direction at a time, which is far faster than all at once:
% first along x, then y, then z.
[nx, ny, nz] = size(F);
G = reshape(Vx.'*reshape(F, nx, []), nx, ny, nz);   % along x
G = pagemtimes(G, Vy);                              % along y
G = reshape(reshape(G, [], nz)*Vz, nx, ny, nz);     % along z
end


% =========================================================================
% PLANE: turn mode amounts back into temperatures on the antenna plane
% =========================================================================
function T = plane(S)
% The opposite of abl_forvibes.proj, but only on the plane y = 0 (the slice through
% the antenna), since that is all the controller looks at.
%   $T(x, 0, z) = T_b + \sum_m c_m \, \phi_m(x, 0, z)$
% Output: T  temperature, one row per x and one column per z  [degC]
on_plane = reshape(sum(S.c.*S.vy, 2), numel(S.x), []);     % mode amounts with the y-shape evaluated at y = 0
T = S.P.T_b + S.Vx*on_plane*S.Vz.';                          % add up the x and z shapes  [degC]
end

end
end
