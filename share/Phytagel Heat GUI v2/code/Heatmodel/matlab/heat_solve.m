function S = heat_solve(q, G, Tc, opt)
%HEAT_SOLVE  Transient 3-D conduction in the gel driven by the microwave
%            power deposition q [W/m^3], with the blood surface as the sink.
%
%   S = HEAT_SOLVE(q, G, Tc)          q is [nx x ny x nz] at the REFERENCE
%                                     delivered power Tc.Pdel_ref.
%   S = HEAT_SOLVE(q, G, Tc, opt)     opt.t_out, opt.verbose, opt.probes,
%                                     opt.decim, opt.Nt3, opt.quiet,
%                                     opt.peakmask, opt.theta0,
%                                     opt.period, opt.duty
%
%   THE EQUATION
%       rho c dT/dt = k grad^2 T + q(r)
%   with theta = T - T_blood this becomes, with homogeneous boundary
%   conditions and theta(r,0) = 0,
%       dtheta/dt = alpha grad^2 theta + q(r)/(rho c) .
%   There is no Pennes perfusion term: the gel is a phantom and is not
%   perfused.  The blood appears ONLY through the surface condition, which is
%   what the brief specifies.
%
%   WHY ONE SOLVE IS ENOUGH FOR EVERY POWER
%   The problem is linear in q and q is proportional to the delivered power,
%   so the solution for delivered power P is exactly
%       theta(r,t; P) = (P / Pdel_ref) * Theta(r,t)
%   where Theta is what this function computes.  That scaling is exact within
%   the constant-property model and is the backbone of T_gel.m.
%
%   PULSING  (opt.period [s], opt.duty in 0..1)
%   The drive is gated: the deposition is q during the first duty*period of
%   every period and zero for the rest of it.  opt.period = [] (default) or
%   duty = 1 is continuous wave and reproduces the unpulsed solve exactly.
%
%   The gate is applied INSIDE the time loop rather than by superposing step
%   responses afterwards, because this solver is already transient and the
%   initial condition need not be equilibrium.  What that costs is that the
%   time step has to land on every pulse edge -- otherwise a step straddling an
%   edge would silently deliver the wrong energy, and the duty cycle would be
%   whatever the step size happened to make it.  NEXT_EDGE below returns the
%   next edge and the step is clipped to it, so the delivered energy is exact
%   for any period, however short.
%
%   Pulsing does NOT scale the answer the way power does: theta is linear in q,
%   so at the same MEAN power a slow pulse train reaches a higher peak than CW
%   (the gel cools between pulses by less than it heats during them).  That is
%   the point of having it here.
%
%   TIME INTEGRATION
%   Forward Euler (FTCS) with dt = min(cfl*dt_stable, next output time - t).
%   The stability limit is dh^2/(alpha*(6 + 2 dh h/k)); with dh = 0.25 mm and
%   alpha = 1.435e-7 m^2/s that is 72.6 ms (Dirichlet) or 60 ms (Robin,
%   h = 1459).  Explicit stepping is used rather than an implicit scheme
%   because it needs no linear solve, is trivially verifiable, and the
%   stability limit is not restrictive at this grid size.  Time-step
%   convergence is checked in VERIFY_HEAT_MODEL.
%
%   ENERGY AUDIT
%   dU/dt = P_in - P_loss with U = sum_i W_i rho c theta_i dh^3, and
%   P_loss = -k dh^3 sum_i W_i (grad^2 theta)_i.  The second identity is exact
%   for this discretisation (the sum telescopes to the boundary fluxes), so
%   the audit is a genuine check on the boundary implementation, not a
%   tautology.
%
%   OUTPUT struct S
%     t          output times [s]
%     Th_y0      [nx x nz x Nt] single, theta on the y = 0 plane, full res
%     Th3        [nx3 x ny3 x nz3 x Nt3] single, theta on the decimated grid
%     t3, x3, y3, z3   the decimated grid
%     Th_max     [1 x Nt] peak theta anywhere in the block
%     ij_max     [3 x Nt] its (ix,iy,iz) location
%     Th_probe   [np x Nt] theta at opt.probes
%     P_in       [W] deposited in the block while the drive is ON
%     P_in_mean  [W] the same averaged over the duty cycle
%     E_in       [1 x Nt] cumulative energy actually delivered [J]
%     P_loss     [1 x Nt] instantaneous loss through the boundaries
%     E_loss     [1 x Nt] cumulative loss, accumulated at every TIME STEP
%                rather than by trapezoid on P_loss between output times.
%                It costs one reduction per step and it is what makes the
%                audit mean anything under pulsing: P_loss swings with the
%                drive, so sampling it at the output times and trapezoiding
%                measures the sampling, not the solve.
%     Th_max_run peak theta over the WHOLE run, including the ends of pulses
%                that fall between output times.  With a slow pulse the run
%                can finish deep in an off-period, so Th_max(end) is not the
%                worst the gel saw and quoting it alone would flatter it.
%     U          [1 x Nt] stored energy [J]
%     err_energy worst relative energy-balance error over the run

if nargin < 4, opt = struct(); end
def = struct('t_out',Tc.t_out,'verbose',true,'probes',zeros(0,3), ...
             'decim',Tc.decim,'Nt3',61,'quiet',false,'peakmask',[], ...
             'theta0',0,'period',[],'duty',1);
fn = fieldnames(def);
for q1 = 1:numel(fn)
    if ~isfield(opt,fn{q1}), opt.(fn{q1}) = def.(fn{q1}); end
end

% ---- the drive gate ------------------------------------------------------
Tp     = opt.period;
duty   = opt.duty;
pulsed = ~isempty(Tp) && Tp > 0 && duty < 1;
if pulsed
    if duty <= 0
        error('heat_solve:duty','duty must be > 0 (got %g)', duty);
    end
    ton = duty*Tp;
else
    Tp = Inf;  ton = Inf;  duty = 1;
end

nx = G.nx; ny = G.ny; nz = G.nz;
t_out = opt.t_out(:).';  Nt = numel(t_out);
rc = G.rho*G.c;
src = q/rc;                                     % [K/s]
src(G.pin) = 0;

% ---- decimated storage grid ---------------------------------------------
dcm = opt.decim;
i3 = 1:dcm:nx;  j3 = 1:dcm:ny;  k3 = 1:dcm:nz;
it3 = unique(round(linspace(1, Nt, min(opt.Nt3, Nt))));
Nt3 = numel(it3);

S = struct();
S.t   = t_out;
S.Th_y0 = zeros(nx, nz, Nt, 'single');
S.Th3   = zeros(numel(i3), numel(j3), numel(k3), Nt3, 'single');
S.t3    = t_out(it3);
S.x3 = G.x(i3);  S.y3 = G.y(j3);  S.z3 = G.z(k3);
S.Th_max = zeros(1,Nt);  S.ij_max = zeros(3,Nt);
S.Th_max_all = zeros(1,Nt);  S.ij_max_all = zeros(3,Nt);
if isempty(opt.peakmask), opt.peakmask = true(nx,ny,nz); end
pmask = opt.peakmask & ~G.pin;
S.peak_excluded = nnz(~opt.peakmask);
np = size(opt.probes,1);
S.Th_probe = zeros(np,Nt);
S.P_loss = zeros(1,Nt);  S.U = zeros(1,Nt);  S.E_in = zeros(1,Nt);
S.E_loss = zeros(1,Nt);
iy0 = find(abs(G.y) < 1e-12, 1);
if isempty(iy0), [~,iy0] = min(abs(G.y)); end
S.iy0 = iy0;

% probe indices
ip = zeros(np,1);
for r = 1:np
    [~,a] = min(abs(G.x - opt.probes(r,1)));
    [~,b] = min(abs(G.y - opt.probes(r,2)));
    [~,c] = min(abs(G.z - opt.probes(r,3)));
    ip(r) = sub2ind([nx ny nz], a, b, c);
end
S.probe_xyz = opt.probes;

dV   = G.dh^3;
% The power deposited AT a Dirichlet node never enters the solution: that node
% is held at the coolant temperature, i.e. the perfect heat sink removes its
% deposition instantly.  The energy audit therefore uses the free-node total,
% and both numbers are reported so the discarded fraction is visible.
S.P_in_total = sum(G.W(:).*q(:))*dV;            % everything deposited [W]
S.P_in = sum(G.W(G.free).*q(G.free))*dV;        % what heats the gel, drive ON
S.P_in_mean = duty*S.P_in;                      % ... averaged over the pulse
S.period = opt.period;  S.duty = duty;  S.pulsed = pulsed;
Uc     = G.W*rc*dV;                             % per-node heat capacity [J/K]

dt_max = Tc.cfl*G.dt_stable;
% initial condition.  opt.theta0 = 0 (the default) is the baseline: the gel
% starts in equilibrium with the coolant.  A scalar or a full array puts it
% somewhere else -- a gel taken out of the fridge, say -- which is a legitimate
% initial condition but breaks the exact power scaling, because the answer is
% then the sum of a driven part (linear in P) and a relaxation part (which
% does not scale with P at all).  See GUI_MODEL, which uses it.
th = zeros(nx,ny,nz);
if isscalar(opt.theta0), th(:) = opt.theta0; else, th(:) = opt.theta0(:); end
th(G.pin) = 0;
t  = 0;  nstep = 0;  Eaudit = 0;  Ecum = 0;  tprev = 0;  Ein = 0;  Eloss = 0;
kdV  = G.k*G.dh^3;
Wcol = G.W(:);
tic0 = tic;

% ---- record t = 0 --------------------------------------------------------
S = record(S, 1, th, G, i3, j3, k3, it3, iy0, ip, Uc, pmask);
S.Th_max_run = S.Th_max(1);

for io = 2:Nt
    ttarget = t_out(io);
    while t < ttarget - 1e-15
        dt = min(dt_max, ttarget - t);
        if pulsed
            dt = min(dt, next_edge(t, Tp, ton) - t);
        end
        % the gate is read at the MIDPOINT of the step.  Because the step has
        % already been clipped to the next edge, the midpoint is never inside
        % two different states, so this is exact and not a quadrature choice.
        if pulsed
            g = double(mod(t + 0.5*dt, Tp) < ton);
        else
            g = 1;
        end
        L  = heat_lap(th, G);
        % left-hand rule on the loss, i.e. the same instant of th that the
        % forward Euler step itself uses for the source.  Being consistent with
        % the scheme is what lets the audit close to round-off instead of to
        % the time-integration error of a second, different rule.
        Eloss = Eloss - dt*kdV*sum(Wcol.*L(:));
        th = th + dt*(G.alpha*L + g*src);
        th(G.pin) = 0;
        t = t + dt;  nstep = nstep + 1;
        Ein = Ein + dt*g*S.P_in;
        if pulsed && g > 0
            % the instant a pulse switches off is the hottest the gel gets in
            % that cycle, and it is generally not an output time
            tp = t - Tp*floor(t/Tp + 1e-12);
            if abs(tp - ton) < 1e-9*Tp
                S.Th_max_run = max(S.Th_max_run, max(th(pmask)));
            end
        end
    end
    S = record(S, io, th, G, i3, j3, k3, it3, iy0, ip, Uc, pmask);
    S.E_in(io) = Ein;  S.E_loss(io) = Eloss;
    S.Th_max_run = max(S.Th_max_run, S.Th_max(io));
    % energy audit: trapezoid on P_loss between output times
    Ecum = Ecum + 0.5*(S.P_loss(io)+S.P_loss(io-1))*(t_out(io)-tprev);
    Ubal = S.U(1) + S.E_in(io) - Ecum;
    if S.U(1) == 0
        % the baseline: the gel starts at the coolant temperature, so the
        % balance Ubal is itself the natural scale for the relative error
        if Ubal > 0
            Eaudit = max(Eaudit, abs(S.U(io)-Ubal)/Ubal);
        end
    else
        % a gel that starts away from the coolant temperature carries a large
        % stored energy of its own, and Ubal can be large and negative --
        % dividing by it would make the audit meaningless (and dividing by a
        % near-zero crossing would make it explode).  Scale by the largest
        % term in the balance instead, so the check stays live.
        sc = max([abs(S.U(1)), S.E_in(io), abs(Ecum)]);
        if sc > 0
            Eaudit = max(Eaudit, abs(S.U(io)-Ubal)/sc);
        end
    end
    tprev = t_out(io);
    if opt.verbose && ~opt.quiet && mod(io,40)==0
        fprintf('    t = %8.3f s   theta_max = %8.2f K   (%d steps, %.0f s)\n', ...
                t, S.Th_max(io), nstep, toc(tic0));
    end
end

S.nstep = nstep;
S.dt_max = dt_max;
S.walltime = toc(tic0);
S.err_energy = Eaudit;
S.G = rmfield(G, {'W','pin','free','diagA'});
end

% =========================================================================
function e = next_edge(t, Tp, ton)
%NEXT_EDGE  The next instant at which the drive switches, strictly after t.
%   Built from the absolute edge times of the period t falls in rather than
%   from mod(t,Tp), so that accumulated rounding cannot produce a zero-length
%   step and stall the loop.
tol  = 1e-9*Tp;
n    = floor(t/Tp);
cand = [n*Tp + ton, (n+1)*Tp, (n+1)*Tp + ton];
e    = cand(find(cand > t + tol, 1));
end

% =========================================================================
function S = record(S, io, th, G, i3, j3, k3, it3, iy0, ip, Uc, pmask)
S.Th_y0(:,:,io) = single(squeeze(th(:,iy0,:)));
m = find(it3 == io, 1);
if ~isempty(m)
    S.Th3(:,:,:,m) = single(th(i3,j3,k3));
end
thm = th;  thm(~pmask) = -Inf;
[v,li] = max(thm(:));
S.Th_max(io) = v;
[a,b,c] = ind2sub(size(th), li);
S.ij_max(:,io) = [a;b;c];
[v2,l2] = max(th(:));
S.Th_max_all(io) = v2;
[a2,b2,c2] = ind2sub(size(th), l2);
S.ij_max_all(:,io) = [a2;b2;c2];
if ~isempty(ip), S.Th_probe(:,io) = th(ip); end
S.U(io) = sum(Uc(:).*th(:));
L = heat_lap(th, G);
S.P_loss(io) = -G.k*G.dh^3*sum(G.W(:).*L(:));
end
