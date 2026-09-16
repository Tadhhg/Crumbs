function R = gel_heating_solve(P, varargin)
%GEL_HEATING_SOLVE  One-off 3-D solve of the gel heating problem for the
%                   constants in P, stored as unit responses so that any
%                   power / water temperature / timing can be composed
%                   instantly (gel_heating_eval, the GUI, the optimiser).
%
%   R = GEL_HEATING_SOLVE()            model_defaults (the photographed setup)
%   R = GEL_HEATING_SOLVE(P)           P from model_defaults, fields edited
%   R = GEL_HEATING_SOLVE(P, 'force', true)   re-solve even if cached
%
%   Every distinct P is solved once and cached in solutions/gel_heating_<key>.mat
%   (key = MD5 of P, model_derived).  ~10 s at dh = 0.5 mm.
%
%   The problem is LINEAR.  With the gel initially in equilibrium with water
%   at T_bath0 (the tank walls to the room), water at T_w from t = -t_pre and
%   microwave power P_set from t = 0:
%
%       T(x,z,t) = T_ss + (T_w - T_bath0) U(t_pre + t) + (P_set/P_ref) Theta(t)
%
%   T_ss   steady state without power (uniform T_room when T_bath0 = T_room)
%   U      response to a unit step of the water temperature (no power)
%   Theta  response to the microwave source at P_ref, water held at T_bath0
%   Both are solved in 3-D and stored on the antenna plane y = 0.
%
%   Chain: E_field_2d (homogeneous lossy gel, the medium starting at the
%   insulation radius b_ins) -> q = sig_eff E_rms^2, axisymmetric about the
%   wire -> FTCS conduction with Robin faces (top: water; sides/bottom: tank
%   wall to room), or a pinned top face when h_top = Inf -> catheter lumped
%   balance (cylinder in cross-flow).  See ../Gel_model_V3/gel_heating_notes.md.
%
%   Fields of R used downstream: x, z, t, T_ss, U, Theta (single, [nx nz nfr]),
%   T_bath0, T_room, P_ref, eta, L_el, h_ax, cath, P_gel, P_water, P_film,
%   P (the constants), D (model_derived), file, created.

if nargin < 1 || isempty(P), P = model_defaults(); end
o = struct('force', false, 'dir', fullfile(fileparts(mfilename('fullpath')), 'solutions'));
for q1 = 1:2:numel(varargin), o.(varargin{q1}) = varargin{q1+1}; end

D = model_derived(P);
file = fullfile(o.dir, sprintf('gel_heating_%s.mat', D.key));
if ~o.force && isfile(file)
    S = load(file);  R = S.R;  return
end
if ~isfolder(o.dir), mkdir(o.dir); end

rho = P.rho;  c = P.c;  k = P.k;  alpha = D.alpha;
x = D.x;  y = D.y;  z = D.z;  nx = D.nx;  ny = D.ny;  nz = D.nz;  dh = P.dh;  iy0 = D.iy0;
[X3, Y3, Z3] = ndgrid(x, y, z);
h_ax = P.h_axis;  L_el = P.L;

% ------------------------------------------ SOURCE FROM THE E-FIELD MODEL --
fprintf('gel_heating_solve: E-field source ... ');
ant = {'f', P.f, 'epsr', D.epsr_gel_c, 'L', L_el, 'eta', P.eta, 'h', h_ax, ...
       'a', P.b_ins, 'L_sheath', P.L_sheath, 'ke', D.ke};
rho_lo  = 0.9*h_ax;                                       % nearest gel node is at rho = h_ax
rho_hi  = 1.05*hypot(P.W_gel/2, P.H_gel + h_ax);
rho_tab = logspace(log10(rho_lo), log10(rho_hi), 200);
[Xt, Rt] = ndgrid(x, rho_tab);
[~, S] = E_field_2d(Xt, h_ax - Rt, 0, P.P_ref, ant{:});
R3 = hypot(Y3, Z3 - h_ax);
q  = interp2(log(rho_tab), x, S.q, log(max(R3, rho_lo)), X3, 'linear');
q(isnan(q)) = 0;
P_gel = sum(q(:))*dh^3;

% catheter: cross-flow film coefficient and the power dissipated in that film
Re_c = D.U_w*D.D_c/P.nu_w;
Nu_c = 0.3 + 0.62*sqrt(Re_c)*P.Pr_w^(1/3)/(1 + (0.4/P.Pr_w)^(2/3))^0.25 ...
           *(1 + (Re_c/282000)^(5/8))^(4/5);              % Churchill-Bernstein
h_c  = Nu_c*P.k_w/D.D_c;
d_t  = P.k_w/h_c;
rf = linspace(P.b_ins, P.b_ins + d_t, 40);  xf = -L_el : 0.25e-3 : L_el;
[Xf, Rf] = ndgrid(xf, rf);
[~, Sf] = E_field_2d(Xf, h_ax - Rf, 0, P.P_ref, ant{:});
P_film  = trapz(xf, trapz(rf, Sf.q.*Rf*2*pi, 2));
P_water = S.P_del - P_gel;
fprintf('%.2f W of %.1f W in the gel\n', P_gel, S.P_del);

% ------------------------------------------ 3-D STEADY STATE (no source) --
g_w = D.g_w;
if abs(P.T_bath0 - P.T_room) < 1e-9
    T_ss = P.T_room*ones(nx, ny, nz);                      % exact: every T_inf equal
else
    g_top_ss = D.g_top;  if D.pinned, g_top_ss = 2*dh*1e8/k; end   % Robin -> Dirichlet
    [Lx, bx] = lap1d(nx, g_w, g_w, P.T_room, P.T_room);
    [Ly, by] = lap1d(ny, g_w, g_w, P.T_room, P.T_room);
    [Lz, bz] = lap1d(nz, g_w, g_top_ss, P.T_room, P.T_bath0);
    Ix = speye(nx);  Iy = speye(ny);  Iz = speye(nz);
    L3 = kron(kron(Iz, Iy), Lx) + kron(kron(Iz, Ly), Ix) + kron(kron(Lz, Iy), Ix);
    B3 = bx(:) + reshape(by, 1, []) + reshape(bz, 1, 1, []);
    wx = ones(nx, 1);  wx([1 end]) = 0.5;                 % half-weights on the faces make
    wy = ones(ny, 1);  wy([1 end]) = 0.5;                 % the ghost-node operator symmetric
    wz = ones(nz, 1);  wz([1 end]) = 0.5;
    W3 = kron(kron(wz, wy), wx);
    A  = -spdiags(W3, 0, numel(W3), numel(W3))*L3;
    Lc = ichol(A);
    [T_ss, flag] = pcg(A, W3.*B3(:), 1e-9, 4000, Lc, Lc');
    if flag ~= 0
        warning('gel_heating_solve:pcg', 'steady state: pcg flag %d, using a direct solve', flag);
        T_ss = A \ (W3.*B3(:));
    end
    T_ss = reshape(T_ss, nx, ny, nz);
end

% ---------------------- TRANSIENT UNIT RESPONSES, homogeneous walls --
dt = D.dt;  npf = D.npf;  nfr = D.nfr;
r   = alpha*dt/dh^2;
src = q/(rho*c)*dt;
t_fr = (0:nfr-1)*P.dt_frame;
Th = zeros(nx, ny, nz);  Uh = zeros(nx, ny, nz);
TgT = zeros(nx+2, ny+2, nz+2);  TgU = TgT;
Theta = zeros(nx, nz, nfr, 'single');  U = Theta;  Theta_max3d = zeros(1, nfr);
fprintf('gel_heating_solve: %d steps of 3-D conduction on %dx%dx%d ... ', D.nsteps, nx, ny, nz);
tic
for f = 2:nfr
    for s = 1:npf
        [Th, TgT] = step(Th, TgT, g_w, D.g_top, r, src, 0, D.pinned);   % power on, water at T_bath0
        [Uh, TgU] = step(Uh, TgU, g_w, D.g_top, r, 0,   1, D.pinned);   % water step of +1 K
    end
    Theta(:, :, f) = single(squeeze(Th(:, iy0, :)));
    U(:, :, f)     = single(squeeze(Uh(:, iy0, :)));
    Theta_max3d(f) = max(Th(:));
end
t_solve = toc;
fprintf('%.0f s\n', t_solve);

A_c   = pi*D.D_c*2*L_el;                                  % wetted area: element + braid arm
G_gap = P.k_w*D.D_c*L_el/D.gap;
R = struct('x', x, 'z', z, 'dh', dh, 't', t_fr, ...
           'T_ss', squeeze(T_ss(:, iy0, :)), 'U', U, 'Theta', Theta, ...
           'Theta_max3d', Theta_max3d, 'P_ref', P.P_ref, 'eta', P.eta, ...
           'T_bath0', P.T_bath0, 'T_room', P.T_room, 'h_top', D.h_top, 'h_wall', D.h_wall, ...
           'P_gel', P_gel, 'P_water', P_water, 'P_film', P_film, ...
           'T_ss_range', [min(T_ss(:)) max(T_ss(:))], 'h_ax', h_ax, 'L_el', L_el, ...
           'ke', D.ke, 'P', P, 'D', D, 'file', file, 't_solve', t_solve, 'created', datetime('now'));
R.cath = struct('D_c', D.D_c, 'Re', Re_c, 'Nu', Nu_c, 'h_c', h_c, 'A_c', A_c, ...
                'G_gap', G_gap, 'd_t', d_t, 'P_film', P_film, 'L_c', 2*L_el);
save(file, 'R');
end

% =========================================================================
function [T, Tg] = step(T, Tg, g_w, g_top, r, src, top_inf, pinned)
% one FTCS step with Robin ghost faces; walls see 0, the top sees top_inf
% (or is held at top_inf when pinned)
Tg(2:end-1, 2:end-1, 2:end-1) = T;
Tg(1,   2:end-1, 2:end-1) = T(2,     :, :) - g_w*T(1,   :, :);
Tg(end, 2:end-1, 2:end-1) = T(end-1, :, :) - g_w*T(end, :, :);
Tg(2:end-1, 1,   2:end-1) = T(:, 2,     :) - g_w*T(:, 1,   :);
Tg(2:end-1, end, 2:end-1) = T(:, end-1, :) - g_w*T(:, end, :);
Tg(2:end-1, 2:end-1, 1  ) = T(:, :, 2    ) - g_w*T(:, :, 1  );
if pinned
    Tg(2:end-1, 2:end-1, end) = T(:, :, end-1);                     % irrelevant: face is overwritten
else
    Tg(2:end-1, 2:end-1, end) = T(:, :, end-1) - g_top*(T(:, :, end) - top_inf);
end
T = T + r*( Tg(1:end-2,2:end-1,2:end-1) + Tg(3:end,2:end-1,2:end-1) ...
          + Tg(2:end-1,1:end-2,2:end-1) + Tg(2:end-1,3:end,2:end-1) ...
          + Tg(2:end-1,2:end-1,1:end-2) + Tg(2:end-1,2:end-1,3:end) - 6*T ) + src;
if pinned, T(:, :, end) = top_inf; end
end

function [L, b] = lap1d(n, g_lo, g_hi, T_lo, T_hi)
% 1-D second difference with Robin ends folded in via the ghost node
e = ones(n, 1);
L = spdiags([e -2*e e], -1:1, n, n);
L(1, 2) = 2;    L(1, 1) = -(2 + g_lo);
L(n, n-1) = 2;  L(n, n) = -(2 + g_hi);
b = zeros(n, 1);  b(1) = g_lo*T_lo;  b(n) = g_hi*T_hi;
end
