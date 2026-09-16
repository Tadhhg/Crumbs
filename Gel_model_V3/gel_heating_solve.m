function R = gel_heating_solve(varargin)
%GEL_HEATING_SOLVE  One-off 3-D solve of the gel heating problem, stored as
%                   unit responses so that any power / water temperature /
%                   timing can be composed instantly (the GUI and the
%                   optimiser both build on this).
%
%   R = GEL_HEATING_SOLVE()          loads gel_heating_solution_gap0.5mm.mat if
%                                    it exists, otherwise solves (~1 min) and saves
%   R = GEL_HEATING_SOLVE('gap',2e-3)     catheter 2 mm above the gel (own file)
%   R = GEL_HEATING_SOLVE('force',true)   re-solve
%   R = GEL_HEATING_SOLVE('L',12e-3,'L_sheath',4e-3,'ke',172-26j,'gap',0.65e-3)
%                                    the antenna of the 50 W photographs (12 mm
%                                    element, short braid arm, insulated-wire
%                                    wavenumber) -- see gel_heating_notes.md s6.
%                                    Non-default antenna options go into the
%                                    file name, so every variant is cached.
%
%   The problem is LINEAR, so with the gel initially in equilibrium with
%   37 degC water (T_ss37) and, from t = -t_pre, water at T_w, and from t = 0
%   microwave power P (generator setting) for t seconds:
%
%       T(x,z,t) = T_ss37 + (T_w - 37) U(t_pre + t) + (P/P_ref) Theta(t)
%
%   U     = response to a unit step of the water temperature (no power)
%   Theta = response to the microwave source at P_ref, water held at 37 degC
%   Both are solved in 3-D (the source is concentrated in y as well as z) and
%   stored on the antenna plane y = 0, which is where every maximum lies.
%
%   Physics, sources and the catheter model: see run_gel_heating.m.
%   Fields returned (plane arrays are [nx x nz], single):
%     x, z, dh, t (frame times), T_ss37, U, Theta, Theta_max3d
%     P_ref, eta, T_room, h_top, h_wall, cath (struct), P_gel, P_water, P_film

p = struct('P_ref',100, 'eta',0.2, 't_end',720, 'dt_frame',2, ...
           'gap',0.5e-3, 'dh',0.5e-3, 'file','', 'force',false, ...
           'L',27e-3, 'L_sheath',27e-3, 'ke',[]);   % antenna, passed to E_field_2d
for q1 = 1:2:numel(varargin), p.(varargin{q1}) = varargin{q1+1}; end
if isempty(p.file)
    p.file = sprintf('gel_heating_solution_gap%.1fmm', p.gap*1e3);
    if p.L ~= 27e-3,        p.file = sprintf('%s_L%.0fmm',  p.file, p.L*1e3);        end
    if p.L_sheath ~= 27e-3, p.file = sprintf('%s_arm%.0fmm', p.file, p.L_sheath*1e3); end
    if ~isempty(p.ke),      p.file = sprintf('%s_ke%.0f%+.0fj', p.file, real(p.ke), imag(p.ke)); end
    p.file = [p.file '.mat'];
end
if ~p.force && isfile(p.file)
    R = load(p.file);  R = R.R;  return
end

% ------------------------------------------ GEL, WATER, WALLS (gel_heat_2d) --
rho = 1000;  c = 4180;  k = 0.60;  alpha = k/(rho*c);
T_water = 37;  T_room = 25;
L_gel = 30e-3;  H_gel = 30e-3;  Q_circ = 5/60000;
nu_w = 6.98e-7;  k_w = 0.628;  Pr_w = 4.61;               % water, 37 degC
U_w   = Q_circ/(L_gel*H_gel);                              % 0.093 m/s
h_top = 0.664*sqrt(U_w*L_gel/nu_w)*Pr_w^(1/3)*k_w/L_gel;   % ~1460 W/(m^2 K)
h_wall = 1/(1/10 + 5e-3/0.19);                             % ~7.9 W/(m^2 K)

% ----------------------------------------------------------- CATHETER --
a     = 0.45e-3;                % [m] conductor radius (E-field model)
D_c   = 2*a;                    % [m] bare element diameter
h_ax  = a + p.gap;              % [m] wire axis above the gel surface
L_el  = p.L;    L_c = 2*L_el;   % [m] element; wetted length incl. braid arm
Re_c  = U_w*D_c/nu_w;
Nu_c  = 0.3 + 0.62*sqrt(Re_c)*Pr_w^(1/3)/(1+(0.4/Pr_w)^(2/3))^0.25 ...
            *(1 + (Re_c/282000)^(5/8))^(4/5);              % Churchill-Bernstein
h_c   = Nu_c*k_w/D_c;
d_t   = k_w/h_c;
A_c   = pi*D_c*L_c;
G_gap = k_w*D_c*L_el/p.gap;

% ----------------------------------------------------------------- GRID --
dh = p.dh;
x = -1.5e-3 : dh : (L_el + 1.5e-3);   y = -15e-3 : dh : 15e-3;   z = -30e-3 : dh : 0;
nx = numel(x);  ny = numel(y);  nz = numel(z);
[X3, Y3, Z3] = ndgrid(x, y, z);
iy0 = find(abs(y) < 1e-9, 1);

% ------------------------------------------ SOURCE FROM THE E-FIELD MODEL --
fprintf('gel_heating_solve: E-field source ... ');
rho_tab = logspace(log10(0.6e-3), log10(70e-3), 200);
[Xt, Rt] = ndgrid(x, rho_tab);
ant = {'h', h_ax, 'eta', p.eta, 'L', p.L, 'L_sheath', p.L_sheath, 'ke', p.ke};
[~, S] = E_field_2d(Xt, h_ax - Rt, 0, p.P_ref, ant{:});
R3 = hypot(Y3, Z3 - h_ax);
q  = interp2(log(rho_tab), x, S.q, log(max(R3, rho_tab(1))), X3, 'linear');
q(isnan(q)) = 0;
P_gel = sum(q(:))*dh^3;
rf = linspace(a, a + d_t, 40);  xf = -L_el : 0.25e-3 : L_el;
[Xf, Rf] = ndgrid(xf, rf);
[~, Sf] = E_field_2d(Xf, h_ax - Rf, 0, p.P_ref, ant{:});
P_film  = trapz(xf, trapz(rf, Sf.q.*Rf*2*pi, 2));
P_water = S.P_del - P_gel;
fprintf('%.1f W of %.1f W in the gel\n', P_gel, S.P_del);

% ------------------------------------------ 3-D STEADY STATE (no source) --
g_top = 2*dh*h_top/k;  g_w = 2*dh*h_wall/k;
[Lx, bx] = lap1d(nx, g_w, g_w, T_room, T_room);
[Ly, by] = lap1d(ny, g_w, g_w, T_room, T_room);
[Lz, bz] = lap1d(nz, g_w, g_top, T_room, T_water);
Ix = speye(nx);  Iy = speye(ny);  Iz = speye(nz);
L3 = kron(kron(Iz, Iy), Lx) + kron(kron(Iz, Ly), Ix) + kron(kron(Lz, Iy), Ix);
B3 = bx(:) + reshape(by, 1, []) + reshape(bz, 1, 1, []);
wx = ones(nx, 1);  wx([1 end]) = 0.5;        % half-weights on the faces make
wy = ones(ny, 1);  wy([1 end]) = 0.5;        % the ghost-node operator symmetric
wz = ones(nz, 1);  wz([1 end]) = 0.5;
W3 = kron(kron(wz, wy), wx);                 % same ordering as L3 (x fastest)
A  = -spdiags(W3, 0, numel(W3), numel(W3))*L3;
Lc = ichol(A);
[T_ss, flag] = pcg(A, W3.*B3(:), 1e-9, 2000, Lc, Lc');
assert(flag == 0, 'steady-state solve did not converge');
T_ss = reshape(T_ss, nx, ny, nz);

% ---------------------- TRANSIENT UNIT RESPONSES, homogeneous walls --
dt_max = 0.9*dh^2/(alpha*(6 + g_top + 2*g_w));
npf = ceil(p.dt_frame/dt_max);  dt = p.dt_frame/npf;       % land on every frame
r  = alpha*dt/dh^2;
src = q/(rho*c)*dt;
nfr = round(p.t_end/p.dt_frame) + 1;
t_fr = (0:nfr-1)*p.dt_frame;
Th = zeros(nx, ny, nz);  Uh = zeros(nx, ny, nz);
TgT = zeros(nx+2, ny+2, nz+2);  TgU = TgT;
Theta = zeros(nx, nz, nfr, 'single');  U = Theta;  Theta_max3d = zeros(1, nfr);
fprintf('gel_heating_solve: %d steps of 3-D conduction ... ', (nfr-1)*npf);
tic
for f = 2:nfr
    for s = 1:npf
        [Th, TgT] = step(Th, TgT, g_w, g_top, r, src, 0);   % power on, water at 37
        [Uh, TgU] = step(Uh, TgU, g_w, g_top, r, 0,   1);   % water step of +1 K
    end
    Theta(:, :, f) = single(squeeze(Th(:, iy0, :)));
    U(:, :, f)     = single(squeeze(Uh(:, iy0, :)));
    Theta_max3d(f) = max(Th(:));
end
fprintf('%.0f s\n', toc);

R = struct('x', x, 'z', z, 'dh', dh, 't', t_fr, ...
           'T_ss37', squeeze(T_ss(:, iy0, :)), 'U', U, 'Theta', Theta, ...
           'Theta_max3d', Theta_max3d, 'P_ref', p.P_ref, 'eta', p.eta, ...
           'T_room', T_room, 'h_top', h_top, 'h_wall', h_wall, ...
           'P_gel', P_gel, 'P_water', P_water, 'P_film', P_film, ...
           'T_ss_range', [min(T_ss(:)) max(T_ss(:))], 'gap', p.gap, 'h_ax', h_ax, ...
           'L_el', L_el, 'L_sheath', p.L_sheath, 'ke', S.ke, 'created', datetime('now'));
R.cath = struct('D_c', D_c, 'Re', Re_c, 'Nu', Nu_c, 'h_c', h_c, 'A_c', A_c, ...
                'G_gap', G_gap, 'd_t', d_t, 'P_film', P_film, 'L_c', L_c);
save(p.file, 'R');
end

% =========================================================================
function [T, Tg] = step(T, Tg, g_w, g_top, r, src, top_inf)
% one FTCS step with Robin ghost faces; walls see 0, the top sees top_inf
Tg(2:end-1, 2:end-1, 2:end-1) = T;
Tg(1,   2:end-1, 2:end-1) = T(2,     :, :) - g_w  *T(1,   :, :);
Tg(end, 2:end-1, 2:end-1) = T(end-1, :, :) - g_w  *T(end, :, :);
Tg(2:end-1, 1,   2:end-1) = T(:, 2,     :) - g_w  *T(:, 1,   :);
Tg(2:end-1, end, 2:end-1) = T(:, end-1, :) - g_w  *T(:, end, :);
Tg(2:end-1, 2:end-1, 1  ) = T(:, :, 2    ) - g_w  *T(:, :, 1  );
Tg(2:end-1, 2:end-1, end) = T(:, :, end-1) - g_top*(T(:, :, end) - top_inf);
T = T + r*( Tg(1:end-2,2:end-1,2:end-1) + Tg(3:end,2:end-1,2:end-1) ...
          + Tg(2:end-1,1:end-2,2:end-1) + Tg(2:end-1,3:end,2:end-1) ...
          + Tg(2:end-1,2:end-1,1:end-2) + Tg(2:end-1,2:end-1,3:end) - 6*T ) + src;
end

function [L, b] = lap1d(n, g_lo, g_hi, T_lo, T_hi)
% 1-D second difference with Robin ends folded in via the ghost node
e = ones(n, 1);
L = spdiags([e -2*e e], -1:1, n, n);
L(1, 2) = 2;    L(1, 1) = -(2 + g_lo);
L(n, n-1) = 2;  L(n, n) = -(2 + g_hi);
b = zeros(n, 1);  b(1) = g_lo*T_lo;  b(n) = g_hi*T_hi;
end
