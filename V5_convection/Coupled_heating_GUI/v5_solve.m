function R = v5_solve(P, progress)
%V5_SOLVE  The four stored responses of the V5 gel model, on the antenna plane.
%
%   R = v5_solve(P)               P from v5_defaults (only the 'solve' fields matter)
%   R = v5_solve(P, @(frac,msg) ...)   progress callback
%
%   The heat equation and its boundary conditions are linear (Convection
%   model.docx s3), so every instant setting is a combination of four
%   responses, each started from zero:
%     D  dielectric: source sigma E_rms^2 for 1 W delivered, switched on at t = 0
%     U  diffusion:  surface water = 1 K over the catheter footprint 0 <= x <= L_c
%     Wr walls:      room air = 1 K outside the tank walls and bottom
%     Z  start:      gel at 1 K at t = 0, everything outside at 0
%   T = T_b + P_del D(*) + alpha Theta(*)U + (T_room - T_b) Wr + (T_0 - T_b) Z
%   where (*) marks the switch-off and the Duhamel sum (v5_eval).
%
%   Method.  Finite volumes on a uniform box grid (nodes on the faces),
%   rho c W dT/dt = -k S T + F.  W and S separate into 1-D matrices, so the
%   generalised eigenvectors of the three 1-D problems diagonalise the whole
%   operator and each mode is integrated EXACTLY in time: no time step,
%   only the grid spacing dh sets the accuracy.  3-D uses the symmetry
%   plane y = 0 (the antenna plane, adiabatic); 2-D is the (x,z) slice with
%   the source taken as uniform in y (overestimates: no sideways spreading).
%   Faces: top = Robin h_s to the water, the others = Robin h_w to room air.

if nargin < 2 || isempty(progress), progress = @(varargin) []; end
t_all = tic;
kap = P.k/(P.rho*P.c);                             % diffusivity (m^2/s)
L_c = P.L_c;  if isnan(L_c), L_c = P.L; end
pin = 1e8;                                          % h/k (1/m) standing in for h = Inf

% ------------------------------------------------------------ GRID ------
nx = max(3, round(P.gel_L/P.dh) + 1);  x = P.x_start + linspace(0, P.gel_L, nx);
nz = max(3, round(P.gel_H/P.dh) + 1);  z = -linspace(0, P.gel_H, nz);      % z(1) = surface
if P.dim == 3
    ny = max(3, round(P.gel_W/2/P.dh) + 1);  y = linspace(0, P.gel_W/2, ny);  % y(1) = antenna plane
else
    ny = 1;  y = 0;
end
hw = P.h_w/P.k;  hs = P.h_s/P.k;  if isinf(hs), hs = pin; end
[Vx, lx, wx] = modes1d(x, hw, hw);
[Vz, lz, wz] = modes1d(-z, hs, hw);                 % node 1 = top face
if P.dim == 3
    [Vy, ly, wy] = modes1d(y, 0, hw);               % y = 0 symmetry plane
else
    Vy = 1;  ly = 0;  wy = 1;
end
N = nx*ny*nz;

% ------------------------------------------------------- SOURCE ---------
progress(0.02, sprintf('E-field: %d x-stations ...', nx));
t1 = tic;
if P.dim == 3
    rmax = hypot(P.gel_W/2, P.d + P.gel_H);
    rho_t = exp(linspace(log(P.d), log(rmax*1.001), 300));
else
    rho_t = P.d - z;                                % exact nodes in 2-D
end
[Et, einfo] = v5_erms(P, x, rho_t);                 % numel(rho_t) x nx, at 1 W delivered
sigma = einfo.sigma;
if P.dim == 3
    F = griddedInterpolant({log(rho_t), x}, log(Et), 'linear');
    [XX, YY, ZZ] = ndgrid(x, y, z);
    Q = sigma*exp(F(log(sqrt(YY.^2 + (P.d - ZZ).^2)), XX)).^2;    % nx x ny x nz (W/m^3 per W)
    clear XX YY ZZ
else
    Q = reshape(sigma*(Et.').^2, nx, 1, nz);
end
t_E = toc(t1);

% --------------------------------------------- PROJECT THE FORCINGS -----
progress(0.35, sprintf('%d nodes: projecting on the modes ...', N));
t2 = tic;
Wq  = Q.*reshape(wx, [], 1, 1).*reshape(wy, 1, [], 1).*reshape(wz, 1, 1, []);
gD  = modal(Wq, Vx, Vy, Vz);
P_gel = sum(Wq(:))*(1 + (P.dim == 3));             % W absorbed in the block per W delivered (x2: both halves)
if P.dim == 2, P_gel = NaN; end                     % per metre of y in 2-D: not a power

% surface water over the footprint: fraction of each top control volume inside [0, L_c]
lo = max(x - [0 diff(x)]/2, x(1));  hi = min(x + [diff(x) 0]/2, x(end));
fx = max(0, min(hi, L_c) - max(lo, 0))./(hi - lo);
gU = kron3(Vx.'*(wx(:).*fx(:)), Vy.'*wy(:), Vz(1,:).'*P.k*hs);

% walls: x faces, far y face, bottom
one = @(V, w) V.'*w(:);
gW = kron3(Vx(1,:).'*P.k*hw + Vx(end,:).'*P.k*hw, one(Vy, wy), one(Vz, wz)) ...
   + kron3(one(Vx, wx), one(Vy, wy), Vz(end,:).'*P.k*hw);
if P.dim == 3
    gW = gW + kron3(one(Vx, wx), Vy(end,:).'*P.k*hw, one(Vz, wz));
end
gZ = kron3(one(Vx, wx), one(Vy, wy), one(Vz, wz));   % initial 1 K

% exact modal time response on the antenna plane y = 0
Lam = kap*(reshape(lx, [], 1, 1) + reshape(ly, 1, [], 1) + reshape(lz, 1, 1, []));
vy1 = reshape(Vy(1,:), 1, [], 1);
gD = gD.*vy1/(P.rho*P.c);  gU = gU.*vy1/(P.rho*P.c);  gW = gW.*vy1/(P.rho*P.c);  gZ = gZ.*vy1;
t_M = toc(t2);

t = 0:P.dt_frame:P.t_end;  nt = numel(t);
[D, U, Wr, Z] = deal(zeros(nz, nx, nt, 'single'));
t3 = tic;
small = Lam < 1e-12;
for n = 1:nt
    if mod(n, 10) == 1
        progress(0.4 + 0.58*n/nt, sprintf('frames: t = %.0f of %.0f s ...', t(n), P.t_end));
    end
    Et_ = exp(-Lam*t(n));
    phi = -expm1(-Lam*t(n))./Lam;  phi(small) = t(n);
    D(:,:,n)  = plane(gD.*phi, Vx, Vz);
    U(:,:,n)  = plane(gU.*phi, Vx, Vz);
    Wr(:,:,n) = plane(gW.*phi, Vx, Vz);
    Z(:,:,n)  = plane(gZ.*Et_, Vx, Vz);
end
t_F = toc(t3);

R = struct();
R.P = P;  R.x = x;  R.z = z;  R.t = t;  R.dim = P.dim;  R.N = N;  R.nx = nx;  R.ny = ny;  R.nz = nz;
R.D = D;  R.U = U;  R.Wr = Wr;  R.Z = Z;
R.Q = reshape(Q(:,1,:), nx, nz).';                  % antenna plane, W/m^3 per W delivered
R.E = sqrt(R.Q/sigma);                              % V/m at 1 W delivered
R.info = einfo;  R.P_gel = P_gel;  R.L_c = L_c;
R.timing = struct('E', t_E, 'modes', t_M, 'frames', t_F, 'total', toc(t_all), ...
                  'nseg', round(P.L/min(0.05e-3, P.L/200)) + round(P.l_f/min(0.05e-3, P.L/200)), ...
                  'nrho', numel(rho_t));
R.key = v5_key(P);
progress(1, 'done');
end

% =========================================================================
function [V, lam, w] = modes1d(s, hk_lo, hk_hi)
% finite-volume 1-D operator on nodes s (uniform), Robin h/k at both ends:
% S V = W V diag(lam),  V' W V = I
n = numel(s);  ds = s(2) - s(1);
w = ds*ones(n, 1);  w([1 n]) = ds/2;
e = ones(n, 1)/ds;
S = spdiags([-e 2*e -e], -1:1, n, n);
S(1,1) = 1/ds + hk_lo;  S(n,n) = 1/ds + hk_hi;
Wm = 1./sqrt(w);
A = full(Wm.*S.*Wm.');  A = (A + A.')/2;
[Qm, Lm] = eig(A);
lam = max(diag(Lm), 0);
V = Wm.*Qm;
end

function G = modal(F, Vx, Vy, Vz)
% G = (Vx (x) Vy (x) Vz)' F, for F an nx x ny x nz array
[nx, ny, nz] = size(F);
G = reshape(Vx.'*reshape(F, nx, []), nx, ny, nz);
G = permute(reshape(Vy.'*reshape(permute(G, [2 1 3]), ny, []), ny, nx, nz), [2 1 3]);
G = reshape(reshape(G, [], nz)*Vz, nx, ny, nz);
end

function G = kron3(a, b, c)
G = reshape(a, [], 1, 1).*reshape(b, 1, [], 1).*reshape(c, 1, 1, []);
end

function T = plane(C, Vx, Vz)
% nodal values on y = 0 (the y-eigenvector row is already folded into C)
Y = reshape(sum(C, 2), size(C, 1), size(C, 3));
T = single((Vx*Y*Vz.').');                         % nz x nx
end
