function D = model_derived(P)
%MODEL_DERIVED  Everything that follows from the constants in P without a
%               solve: wavenumbers, film coefficients, the grid, the time
%               step, the cost estimate and the cache key.  The GUI shows
%               these on the Constants tab; gel_heating_solve uses them.

eps0 = 8.8541878128e-12;  mu0 = 1.25663706212e-6;
w  = 2*pi*P.f;  k0 = w*sqrt(mu0*eps0);
D  = struct();

% ------------------------------------------------------------- GEL, EM --
D.alpha      = P.k/(P.rho*P.c);                          % [m^2/s]
D.epsr_gel_c = P.epsr_gel - 1j*P.epsi_gel;
kc = k0*sqrt(D.epsr_gel_c);  if imag(kc) > 0, kc = conj(kc); end
D.kc         = kc;
D.sig_eff    = w*eps0*P.epsi_gel;                        % [S/m]
D.lambda_gel = 2*pi/real(kc);                            % [m]
D.d_1e_power = 1/(2*(-imag(kc)));                        % [m] power 1/e depth

% -------------------------------------------------------- CURRENT WAVE --
if isempty(P.ke)
    [D.ke, D.king] = king_insulated_ke(P.a, P.b_ins, P.epsr_ins, P.epsr_w, P.f);
    D.ke_source = 'King (insulated wire in water)';
else
    D.ke = P.ke;  D.king = [];  D.ke_source = 'manual';
end
D.lambda_e = 2*pi/real(D.ke);
D.L_over_lambda = P.L/D.lambda_e;

% ------------------------------------------------------------ TOP FACE --
L_plate    = P.x_max - P.x_min;
D.U_w      = P.Q_circ/P.A_flow;                          % [m/s]
D.Re_L     = D.U_w*L_plate/P.nu_w;
D.h_top_corr = 0.664*sqrt(D.Re_L)*P.Pr_w^(1/3)*P.k_w/L_plate;  % laminar plate
if isempty(P.h_top),     D.h_top = D.h_top_corr;  D.h_top_source = 'laminar plate';
elseif isinf(P.h_top),   D.h_top = Inf;           D.h_top_source = 'pinned (Dirichlet)';
else,                    D.h_top = P.h_top;       D.h_top_source = 'manual';
end
D.pinned = isinf(D.h_top);
D.Bi_top = D.h_top*P.H_gel/P.k;

% --------------------------------------------------------------- WALLS --
D.h_wall_series = 1/(1/P.h_air + P.t_wall/P.k_wall);
if isempty(P.h_wall), D.h_wall = D.h_wall_series;  D.h_wall_source = 'air + acrylic in series';
else,                 D.h_wall = P.h_wall;         D.h_wall_source = 'manual';
end
D.Bi_wall = D.h_wall*P.H_gel/P.k;

% ---------------------------------------------------------------- GRID --
dh = P.dh;
D.nx = round((P.x_max - P.x_min)/dh) + 1;  D.x = P.x_min + (0:D.nx-1)*dh;
D.ny = 2*round(P.W_gel/(2*dh)) + 1;        D.y = ((0:D.ny-1) - (D.ny-1)/2)*dh;   % always holds y = 0
D.nz = round(P.H_gel/dh) + 1;              D.z = -((D.nz-1):-1:0)*dh;            % z(end) = 0, the surface
D.N  = D.nx*D.ny*D.nz;
D.iy0 = (D.ny + 1)/2;

% ---------------------------------------------------------- TIME STEP --
g_top = 0;  if ~D.pinned, g_top = 2*dh*D.h_top/P.k; end
g_w   = 2*dh*D.h_wall/P.k;
D.g_top = g_top;  D.g_w = g_w;
D.dt_max = 0.9*dh^2/(D.alpha*(6 + g_top + 2*g_w));
D.npf    = ceil(P.dt_frame/D.dt_max);
D.dt     = P.dt_frame/D.npf;
D.nfr    = round(P.t_end/P.dt_frame) + 1;
D.nsteps = (D.nfr - 1)*D.npf;
D.t_est  = 3.2e-8*D.N*D.nsteps;                          % [s] measured: 61^3 nodes x 3960 steps = 29 s

% ------------------------------------------------------------ CATHETER --
D.gap = P.h_axis - P.b_ins;                              % [m] water between insulation and gel
D.D_c = 2*P.b_ins;                                       % wetted diameter

% ------------------------------------------------------------- CACHE KEY --
D.key = param_key(P);
end

% =========================================================================
function key = param_key(P)
% MD5 of every field of P (name, value, complex parts), 12 hex characters
fn = fieldnames(P);  s = '';
for i = 1:numel(fn)
    v = P.(fn{i});
    if isempty(v)
        s = [s fn{i} '=[];']; %#ok<AGROW>
    else
        s = [s fn{i} '=' sprintf('%.10g%+.10gj,', [real(v(:)) imag(v(:))].') ';']; %#ok<AGROW>
    end
end
md = java.security.MessageDigest.getInstance('MD5');
md.update(uint8(s));
h = typecast(md.digest(), 'uint8');
key = lower(reshape(dec2hex(double(h), 2).', 1, []));
key = key(1:12);
end
