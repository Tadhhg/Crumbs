function R = gui_model(varargin)
%GUI_MODEL  One fast forward solve -- parameters in, SAR map and temperature
%           map out.  This is the engine behind HEAT_GUI.
%
%   R = GUI_MODEL('name',value,...)
%
%   WHY THIS EXISTS.  RUN_HEAT_MODEL solves the baseline once, at 0.25 mm and
%   for 600 s, and takes about a quarter of an hour; everything downstream
%   (T_GEL, HEAT_SCENARIO) then re-uses that single stored solution through
%   the exact power scaling.  That is the right architecture as long as the
%   PHYSICS is fixed -- but it means nothing downstream can change the gel
%   permittivity, the coolant, or the catheter, because those change Theta
%   itself.  This function re-solves everything, from the electromagnetics
%   down, in a few seconds, so those parameters can be swept interactively.
%   It is the same physics and the same code -- SAR_VOLUME, HEAT_GRID,
%   HEAT_LAP, HEAT_SOLVE and CONVECTION_H are all called unchanged -- run on a
%   coarser grid and only as far in time as the exposure lasts.
%   VERIFY_GUI_MODEL measures what that costs against the stored baseline.
%
%   PARAMETERS  (name/value; every one of them is live in the GUI)
%
%     'eps_r'      gel relative permittivity, REAL part          [-]
%     'sigma'      gel effective conductivity                    [S/m]
%                  Whichever of the two is not supplied comes from
%                  SALINE_STOGRYN at the gel temperature.  Passing either one
%                  explicitly turns 'link_eps' off unless you set it yourself,
%                  so gui_model('sigma',1.8) uses 1.8 rather than quietly
%                  overriding it.
%     'T_gel'      INITIAL gel temperature                       [degC]
%     'T_water'    coolant temperature: the sink, and the temperature the
%                  surface is held or driven towards              [degC]
%     'Qflow_Lmin' coolant volumetric flow rate                  [L/min]
%     'L_cath_mm'  length of catheter lying on the surface, measured back
%                  from the antenna tip at x = +12 mm            [mm]
%     'Pgen'       generator FORWARD power                       [W]
%     't_stim'     exposure duration                             [s]
%
%     'link_eps'   true (default) -> eps_r and sigma both come from Stogryn at
%                  T_gel;  false -> the values passed in are used as given.
%                  Turns itself off when eps_r or sigma is passed explicitly.
%     'surf_bc'    'robin' (default here) | 'dirichlet' | 'adiabatic'
%     'cath_bc'    'cooled' (default here) | 'transparent' | 'insulating'
%     'dh'         grid spacing [m]: 1e-3 | 0.5e-3 (default) | 0.25e-3
%     'Nfft'       FFT points per transverse direction, default 512
%     'Nt'         number of stored output times, default 41
%     'cfl'        safety factor on the explicit stability limit, default 0.90.
%                  Forward Euler is first order in time, so this is the knob
%                  that trades wall time against the O(dt) error -- see test 7
%                  of VERIFY_GUI_MODEL, which uses it to show that the residual
%                  in the energy audit IS that error and nothing else.
%     'em'         a previous R.em, re-used if the EM parameters are unchanged
%     'verbose'    print timings
%
%   THE TWO DEFAULTS THAT DIFFER FROM THE WRITE-UP BASELINE, AND WHY
%   The write-up baseline is surf_bc = 'dirichlet' (the surface pinned at
%   37 degC, which is what was specified) with cath_bc = 'transparent'.  Under
%   those two conditions the flow rate and the catheter length have NO EFFECT
%   AT ALL -- an infinite heat sink does not care how fast the coolant moves,
%   and a catheter resting on an already-pinned surface changes nothing.  Both
%   are parameters the GUI is asked to sweep, so the GUI defaults to the
%   physical pair instead: a convective surface with h from CONVECTION_H, and
%   a catheter that acts as a second cooled surface.  Both baselines remain
%   selectable, and the Dirichlet answer is still the lower bound on
%   temperature (write-up Sec. 3.2, 10B).
%
%   OUTPUT R
%     x, y, z       grid coordinates [m]  (z <= 0 in the gel)
%     t             output times [s]
%     SAR           [nx x nz] time-averaged SAR on the y = 0 plane [W/kg]
%     q_y0          the same as volumetric deposition [W/m^3]
%     T_plane       [nx x nz x Nt] temperature on the y = 0 plane [degC]
%     T_max         peak temperature vs time, feed sphere excluded [degC]
%     T_max_all     the same including it (not a prediction -- Sec. 6.2)
%     xyz_max       location of the masked peak at each time [m]
%     pmask_y0      false where the feed singularity is blanked
%     em            the cached EM solution -- pass it back in as 'em'
%     info          every derived quantity worth displaying (h, eps, lambda,
%                   penetration depth, absorbed power, timings, flags)
%
%   See also HEAT_GUI, VERIFY_GUI_MODEL, RUN_HEAT_MODEL, HEAT_SCENARIO.

addpath(fileparts(mfilename('fullpath')));   % the rest of the heating model
hm_addpath();                                % and the field model, wherever it is

% ------------------------------------------------------------- options ---
% T_gel defaults to the coolant temperature, so that out of the box the gel
% starts in equilibrium and the exact power scaling of Sec. 11.1 holds.  Set it
% lower for a phantom brought straight from the bench.
o = struct('eps_r',[], 'sigma',[], 'T_gel',37, 'T_water',37, ...
           'Qflow_Lmin',5, 'L_cath_mm',32, 'Pgen',12, 't_stim',60, ...
           'link_eps',true, 'surf_bc','robin', 'cath_bc','cooled', ...
           'dh',0.5e-3, 'Nfft',512, 'Nt',41, 'cfl',0.90, 'em',[], ...
           'verbose',false);
known = fieldnames(o);
given = false(size(known));
for a = 1:2:numel(varargin)
    nm = varargin{a};
    if ~ischar(nm) && ~isstring(nm)
        error('gui_model:opt','option %d is not a name', (a+1)/2);
    end
    j = find(strcmpi(known, nm), 1);
    if isempty(j)
        error('gui_model:opt','unknown option ''%s''.  Known: %s', ...
              char(nm), strjoin(known(:).', ', '));
    end
    o.(known{j}) = varargin{a+1};
    given(j) = true;
end

% If a permittivity or a conductivity is passed explicitly, the caller means
% it.  link_eps defaults to true, so without this a call like
%     gui_model('sigma',1.8)
% would silently discard the 1.8 and use Stogryn instead -- the worst kind of
% surprise, since the run would look fine and answer a different question.
% An explicit link_eps still wins, both ways.
if ~given(strcmp(known,'link_eps')) && ...
        (given(strcmp(known,'eps_r')) || given(strcmp(known,'sigma')))
    o.link_eps = false;
end

wall = tic;
eps0 = 8.8541878128e-12;                        % [F/m] CODATA 2018
f    = 2.45e9;                                  % [Hz]  fixed by the brief
w    = 2*pi*f;

% ------------------------------------------------- the gel's permittivity --
% Two routes, and the GUI can use either.  Linked: the gel is 0.9 % saline, so
% Stogryn gives eps' and sigma at whatever temperature it sits at.  Unlinked:
% eps' and sigma are whatever the user types -- the honest option once the gel
% has actually been measured with a coaxial probe (write-up Sec. 4, 11).
% Whatever is not supplied comes from Stogryn, so a caller can pin one of the
% two and let the other follow the temperature.
Nnorm  = 0.9*10/58.443;                         % [mol/L] 0.9 % w/v NaCl
use_eg = o.link_eps || isempty(o.eps_r) || isempty(o.sigma);
if use_eg
    Ts = o.T_gel;
    if Ts < 0 || Ts > 40
        warning('gui_model:stogryn', ...
            ['Stogryn (1971) is fitted for 0-40 degC; %.1f degC is an ' ...
             'extrapolation.  Enter eps_r and sigma by hand instead.'], Ts);
        Ts = min(max(Ts,0),40);
    end
    [eg, ginfo] = saline_stogryn(f, Ts, Nnorm);
end
if o.link_eps
    eps_r = real(eg);   sigma = ginfo.sigma_eff_Spm;
    eps_src = sprintf('Stogryn, 0.9%% NaCl at %.1f degC', o.T_gel);
else
    if isempty(o.eps_r), eps_r = real(eg);            else, eps_r = o.eps_r; end
    if isempty(o.sigma), sigma = ginfo.sigma_eff_Spm; else, sigma = o.sigma; end
    src = {'given','given'};
    if isempty(o.eps_r), src{1} = 'Stogryn'; end
    if isempty(o.sigma), src{2} = 'Stogryn'; end
    if strcmp(src{1},src{2}) && strcmp(src{1},'given')
        eps_src = 'entered by hand';
    else
        eps_src = sprintf('eps'' %s, sigma %s', src{1}, src{2});
    end
end
if eps_r <= 1 || sigma < 0
    error('gui_model:eps', ...
        'need eps_r > 1 and sigma >= 0; got eps_r = %g, sigma = %g', eps_r, sigma);
end
epsr_gel = eps_r - 1j*sigma/(w*eps0);

% --------------------------------------------------- the thermal settings --
% water_properties is tabulated over 16.85..96.85 degC and CONVECTION_H needs
% it, so the coolant temperature is bounded by the property table, not by the
% physics.  The bottom of the table exists so that room-temperature coolant
% (~20-25 degC) is inside it, not an extrapolation.
if o.T_water < 17 || o.T_water > 96
    error('gui_model:Twater', ...
        ['coolant temperature %.1f degC is outside the tabulated water ' ...
         'property range 17..96 degC'], o.T_water);
end
if o.Qflow_Lmin <= 0
    error('gui_model:Q','flow rate must be positive (got %g L/min)', o.Qflow_Lmin);
end
x_tip = 12e-3;                                  % the catheter runs BACK from
                                                % the antenna tip
Tc = thermal_config('T_blood', o.T_water, ...
                    'T_init',  o.T_gel, ...
                    'Qflow',   o.Qflow_Lmin/60000, ...
                    'surf_bc', o.surf_bc, ...
                    'cath_bc', o.cath_bc, ...
                    'cath_x0', x_tip - o.L_cath_mm*1e-3, ...
                    'cath_L',  x_tip, ...
                    'cfl',     o.cfl, ...
                    'dh',      o.dh);
Pdel = Tc.eta*o.Pgen;                           % [W] delivered to the gel

xg = (Tc.xlim(1):o.dh:Tc.xlim(2)+1e-12).';
yg = (Tc.ylim(1):o.dh:Tc.ylim(2)+1e-12).';
zg = -(0:o.dh:Tc.d+1e-12).';

% ================= 1. the electromagnetic source ==========================
% q is linear in the delivered power, so the solve is done ONCE at 1 W
% delivered and cached; changing the power afterwards is a multiplication.
% The cache key is every input the EM solve actually depends on.
key = [eps_r, sigma, o.dh, o.Nfft, f];
if ~isempty(o.em) && isequal(o.em.key, key)
    em   = o.em;
    t_em = 0;
else
    t0 = tic;
    C = model_config('Pdel', 1);
    cfg = struct('f',C.f, 'epsr',[C.epsr_air, epsr_gel, C.epsr_air], ...
                 'd',C.d, 'h',C.h, 'L',C.L, 'ke',C.ke, 'Idist',C.Idist, ...
                 'Nfft',o.Nfft, 'dx',C.dx, 'Pdel',1, ...
                 'zn',zg, 'xlim',Tc.xlim, 'ylim',Tc.ylim);
    s = sar_volume(cfg);
    % the EM solve always runs on the 0.25 mm FFT grid; subsample onto the
    % thermal grid, exactly as RUN_HEAT_MODEL does
    step = round(o.dh/C.dx);
    if step > 1
        s.q = s.q(1:step:end, 1:step:end, :);
        s.x = s.x(1:step:end);  s.y = s.y(1:step:end);
    end
    assert(numel(s.x)==numel(xg) && numel(s.y)==numel(yg) && ...
           size(s.q,3)==numel(zg), 'EM block does not match the thermal grid');
    em = struct('key',key, 'q_unit',s.q, 'Pgel',s.Pgel, 'Pup',s.Pup, ...
                'Pdn',s.Pdn, 'Pblock',s.Pblock, 'sigma1',s.sigma1, ...
                'k1',s.k(2), 'L',C.L, 'hwire',C.h);
    t_em = toc(t0);
end
q = em.q_unit*Pdel;                             % [W/m^3] at this power

% the feed singularity: mask it out of the peak reporting and out of the maps
[q, finfo] = feed_regularise(q, xg, yg, zg, Tc);

% ================= 2. the transient conduction solve ======================
t0 = tic;
G = heat_grid(Tc, xg, yg, zg);
tmin  = max(o.t_stim/1e4, 1e-3);
t_out = unique([0, logspace(log10(tmin), log10(o.t_stim), o.Nt-1)]);
S = heat_solve(q, G, Tc, struct('t_out',t_out, 'verbose',false, 'quiet',true, ...
                                'decim',max(1,round(2e-3/o.dh)), 'Nt3',2, ...
                                'peakmask',finfo.mask_peak, ...
                                'theta0',o.T_gel - o.T_water));
t_ht = toc(t0);

% ================= 3. assemble what the GUI draws =========================
iy0 = S.iy0;
R = struct();
R.x = xg;  R.y = yg;  R.z = zg;  R.t = S.t;
R.q_y0 = squeeze(em.q_unit(:,iy0,:))*Pdel;      % [W/m^3], unregularised
R.SAR  = R.q_y0/Tc.rho_gel;                     % [W/kg]
R.T_plane   = o.T_water + double(S.Th_y0);      % [degC]
R.T_max     = o.T_water + S.Th_max;
R.T_max_all = o.T_water + S.Th_max_all;
R.xyz_max   = [xg(S.ij_max(1,:)).'; yg(S.ij_max(2,:)).'; zg(S.ij_max(3,:)).'];
R.pmask_y0  = squeeze(finfo.mask_peak(:,iy0,:));
R.em  = em;
R.Tc  = Tc;
R.opt = o;
% the raw material of the energy audit, so VERIFY_GUI_MODEL can take it apart
R.audit = struct('t',S.t, 'U',S.U, 'P_loss',S.P_loss, 'P_in',S.P_in, ...
                 'err',S.err_energy);

% ------------------------------------------------------------ diagnostics -
k1 = em.k1;
i = struct();
i.eps_r      = eps_r;
i.sigma      = sigma;
i.eps_src    = eps_src;
i.epsr_gel   = epsr_gel;
i.lambda_gel = 2*pi/real(k1);                   % [m] wavelength in the gel
i.delta_pw   = 1/abs(imag(k1));                 % [m] plane-wave skin depth
i.h_surf     = Tc.h_surf;
i.h_cath     = Tc.h_cath;
i.U_flow     = Tc.conv.U;
i.ReL        = Tc.conv.ReL;
i.dT_coolant = Tc.conv.dT_per_W*Pdel;           % [K] bulk coolant rise
i.Pdel       = Pdel;
i.eta        = Tc.eta;
i.Pgel       = em.Pgel*Pdel;                    % [W] absorbed in the whole gel
i.frac_gel   = em.Pgel/(em.Pgel + em.Pup + em.Pdn);
i.frac_block = em.Pblock/em.Pgel;
i.P_in       = S.P_in;                          % [W] heating the free nodes
i.SAR_max    = max(R.SAR(R.pmask_y0));
i.T_end      = R.T_max(end);
i.T_end_all  = R.T_max_all(end);
i.xyz_end    = R.xyz_max(:,end).';

% The masked peak usually sits ON the feed-exclusion sphere, because the
% deposition still rises towards the feed just outside it (write-up Sec. 6.2,
% 9.2).  It is the right number to quote as a bound, but it is not a number a
% thermocouple can be put on, so report the canonical bench point as well:
% (6, 0, -2) mm, under the middle of the element, which is what Sec. 11.5 uses
% for comparing against a measurement.
[~,ix6] = min(abs(xg - 6e-3));
[~,iz2] = min(abs(zg + 2e-3));
i.probe_xyz = [xg(ix6) 0 zg(iz2)];
i.T_probe   = squeeze(R.T_plane(ix6,iz2,:)).';
i.T_probe_end = i.T_probe(end);
i.err_energy = S.err_energy;
i.dt         = S.dt_max;
i.nstep      = S.nstep;
i.t_em       = t_em;
i.t_heat     = t_ht;
i.t_total    = toc(wall);
i.nodes      = G.nx*G.ny*G.nz;

% How deep the heating reaches, under the middle of the element.  Two
% definitions, because the obvious one is easy to misread:
%
%   d_1e  the depth at which q has fallen to 1/e of its SURFACE value.  Note
%         this is q, not |E|.  Since q = 0.5 sigma |E|^2, the 1/e depth of q is
%         HALF that of |E| for a decaying exponential -- so this number is
%         about half the "field falls to 1/e by 2.7 mm" of the field write-up,
%         and they do not contradict each other.
%   d50, d90  the depths above which half and nine tenths of the power
%         deposited in that vertical column lie.  These need no reference
%         value, so they are the honest answer to "how deep does it heat",
%         and they match the framing the field write-up uses (Sec. 9.4).
prof   = R.q_y0(ix6,:).';                       % ix6 found above
i.d_1e = first_crossing(-zg, prof, prof(1)/exp(1));
cum    = cumtrapz(-zg, prof);
if cum(end) > 0
    frac  = cum/cum(end);
    [fu, iu] = unique(frac);
    i.d50 = interp1(fu, -zg(iu), 0.50);
    i.d90 = interp1(fu, -zg(iu), 0.90);
else
    i.d50 = NaN;  i.d90 = NaN;
end

% ------------------------------------------------------------- the flags --
fl = {};
if i.T_end > 100
    fl{end+1} = sprintf(['peak %.0f degC is above boiling -- there is no ' ...
        'phase change in this model, so that is a warning, not a temperature'], ...
        i.T_end);
elseif i.T_end > 60
    fl{end+1} = sprintf(['peak %.0f degC: above 60 degC the constant-property ' ...
        'assumption costs 10-15 %%, and it UNDER-predicts (Sec. 10K)'], i.T_end);
end
if i.T_end > 50
    fl{end+1} = ['above ~50 degC low-acyl gellan may be softening -- the ' ...
                 'phantom may stop being a phantom (Sec. 12.5)'];
end
if strcmpi(o.surf_bc,'dirichlet')
    fl{end+1} = ['Dirichlet surface: the flow rate has NO effect (perfect ' ...
                 'sink).  This is the lower bound on temperature.'];
end
if strcmpi(o.cath_bc,'transparent')
    fl{end+1} = 'transparent catheter: the catheter length has NO effect';
end
if i.dT_coolant > 1
    fl{end+1} = sprintf(['the coolant itself would warm by %.1f K at this ' ...
        'power -- the constant-sink assumption is getting thin (Sec. 3.3)'], ...
        i.dT_coolant);
end
if abs(o.T_gel - o.T_water) > 0.05
    fl{end+1} = sprintf(['the gel starts %.1f K from the coolant, so the ' ...
        'first seconds are a relaxation as well as a heating transient'], ...
        o.T_gel - o.T_water);
end
% Judge conservation on the balance at the FINAL time, not on S.err_energy.
% S.err_energy is the worst over all output times, and the worst is always one
% of the earliest samples, where the stored energy it divides by is still
% almost nothing -- it reads as 1e-3 on a perfectly healthy run.  The final
% balance is the meaningful test.  VERIFY_GUI_MODEL test 7 separates the two
% things that limit it: the trapezoid used on P_loss between output times, and
% the O(dt) error of forward Euler, which is the one that dominates here.
Ubal_end = S.U(1) + S.P_in*S.t(end) - trapz(S.t, S.P_loss);
if abs(Ubal_end) > 0
    i.err_energy_final = abs(S.U(end) - Ubal_end)/abs(Ubal_end);
else
    i.err_energy_final = 0;
end
if i.err_energy_final > 5e-3
    fl{end+1} = sprintf(['the energy balance closes only to %.1e at the end ' ...
        'of the run -- reduce the time step (cfl) or the grid'], i.err_energy_final);
end
i.flags = fl;
R.info  = i;

if o.verbose
    fprintf(['gui_model: %.4g W for %.4g s, eps_r %.2f, sigma %.3f S/m, ' ...
             'h %.0f W/(m^2 K)\n'], o.Pgen, o.t_stim, eps_r, sigma, Tc.h_surf);
    fprintf('  peak %.1f degC at (%.1f, %.1f, %.2f) mm; SAR max %.3e W/kg\n', ...
            i.T_end, 1e3*i.xyz_end(1), 1e3*i.xyz_end(2), -1e3*i.xyz_end(3), ...
            i.SAR_max);
    fprintf('  EM %.2f s + heat %.2f s = %.2f s (%d steps, %.2f M nodes)\n', ...
            t_em, t_ht, i.t_total, i.nstep, i.nodes/1e6);
end
end

% =========================================================================
function d = first_crossing(depth, prof, level)
%FIRST_CROSSING  Depth at which prof first falls below level, interpolated.
j = find(prof < level, 1);
if isempty(j) || j == 1, d = NaN; return; end
d = interp1(prof([j-1 j]), depth([j-1 j]), level);
end
