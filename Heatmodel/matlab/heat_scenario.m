function S = heat_scenario(P, period, t_stim, varargin)
%HEAT_SCENARIO  Run one bench scenario and draw the temperature profile.
%
%   S = HEAT_SCENARIO(P, period, t_stim)
%
%     P       generator FORWARD power while the microwave is on          [W]
%     period  pulse repetition period [s].  Use [] (or Inf, or 0) for
%             continuous wave.  If a period is given and no duty cycle is,
%             50 %% is assumed and said so on the console.
%     t_stim  total stimulation time [s] -- how long the run lasts
%
%   Draws two figures and writes a CSV of the traces, so a measured
%   thermocouple or fibre-optic record can be overlaid on the prediction.
%
%   OPTIONS  (name/value)
%     'duty'      duty cycle in (0,1]                       default 0.5 if a
%                                                           period is given
%     'ton'       on-time per period [s]                    (instead of duty)
%     'model'     'dirichlet' (default, the surface held at 37 degC as
%                 specified) | 'robin' (the convective surface, h = 1460)
%     'probes'    n x 3 matrix of (x,y,z) probe positions [m].  Default is a
%                 depth ladder at x = 6 mm: 1, 2, 5 and 10 mm deep.
%     'measured'  [t, T1, T2, ...] measured record to overlay: first column
%                 time [s], remaining columns degC, one per probe
%     'tview'     how long to plot for [s]                  default 1.5*t_stim
%     'name'      tag used in the output file names
%     'delivered' true -> P is the DELIVERED power, not the generator setting
%     'save'      false -> draw but write nothing            default true
%
%   EXAMPLES
%     heat_scenario(12, [], 60)                     % 12 W continuous, 60 s
%     heat_scenario(40, 1, 60, 'duty', 0.3)         % 40 W, 30 %% at 1 Hz, 60 s
%     heat_scenario(12, [], 60, 'model','robin')    % convective surface
%     heat_scenario(12, [], 60, 'probes', [6e-3 0 -2e-3; 6e-3 0 -6e-3])
%     heat_scenario(12, [], 60, 'measured', myLog)  % overlay a bench record
%
%   HOW IT IS COMPUTED -- no new solve is needed
%   The heat equation is linear and time invariant and the deposition is
%   proportional to the delivered power, so the response to ANY piecewise-
%   constant power history is a superposition of the stored step response
%   Theta(r,t):
%
%       theta(r,t) = (1/P_ref) SUM_j  dP_j * Theta(r, t - tau_j)
%
%   where the drive switches by dP_j at time tau_j.  For continuous wave that
%   is exactly two terms -- switch on at t = 0, switch off at t = t_stim --
%   so the cool-down after the exposure is free and exact.  For a pulse train
%   the same sum is evaluated as a discrete convolution on a uniform grid.
%   Everything is therefore as accurate as the single solve in
%   RUN_HEAT_MODEL, and takes a few seconds.
%
%   WHAT TO COMPARE ON THE BENCH
%   The CSV holds the drive waveform, the peak temperature and every probe
%   trace.  Two warnings before you compare:
%     * The BASELINE model holds the gel surface at exactly 37 degC.  That is
%       what was specified, and it is an optimistic bound: it under-predicts.
%       Run 'model','robin' as well and expect the measurement to fall between
%       the two (write-up Sec. 3.2, 10B).
%     * A probe within 2 mm of the feed (x = 0, z = 0) sits in the feed
%       singularity of the field model and is not predictable (Sec. 6.2).
%       Put your probes under the middle of the element, not behind it.
%
%   See also T_GEL, P_FOR_T, TIME_TO_T, RUN_HEAT_MODEL.

here = fileparts(mfilename('fullpath'));
root = fileparts(here);

% ------------------------------------------------------------- options ---
o = struct('duty',[], 'ton',[], 'model','dirichlet', 'probes',[], ...
           'measured',[], 'tview',[], 'name','', 'delivered',false, ...
           'save',true, 'nmax_trace',8000, 'nmax_map',700);

% ---- check the three positional arguments before anything else ----------
% The commonest mistake is to slide an option in front of t_stim, e.g.
%   heat_scenario(20, [], 'model','robin', 360)      <- WRONG
%   heat_scenario(20, [], 360, 'model','robin')      <- right
% Without this check that lands a char array in t_stim and fails much later
% with an unrelated message about array sizes.
if ischar(P) || isstring(P) || ~isnumeric(P) || ~isscalar(P) || ~isfinite(P) || P < 0
    error('heat_scenario:P', ['first argument must be the generator power in ' ...
        'watts (a non-negative scalar).\n  Usage: heat_scenario(P, period, ' ...
        't_stim, ''name'',value, ...)']);
end
if ~isempty(period) && (ischar(period) || isstring(period) || ...
                        ~isnumeric(period) || ~isscalar(period))
    error('heat_scenario:period', ['second argument must be the pulse period ' ...
        'in seconds, or [] for continuous wave.\n  Usage: heat_scenario(P, ' ...
        'period, t_stim, ''name'',value, ...)']);
end
if ischar(t_stim) || isstring(t_stim)
    error('heat_scenario:t_stim', ['third argument must be the stimulation ' ...
        'time in seconds, but got the option name ''%s''.\n' ...
        '  The three positional arguments come FIRST, options after:\n' ...
        '      heat_scenario(%g, %s, t_stim, ''%s'', ...)\n' ...
        '  e.g. heat_scenario(20, [], 360, ''model'', ''robin'')'], ...
        char(t_stim), P, ternary(isempty(period),'[]',num2str(period)), char(t_stim));
end
if ~isnumeric(t_stim) || ~isscalar(t_stim) || ~isfinite(t_stim) || t_stim <= 0
    error('heat_scenario:t_stim', ['third argument must be the stimulation ' ...
        'time in seconds (a positive scalar).\n  Usage: heat_scenario(P, ' ...
        'period, t_stim, ''name'',value, ...)']);
end

% ---- then the name/value options ----------------------------------------
known = fieldnames(o);
if mod(numel(varargin), 2) ~= 0
    error('heat_scenario:pairs', ['options must come in name/value pairs; ' ...
        'got %d trailing argument(s).\n  Usage: heat_scenario(P, period, ' ...
        't_stim, ''name'',value, ...)'], numel(varargin));
end
for q = 1:2:numel(varargin)
    nm = varargin{q};
    if ~(ischar(nm) || isstring(nm))
        error('heat_scenario:optname', ...
            'option %d is not a name; expected one of: %s', ...
            (q+1)/2, strjoin(known.', ', '));
    end
    nm = char(nm);
    hit = known(strcmpi(known, nm));
    if isempty(hit)
        error('heat_scenario:option', ['unknown option ''%s''.  Valid ' ...
            'options are:\n    %s'], nm, strjoin(known.', ', '));
    end
    o.(hit{1}) = varargin{q+1};
end

M = load_model(o.model, root);
Tb = M.T_blood;  eta = M.eta;  Pref = M.Pdel_ref;

if isempty(o.probes)
    o.probes = [6e-3 0 -1e-3; 6e-3 0 -2e-3; 6e-3 0 -5e-3; 6e-3 0 -10e-3];
end
if isempty(o.tview), o.tview = 1.5*t_stim; end
if o.tview > M.t(end)
    warning('heat_scenario:tview', ...
        'tview %.0f s exceeds the computed history (%.0f s); clamped.', ...
        o.tview, M.t(end));
    o.tview = M.t(end);
end

% ------------------------------------------------- the drive waveform ----
cw = isempty(period) || ~isfinite(period) || period <= 0;
if cw
    duty = 1;  period = Inf;
else
    if ~isempty(o.ton)
        duty = o.ton/period;
    elseif ~isempty(o.duty)
        duty = o.duty;
    else
        duty = 0.5;
        fprintf(['  NOTE: a pulse period was given but no duty cycle; ' ...
                 'assuming 50 %%.\n        Pass ''duty'' or ''ton'' to set it.\n']);
    end
    if duty <= 0 || duty > 1, error('heat_scenario:duty','duty must be in (0,1]'); end
    if duty == 1, cw = true; period = Inf; end
end
Pdel_on = P;  if ~o.delivered, Pdel_on = eta*P; end
sON     = Pdel_on/Pref;                       % scale factor while on
Pmean   = P*duty;

if isempty(o.name)
    if cw
        o.name = sprintf('P%gW_cw_%gs', P, t_stim);
    else
        o.name = sprintf('P%gW_%gHz_%gpc_%gs', P, 1/period, round(100*duty), t_stim);
    end
    if strcmpi(o.model,'robin'), o.name = [o.name '_robin']; end
end
o.name = regexprep(o.name, '[^A-Za-z0-9_.-]', '_');

% ------------------------------------------------------ the time bases ---
% TWO grids, for two different jobs.  The temperature MAPS live on a coarse
% grid, because a 121 x 121 plane at many times is a large dense object.  The
% TRACES -- peak and probes -- live on a fine grid, because a pulse train has
% to be resolved within each pulse or the peak is smeared.  For continuous
% wave the two coincide and the response is the exact two-step superposition
% (switch on at t = 0, switch off at t_stim), so no grid is needed at all.
if cw
    % the stored log grid resolves the first millisecond; add the switch-off
    tv = unique([M.t(M.t <= o.tview), o.tview, t_stim, ...
                 t_stim + [0 1e-3 1e-2 0.1 0.3 1 3 10 30 100]]);
    tv = tv(tv >= 0 & tv <= o.tview);
    tmap = tv;
else
    ton  = duty*period;  toff = (1-duty)*period;
    dtT  = min([ton/20, toff/20, o.tview/400]);
    NT   = min(o.nmax_trace, max(400, ceil(o.tview/dtT)));
    tv   = (0:NT)*(o.tview/NT);
    NP   = min(o.nmax_map, NT);
    tmap = (0:NP)*(o.tview/NP);
end
Pdrive = drive_power(tv, P, period, duty, t_stim);
pfrac  = drive_power(tmap, P, period, duty, t_stim);
if ~o.delivered, sc = eta/Pref; else, sc = 1/Pref; end

% ------------------------------------------ kernels: plane and probes ----
nx = numel(M.x);  nz = numel(M.z);  npt = nx*nz;
Kplane = reshape(double(M.Th_y0), npt, []);            % [npt x Nt]
np = size(o.probes,1);
Kprobe = zeros(np, numel(M.t));
for r = 1:np
    Kprobe(r,:) = probe_kernel(M, o.probes(r,:));
end

% ------------------------------------- pass 1: the plane, coarse in time --
tic0 = tic;
if cw
    Yp = sON*( kern_at(M.t, Kplane, tmap) - kern_at(M.t, Kplane, tmap - t_stim) );
else
    Yp = duhamel_matrix(M.t, Kplane, tmap, pfrac*sc);
end
Tplane = reshape(Tb + Yp, nx, nz, numel(tmap));
clear Yp

% peak over the plane at the map times, feed sphere excluded, to find every
% location that is ever the hot spot
Tm = reshape(Tplane, npt, []);
Tm(~M.pmask_y0(:), :) = -Inf;
[~, licand] = max(Tm, [], 1);
cand = unique(licand);
clear Tm

% ------------------------------- pass 2: candidates and probes, fine in time
Ktr = [Kplane(cand,:); Kprobe];
if cw
    Ytr = sON*( kern_at(M.t, Ktr, tv) - kern_at(M.t, Ktr, tv - t_stim) );
else
    Ytr = duhamel_matrix(M.t, Ktr, tv, Pdrive*sc);
end
Ttr = Tb + Ytr;
[Tpk, ib] = max(Ttr(1:numel(cand),:), [], 1);
[ipx, ipz] = ind2sub([nx nz], cand(ib));
Tprobe = Ttr(numel(cand)+1:end, :);
walltime = toc(tic0);

% thermal dose, on the map grid
Dose = cem43(tmap, reshape(Tplane, npt, []));
Dose = reshape(Dose, nx, nz);
Dose(~M.pmask_y0) = 0;

% ---------------------------------------------------------- summary -----
[Tmax, imx] = max(Tpk);
S = struct();
S.t = tv(:);  S.P_drive = Pdrive(:);  S.T_peak = Tpk(:);
S.T_probe = Tprobe.';  S.probes = o.probes;
S.t_map = tmap(:);  S.n_hotspots = numel(cand);
S.T_plane = Tplane;  S.x = M.x;  S.z = M.z;  S.CEM43 = Dose;
S.P_generator = P;  S.P_delivered = Pdel_on;  S.duty = duty;
S.period = period;  S.t_stim = t_stim;  S.model = o.model;
S.P_mean_generator = Pmean;  S.energy_delivered = Pdel_on*duty*t_stim;
S.T_max = Tmax;  S.t_at_T_max = tv(imx);
S.hotspot = [M.x(ipx(imx)) 0 M.z(ipz(imx))];
S.name = o.name;

print_summary(S, M, o, cw, walltime);

% ------------------------------------------------------------ output ----
make_scenario_figures(S, M, o, cw);
if o.save
    write_csv(S, root);
    fprintf('  wrote figures/scenario_%s_traces.png, _maps.png\n', o.name);
    fprintf('        results/scenario_%s.csv\n', o.name);
end
end

% =========================================================================
function M = load_model(name, root)
persistent cache
key = matlab.lang.makeValidName(lower(name));
if ~isempty(cache) && isfield(cache,key), M = cache.(key); return; end
switch lower(name)
    case 'dirichlet', f = 'heat_model.mat';
    case 'robin',     f = 'heat_model_robin.mat';
    otherwise, error('heat_scenario:model','model must be dirichlet or robin');
end
fp = fullfile(root,'results',f);
if ~exist(fp,'file')
    error('heat_scenario:noModel','%s not found. Run run_heat_model first.', fp);
end
M = load(fp);
M.Th3 = double(M.Th3);  M.Th_ss = double(M.Th_ss);
cache.(key) = M;
end

function p = drive_power(t, P, period, duty, t_stim)
%DRIVE_POWER  Generator forward power at time t [W].
p = zeros(size(t));
on = t >= 0 & t < t_stim;
if ~isfinite(period)
    p(on) = P;
else
    p(on) = P*(mod(t(on), period) < duty*period);
end
end

function Y = kern_at(tk, K, tq)
%KERN_AT  K interpolated in time, zero for tq <= 0.
tq = tq(:).';
Y = zeros(size(K,1), numel(tq));
m = tq > 0;
if any(m)
    Y(:,m) = interp1(tk, K.', min(tq(m), tk(end)), 'linear').';
    if nnz(m) == 1, Y(:,m) = reshape(Y(:,m), [], 1); end
end
end

function Y = duhamel_matrix(tk, K, tu, pgen_over_Pref)
%DUHAMEL_MATRIX  theta on a uniform grid for a piecewise-constant drive.
%   theta(t_n) = SUM_{m=1..n} pbar_m * [Theta(t_{n-m+1}) - Theta(t_{n-m})]
%   evaluated as one dense matrix product, which is exact for a drive that is
%   constant on each interval of the grid.
N  = numel(tu) - 1;
Ku = interp1(tk, K.', min(tu, tk(end)), 'linear');      % [(N+1) x nrow]
dK = diff(Ku, 1, 1);                                    % [N x nrow]
pm = 0.5*(pgen_over_Pref(1:end-1) + pgen_over_Pref(2:end));
C  = tril(toeplitz(pm(:), [pm(1) zeros(1,N-1)]));       % [N x N]
Y  = [zeros(size(K,1),1), (C*dK).'];                    % [nrow x (N+1)]
end

function k = probe_kernel(M, xyz)
%PROBE_KERNEL  Theta(t) at one point, on the stored time base.
x = xyz(1); y = xyz(2); z = xyz(3);
if any(x < M.x(1)-1e-9) || any(x > M.x(end)+1e-9) || ...
   any(y < M.y(1)-1e-9) || any(y > M.y(end)+1e-9) || z > 1e-9 || z < M.z(end)-1e-9
    error('heat_scenario:probe', ...
        'probe (%.1f, %.1f, %.1f) mm is outside the region of interest', 1e3*xyz);
end
if abs(y) < 1e-9
    % full resolution from the stored y = 0 plane
    k = interp2_plane(M, x, -z);
else
    k = zeros(1, numel(M.t3));
    for m = 1:numel(M.t3)
        k(m) = interpn(M.x3(:).', M.y3(:).', -M.z3(:).', M.Th3(:,:,:,m), ...
                       x, y, -z, 'linear');
    end
    k = interp1(M.t3, k, M.t, 'linear', 'extrap');
end
end

function k = interp2_plane(M, xq, sq)
%INTERP2_PLANE  bilinear in (x, depth) at every stored time, full resolution.
xg = M.x(:).';  sg = -M.z(:).';
i0 = max(1, min(numel(xg)-1, find(xg <= xq, 1, 'last')));
j0 = max(1, min(numel(sg)-1, find(sg <= sq, 1, 'last')));
a = (xq - xg(i0))/(xg(i0+1)-xg(i0));
b = (sq - sg(j0))/(sg(j0+1)-sg(j0));
T = double(M.Th_y0);
k = (1-a)*(1-b)*squeeze(T(i0,j0,:)).' + a*(1-b)*squeeze(T(i0+1,j0,:)).' + ...
    (1-a)*b    *squeeze(T(i0,j0+1,:)).' + a*b    *squeeze(T(i0+1,j0+1,:)).';
end

function print_summary(S, M, o, cw, walltime)
fprintf('\n=========== SCENARIO: %s ===========\n', S.name);
if cw
    fprintf('  drive       : %.3f W generator, CONTINUOUS, for %.4g s\n', ...
            S.P_generator, S.t_stim);
else
    fprintf('  drive       : %.3f W generator, %.1f %% duty at %.4g s period ', ...
            S.P_generator, 100*S.duty, S.period);
    fprintf('(%.4g Hz), for %.4g s\n', 1/S.period, S.t_stim);
    fprintf('                mean generator power %.3f W\n', S.P_mean_generator);
end
fprintf('  delivered   : %.3f W while on (eta = %.3f), %.1f J total\n', ...
        S.P_delivered, M.eta, S.energy_delivered);
fprintf('  surface     : %s\n', ternary(strcmpi(o.model,'robin'), ...
        'CONVECTIVE, h = 1460 W/(m^2 K)', 'held at 37 degC (as specified; a lower bound)'));
fprintf('  computed in %.1f s by superposition -- no new solve\n', walltime);
fprintf('                traces on %d time points, maps on %d; %d distinct hot spots\n', ...
        numel(S.t), numel(S.t_map), S.n_hotspots);

fprintf('\n  peak temperature %.2f degC at t = %.4g s, at (%.2f, %.2f, %.2f) mm\n', ...
        S.T_max, S.t_at_T_max, 1e3*S.hotspot(1), 1e3*S.hotspot(2), 1e3*S.hotspot(3));
fprintf('  temperature at the end of stimulation (t = %.4g s):\n', S.t_stim);
[~,ie] = min(abs(S.t - S.t_stim));
fprintf('     peak %.2f degC\n', S.T_peak(ie));
for r = 1:size(S.probes,1)
    fprintf('     probe %d  (%5.1f, %4.1f, %6.2f) mm : %6.2f degC  (max %6.2f degC)\n', ...
        r, 1e3*S.probes(r,1), 1e3*S.probes(r,2), 1e3*S.probes(r,3), ...
        S.T_probe(ie,r), max(S.T_probe(:,r)));
end

d240 = zone(S.CEM43, 240, S.x, S.z);
d60  = zone(S.CEM43,  60, S.x, S.z);
fprintf('\n  thermal dose (tissue equivalent, CEM43):\n');
fprintf('     >= 240 min : %s\n', d240);
fprintf('     >=  60 min : %s\n', d60);

fprintf('\n  validity:\n');
if S.T_max > 100
    fprintf('     *** T exceeds 100 degC -- the gel boils and the model does not know it\n');
end
if S.T_max > 60
    fprintf('     *  T exceeds 60 degC -- permittivity and conductivity drift ~10 %%,\n');
    fprintf('        so the true rise is 10-15 %% ABOVE this (write-up Sec. 10K)\n');
end
if S.T_max > 50
    fprintf('     *  T exceeds 50 degC -- check the gel-sol transition of the Phytagel\n');
end
if strcmpi(o.model,'dirichlet')
    fprintf('     *  surface pinned at 37 degC: this is a LOWER BOUND.  Re-run with\n');
    fprintf('        ''model'',''robin'' and expect the measurement between the two.\n');
end
if S.T_max <= 50
    fprintf('     nothing above 50 degC; the constant-property model is comfortable here\n');
end
end

function s = zone(D, lev, x, z)
m = D >= lev;
if ~any(m(:)), s = 'none'; return; end
[ii,jj] = find(m);
s = sprintf('%.2f mm deep, x from %.1f to %.1f mm (%.1f mm long)', ...
    max(-z(jj))*1e3, min(x(ii))*1e3, max(x(ii))*1e3, ...
    (max(x(ii))-min(x(ii)))*1e3);
end

function v = ternary(c,a,b), if c, v = a; else, v = b; end, end

function write_csv(S, root)
f = fullfile(root,'results',sprintf('scenario_%s.csv', S.name));
fid = fopen(f,'w');
fprintf(fid, '# scenario %s\n', S.name);
fprintf(fid, '# generator %.4f W, duty %.4f, period %.6g s, stimulation %.6g s, model %s\n', ...
        S.P_generator, S.duty, S.period, S.t_stim, S.model);
fprintf(fid, '# delivered %.4f W while on, %.2f J total\n', S.P_delivered, S.energy_delivered);
fprintf(fid, 't_s,P_generator_W,T_peak_C');
for r = 1:size(S.probes,1)
    fprintf(fid, ',T_probe%d_C_at_%.1f_%.1f_%.1f_mm', r, 1e3*S.probes(r,:));
end
fprintf(fid, '\n');
A = [S.t, S.P_drive, S.T_peak, S.T_probe];
fprintf(fid, [repmat('%.6g,',1,size(A,2)-1) '%.6g\n'], A.');
fclose(fid);
end
