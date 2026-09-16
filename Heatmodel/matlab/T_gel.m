function [T, S] = T_gel(P, t, varargin)
%T_GEL  THE DELIVERABLE.  Temperature of the gel as a function of the
%       microwave generator power and time.
%
%   T = T_GEL(P, t)              peak temperature anywhere in the 30 mm cube
%                                below the antenna, at time t   [degC]
%                                -- EXCLUDING a 2 mm sphere around the feed
%                                point, where the prescribed antenna current
%                                of the field model is singular.  See THE
%                                FEED EXCLUSION below; it is not cosmetic.
%   T = T_GEL(P, t, x, z)        temperature at a point on the y = 0 plane
%   T = T_GEL(P, t, x, y, z)     temperature at a point in the block
%   [T, S] = T_GEL(...)          plus a struct: rise above blood, thermal
%                                dose CEM43, where the hot spot is, and every
%                                validity flag that applies
%   T = T_GEL(..., 'delivered')  P is the DELIVERED power, not the generator
%                                forward power (see UNITS below)
%
%   P may also be a FUNCTION HANDLE @(t) giving the generator power at time t,
%   for pulsed or ramped drive; see PULSED OPERATION below.
%
%   UNITS.  x, y, z in metres, t in seconds, P in watts, T in degrees Celsius.
%   z is the height, so z is NEGATIVE inside the gel (surface z = 0, floor
%   z = -30 mm), exactly as in E_gel.m.  x runs along the antenna, feed at
%   x = 0, tip at x = 12 mm.
%
%   WHICH POWER?  P is by default the GENERATOR FORWARD power, because that is
%   the number on the front panel.  The gel only sees the DELIVERED power
%       P_delivered = eta * P_generator,   eta = 0.75 baseline,
%   i.e. 200 W forward minus the 50 W maximum reflected quoted for the
%   generator.  eta is the single most important number to check against your
%   own measurement: if your match is better than the worst case, the gel gets
%   more power and heats proportionally faster.  Pass 'delivered' to bypass
%   eta and give the delivered power directly.
%
%   THE FUNCTION
%   The heat equation with a fixed spatial source shape is LINEAR, and the
%   deposition is proportional to the delivered power, so
%
%       T(r,t) = T_blood + (P_delivered / P_ref) * Theta(r,t)
%
%   EXACTLY, where Theta is the temperature rise computed once by
%   RUN_HEAT_MODEL at P_ref = 150 W delivered, and T_blood = 37 degC.  There
%   is no approximation in the power scaling itself; its limits are the
%   constant-property assumption (see VALIDITY).
%
%   Note the contrast with the field model: |E| scales as sqrt(P), because the
%   field is an amplitude; temperature rise scales as P, because heating goes
%   with the square of the field.
%
%   PULSED OPERATION.  Because the system is linear and time invariant, the
%   response to ANY generator power history P(t) is a Duhamel superposition of
%   the step response Theta:
%
%       theta(r,t) = (1/P_ref) Integral_0^t P_del(tau) dTheta/dt (r, t-tau) dtau
%
%   Pass a function handle for P and this is evaluated by discrete
%   superposition of step responses, which is exact for a piecewise-constant
%   power history.  Example, 60 s at 40 W with a 30 %% duty cycle at 1 Hz:
%       Pw = @(tt) 40*(mod(tt,1) < 0.3);
%       T  = T_gel(Pw, linspace(0,60,2000));
%
%   t = Inf is accepted and returns the steady state.
%
%   THE FEED EXCLUSION
%   The field model prescribes the antenna current as sin(k_e (L-x)) on
%   0 <= x <= L.  That current is at its MAXIMUM where it is truncated, at the
%   feed, which puts a fictitious point charge there and makes the deposition
%   -- and hence the temperature -- singular at (0,0,0).  The peak returned by
%   T_gel(P,t) therefore excludes a 2 mm sphere about the feed.  The excluded
%   region is genuinely not predictable with this source model: the peak
%   inside it does not converge under grid refinement, whereas everything more
%   than about 3 mm away is unaffected to better than 1 %%.  A real feed
%   junction IS a hot spot, so watch it on the bench; this model just cannot
%   tell you how hot.  Full diagnosis in feed_regularise.m and the write-up,
%   Sec. 6.2.  Point queries (x,z) and (x,y,z) are NOT masked -- ask for a
%   point inside the sphere and you get the singular value.
%
%   VALIDITY -- read this before quoting a number
%     1. Above about 100 degC the gel boils.  Nothing in this model knows
%        that; it will happily report 900 degC.  S.flags tells you when the
%        result has left the physical range.
%     2. Above about 60 degC the permittivity and conductivity of the medium
%        drift by order 10 %%, so the deposition is no longer exactly
%        proportional to power.  See the write-up, Sec. 9F.
%     3. Low-acyl gellan gels soften and melt well below boiling.  The phantom
%        may stop being a phantom before the model stops being arithmetic.
%     4. The surface is held at 37 degC (perfect heat sink).  That is a LOWER
%        BOUND on temperature; a finite film coefficient gives more.  The
%        convective case is quantified in the write-up, Sec. 9B.
%     5. Everything inherits the field model assumptions, above all the 0.5 mm
%        standoff between the conductor and the gel surface.
%
%   See also P_FOR_T, TIME_TO_T, CEM43, RUN_HEAT_MODEL, E_GEL.

persistent M
if isempty(M)
    here = fileparts(mfilename('fullpath'));
    f = fullfile(fileparts(here),'results','heat_model.mat');
    if ~exist(f,'file')
        error('T_gel:noModel', ...
            'Heatmodel/results/heat_model.mat not found.  Run  run_heat_model  first.');
    end
    M = load(f);
    M.s   = -M.z(:).';                    % depth, ascending, for interpolation
    M.s3  = -M.z3(:).';
    M.Th_y0 = double(M.Th_y0);            % [nx x nz x Nt]
    M.Th3   = double(M.Th3);
    M.Th_ss_y0 = double(M.Th_ss_y0);
    M.Th_ss = double(M.Th_ss);
end

% ---- trailing flags ------------------------------------------------------
args = varargin;
delivered = false;
while ~isempty(args) && (ischar(args{end}) || isstring(args{end}))
    switch lower(char(args{end}))
        case 'delivered', delivered = true;
        otherwise, error('T_gel:flag','unknown option ''%s''', char(args{end}));
    end
    args(end) = [];
end
nc = numel(args);
if ~ismember(nc, [0 2 3])
    error('T_gel:nargin', ...
        'give no coordinates (peak), (x,z) on y = 0, or (x,y,z)');
end

t   = double(t);
sz  = size(t);
tv  = t(:).';
inf_mask = isinf(tv) & tv > 0;
tf  = tv;  tf(inf_mask) = M.t(end);
if any(tf > M.t(end)+1e-9)
    warning('T_gel:tmax', ...
        ['t beyond the computed history (%.0f s) is clamped to it. ' ...
         'Use t = Inf for the steady state.'], M.t(end));
    tf = min(tf, M.t(end));
end
if any(tf < 0), error('T_gel:tneg','t must be >= 0'); end

% ---- the kernel: Theta(t) for the requested location(s) ------------------
switch nc
    case 0                                   % peak over the block
        [kern, kss, loc] = peak_kernels(M, isa(P,'function_handle'));
    case 2
        xq = args{1};  zq = args{2};
        [kern, kss] = plane_kernel(M, xq, zq);
        loc = struct('x',xq,'y',0,'z',zq);
    case 3
        xq = args{1};  yq = args{2};  zq = args{3};
        [kern, kss] = block_kernel(M, xq, yq, zq);
        loc = struct('x',xq,'y',yq,'z',zq);
end
% kern : [ncand x Nt] rows are Theta(t) at candidate points
% kss  : [ncand x 1]  their steady-state values

% ---- the power history ---------------------------------------------------
if isa(P,'function_handle')
    theta_rows = duhamel(kern, M.t, P, tf, delivered, M.eta, M.Pdel_ref);
    Pdel_report = NaN;  scale = NaN;
else
    if ~isscalar(P), error('T_gel:P','P must be a scalar or a function handle'); end
    Pdel_report = P;  if ~delivered, Pdel_report = M.eta*P; end
    scale = Pdel_report/M.Pdel_ref;
    theta_rows = scale*interp1(M.t, kern.', tf, 'linear').';
    if size(kern,1) == 1, theta_rows = reshape(theta_rows, 1, []); end
    if any(inf_mask)
        theta_rows(:,inf_mask) = repmat(scale*kss, 1, nnz(inf_mask));
    end
end

theta = max(theta_rows, [], 1);
if nc == 0
    [~,ib] = max(theta_rows, [], 1);
else
    ib = ones(size(theta));
end

T = reshape(M.T_blood + theta, sz);

if nargout > 1
    S = struct();
    S.theta   = reshape(theta, sz);
    S.T       = T;
    S.T_blood = M.T_blood;
    S.Pgen      = P;
    S.Pdelivered= Pdel_report;
    S.eta       = M.eta;
    S.scale     = scale;
    S.Pdel_ref  = M.Pdel_ref;
    S.t         = t;
    if nc == 0
        if numel(loc) == numel(M.t)
            % constant power: the stored per-time hot-spot location applies
            [~,itl] = min(abs(M.t(:) - tf(:).'), [], 1);
            S.location = loc(itl);
        else
            S.location = loc(ib);
        end
        S.note_location = 'hot spot, (x,y,z) in m, feed sphere excluded';
    else
        S.location = loc;
    end
    if numel(tv) > 1
        [tsorted, isort] = sort(tv);
        S.CEM43 = cem43(tsorted, M.T_blood + theta(isort));
    else
        S.CEM43 = NaN;
    end
    S.flags = {};
    if max(T(:)) > 100
        S.flags{end+1} = 'T exceeds 100 degC: the gel boils, the model does not';
    end
    if max(T(:)) > 60
        S.flags{end+1} = 'T exceeds 60 degC: permittivity/conductivity drift ~10 %, power scaling no longer exact';
    end
    if max(T(:)) > 50
        S.flags{end+1} = 'T exceeds 50 degC: check the gel-sol transition of the Phytagel itself';
    end
    S.flags{end+1} = 'surface pinned at 37 degC: this is a LOWER BOUND on temperature';
end
end

% =========================================================================
function [kern, kss, loc] = peak_kernels(M, waveform)
%PEAK_KERNELS  Theta(t) for the hottest point in the block, feed sphere
%              excluded.
%
%   CONSTANT POWER (the usual case).  RUN_HEAT_MODEL already searched the
%   whole 3-D block at every output time and stored the masked maximum in
%   M.Th_max, with its location in M.ij_max.  Since theta scales exactly with
%   power, scaling that curve IS the peak -- no candidate set, no
%   approximation.  This branch is exact.
%
%   ARBITRARY POWER HISTORY.  Duhamel superposition has to be applied at a
%   FIXED point, so the maximum is taken over the set of locations that are
%   ever the hot spot.  Those are read off M.ij_max and restricted to the
%   symmetry plane y = 0, where the full-resolution history is stored.  The
%   hot spot does leave that plane briefly (around 2 s at the reference
%   power), but only because the feed-exclusion sphere cuts across it, and by
%   about one millimetre; the resulting error in a waveform peak is well under
%   a percent.  The candidate set is filtered by M.pmask_y0 so that no point
%   inside the exclusion sphere can be selected.
if nargin < 2, waveform = false; end

if ~waveform
    kern = M.Th_max(:).';
    kss  = M.Th_max_ss;
    loc  = struct('x',{},'y',{},'z',{});
    ij = M.ij_max;
    for r = 1:numel(M.t)
        if ij(1,r) == 0
            loc(r) = struct('x',NaN,'y',NaN,'z',NaN);
        else
            loc(r) = struct('x',M.x(ij(1,r)),'y',M.y(ij(2,r)),'z',M.z(ij(3,r)));
        end
    end
    return
end

ij = M.ij_max(:,2:end);
on = ij(2,:) == M.iy0;
pairs = unique([ij([1 3],on).'; M.ij_ss([1 3]).'], 'rows');
ok = false(size(pairs,1),1);
for r = 1:size(pairs,1)
    ok(r) = M.pmask_y0(pairs(r,1), pairs(r,2));
end
pairs = pairs(ok,:);
if isempty(pairs)                    % should not happen; fall back to the curve
    kern = M.Th_max(:).';  kss = M.Th_max_ss;
    loc  = struct('x',NaN,'y',NaN,'z',NaN);
    return
end
n = size(pairs,1);
kern = zeros(n, numel(M.t));
kss  = zeros(n,1);
loc  = struct('x',cell(1,n),'y',cell(1,n),'z',cell(1,n));
for r = 1:n
    kern(r,:) = squeeze(M.Th_y0(pairs(r,1), pairs(r,2), :)).';
    kss(r)    = M.Th_ss_y0(pairs(r,1), pairs(r,2));
    loc(r) = struct('x',M.x(pairs(r,1)),'y',0,'z',M.z(pairs(r,2)));
end
end

function [kern, kss] = plane_kernel(M, xq, zq)
%PLANE_KERNEL  Theta(t) at (xq, 0, zq), full-resolution y = 0 plane.
check_in(M, xq, 0, zq);
w = bilinear_weights(M.x, M.s, xq, -zq);
kern = zeros(1, numel(M.t));
for c = 1:4
    kern = kern + w.c(c)*squeeze(M.Th_y0(w.i(c), w.j(c), :)).';
end
kss = 0;
for c = 1:4
    kss = kss + w.c(c)*M.Th_ss_y0(w.i(c), w.j(c));
end
end

function [kern, kss] = block_kernel(M, xq, yq, zq)
%BLOCK_KERNEL  Theta(t) at (xq,yq,zq) from the 3-D history.
%   The 3-D history is stored decimated (0.5 mm) to keep the file usable; the
%   y = 0 plane is kept at the full 0.25 mm.  For points near the surface,
%   where the gradient is steepest, prefer the (x,z) form.
check_in(M, xq, yq, zq);
nt3  = numel(M.t3);
k3   = zeros(1, nt3);
for m = 1:nt3
    k3(m) = interpn(M.x3(:).', M.y3(:).', M.s3, M.Th3(:,:,:,m), ...
                    xq, yq, -zq, 'linear');
end
kern = interp1(M.t3, k3, M.t, 'linear', 'extrap');
% the steady state is stored at full resolution, so use it directly
kss  = interpn(M.x(:).', M.y(:).', -M.z(:).', M.Th_ss, xq, yq, -zq, 'linear');
end

function check_in(M, xq, yq, zq)
if any(xq(:) < M.x(1)-1e-9) || any(xq(:) > M.x(end)+1e-9) || ...
   any(yq(:) < M.y(1)-1e-9) || any(yq(:) > M.y(end)+1e-9) || ...
   any(zq(:) > 1e-9)        || any(zq(:) < M.z(end)-1e-9)
    error('T_gel:outside', ...
        ['point outside the region of interest: x %.1f..%.1f mm, ' ...
         'y %.1f..%.1f mm, z %.1f..0 mm'], ...
        1e3*M.x(1), 1e3*M.x(end), 1e3*M.y(1), 1e3*M.y(end), 1e3*M.z(end));
end
end

function w = bilinear_weights(xg, sg, xq, sq)
xg = xg(:).'; sg = sg(:).';
i0 = max(1, min(numel(xg)-1, find(xg <= xq, 1, 'last')));
j0 = max(1, min(numel(sg)-1, find(sg <= sq, 1, 'last')));
if isempty(i0), i0 = 1; end
if isempty(j0), j0 = 1; end
a = (xq - xg(i0))/(xg(i0+1)-xg(i0));
b = (sq - sg(j0))/(sg(j0+1)-sg(j0));
w.i = [i0 i0+1 i0   i0+1];
w.j = [j0 j0   j0+1 j0+1];
w.c = [(1-a)*(1-b), a*(1-b), (1-a)*b, a*b];
end

function rows = duhamel(kern, tk, Pfun, tq, delivered, eta, Pref)
%DUHAMEL  theta(t) for an arbitrary generator power history.
%   Superposition of step responses, exact for piecewise-constant power:
%       theta(t_n) = sum_m  Pbar_m/Pref * [ Theta(t_n-t_{m-1}) - Theta(t_n-t_m) ]
%   which is the discrete convolution of the mean power in each interval with
%   the increments of the step response.
tmax = max(tq(:));
if tmax <= 0, rows = zeros(size(kern,1), numel(tq)); return; end
nu   = 4000;
dtu  = tmax/nu;
tu   = (0:nu)*dtu;
Pm   = 0.5*(Pfun(tu(1:end-1)) + Pfun(tu(2:end)));    % mean power per interval
if ~delivered, Pm = eta*Pm; end
p    = Pm(:).'/Pref;
rows = zeros(size(kern,1), numel(tq));
for r = 1:size(kern,1)
    g  = interp1(tk, kern(r,:), tu, 'linear');
    g(1) = 0;
    dg = diff(g);
    y  = conv(p, dg);
    y  = [0 y(1:nu)];                                 % theta at tu
    rows(r,:) = interp1(tu, y, min(tq,tmax), 'linear');
end
end
