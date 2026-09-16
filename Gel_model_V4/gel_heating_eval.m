function [T, M] = gel_heating_eval(R, P_set, T_w, t_pre, t, T_protect, T_lesion)
%GEL_HEATING_EVAL  Temperature on the antenna plane for one scenario, and the
%                  surface-sparing metrics.
%
%   [T, M] = GEL_HEATING_EVAL(R, P_set, T_w, t_pre, t)
%     R      solution struct from gel_heating_solve
%     P_set  [W]    generator setting (delivered = R.eta * P_set)
%     T_w    [degC] water temperature from t = -t_pre onwards (the gel starts
%                   in equilibrium with water at R.T_bath0)
%     t_pre  [s]    pre-cooling time before the power is switched on
%     t      [s]    time since the power was switched on
%   T is [nx x nz] degC (z increases toward the surface, which is z = 0).
%   M: D        [m]  protected depth: the layer below the surface in which
%                    T <= T_protect everywhere (min over x of the first
%                    crossing, interpolated); Inf if nothing exceeds T_protect
%      T_les    max T at depth >= D  (the lesion), NaN if D is Inf
%      T_peak, ix_peak, iz_peak
%      T_top3, T_below3   max T within / below the first 3 mm
%      T_cath   catheter temperature (lumped cross-flow model)
%      D_les    depth at which T first exceeds T_lesion (min over x), Inf if none

if nargin < 6, T_protect = 45; end
if nargin < 7, T_lesion  = 50; end

T = double(R.T_ss) + (T_w - R.T_bath0)*frame(R.U, t_pre + t, R) ...
                   + (P_set/R.P_ref)*frame(R.Theta, t, R);
if nargout < 2, return; end

depth = -R.z(:).';                              % [1 x nz], 0 at the surface
M = struct();
[M.T_peak, i] = max(T(:));  [M.ix_peak, M.iz_peak] = ind2sub(size(T), i);
M.D     = crossing_depth(T, depth, T_protect);
M.D_les = crossing_depth(T, depth, T_lesion);
if isfinite(M.D)
    M.T_les = max(T(:, depth >= M.D - 1e-12), [], 'all');
else
    M.T_les = NaN;
end
M.T_top3   = max(T(:, depth <= 3e-3 + 1e-12), [], 'all');
M.T_below3 = max(T(:, depth >  3e-3 + 1e-12), [], 'all');
sel = R.x >= 0 & R.x <= R.L_el;
T_s = mean(T(sel, end));                        % gel surface under the element
c = R.cath;
M.T_cath = T_w + ((P_set/R.P_ref)*c.P_film + c.G_gap*(T_s - T_w))/(c.h_c*c.A_c + c.G_gap);
end

function F = frame(A, tq, R)
% linear interpolation between stored frames
dtf = R.t(2) - R.t(1);
fi  = min(max(tq, 0), R.t(end))/dtf + 1;
i0  = floor(fi);  w = fi - i0;  i1 = min(i0 + 1, numel(R.t));
F   = (1 - w)*double(A(:, :, i0)) + w*double(A(:, :, i1));
end

function D = crossing_depth(T, depth, Tc)
% for each column the depth at which T first exceeds Tc going down from the
% surface (interpolated between nodes); D = min over columns.  Inf if never.
nz = numel(depth);
D = Inf;
above = T > Tc;                                  % [nx x nz]
cols = find(any(above, 2)).';
for ix = cols
    k = find(above(ix, :), 1, 'last');           % nearest-to-surface node above Tc
    if k == nz                                   % exceeded at the surface node itself
        d = 0;
    else
        % node k is above Tc, node k+1 (closer to the surface) is not
        f = (Tc - T(ix, k+1))/(T(ix, k) - T(ix, k+1));
        d = depth(k+1) + f*(depth(k) - depth(k+1));
    end
    D = min(D, d);
end
end
