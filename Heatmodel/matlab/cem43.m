function [D, Dcum] = cem43(t, T)
%CEM43  Cumulative equivalent minutes at 43 degC (thermal dose).
%
%   D            = CEM43(t, T)   total dose over the whole record [min]
%   [D, Dcum]    = CEM43(t, T)   also the running cumulative dose
%
%   t  [s]   time samples, increasing.  Row or column.
%   T  [degC] temperature at those times.  If T is a matrix, time must be the
%            LAST dimension and the dose is computed independently for every
%            other element (so a whole map can be dosed in one call).
%
%   THE DEFINITION
%       CEM43 = Integral R^(43 - T(t)) dt ,  t in MINUTES
%       R = 0.5 for T >= 43 degC,   R = 0.25 for T < 43 degC
%   [S. A. Sapareto and W. C. Dewey, "Thermal dose determination in cancer
%    therapy", Int. J. Radiat. Oncol. Biol. Phys. 10(6):787-800, 1984,
%    PMID 6547421.]  R is the reciprocal of the factor by which the time to a
%   given biological effect changes per degree; Sapareto and Dewey note R may
%   lie between 0.4 and 0.8 above 43 degC, but 0.5/0.25 is the universally
%   used convention.
%
%   THRESHOLDS USED IN THIS WORK
%     CEM43 = 240 min   the conventional threshold for thermal ablation /
%                       complete necrosis in muscle and most soft tissue
%     CEM43 = 60 min    a commonly quoted lower bound for irreversible damage
%   These are biological thresholds.  A Phytagel phantom has no biology, so in
%   this model they are reported as a TRANSLATION of the phantom result into
%   the tissue equivalent it stands in for -- see the write-up, Sec. 8.4.
%
%   NUMERICAL NOTE
%   The integrand spans many orders of magnitude (at 100 degC it is 1.4e17
%   min/min), so a plain trapezoid on the integrand is useless.  Each interval
%   is integrated EXACTLY under the assumption that T varies linearly across
%   it, which makes the integrand a pure exponential in t:
%       Int_{t1}^{t2} f dt = (f2 - f1) (t2 - t1) / ln(f2/f1) ,   f = R^(43-T)
%   reducing to f*(t2-t1) when f2 == f1.

t = t(:).';
sz = size(T);
if isvector(T), T = T(:).'; sz = size(T); end
nt = sz(end);
if numel(t) ~= nt
    error('cem43:size','numel(t) = %d must equal size(T,end) = %d', numel(t), nt);
end
Tm = reshape(T, [], nt);

R = 0.25 + 0.25*(Tm >= 43);                 % 0.25 below 43 degC, 0.5 above
lf = (43 - Tm).*log(R);                     % ln of the integrand
lf = min(lf, 700);                          % guard the exp against overflow
f  = exp(lf);

dt_min = diff(t)/60;                        % minutes
f1 = f(:,1:end-1);  f2 = f(:,2:end);
r  = f2./f1;
inc = f1.*repmat(dt_min, size(f1,1), 1);    % limit for f2 == f1
lr  = log(r);
use = abs(lr) > 1e-8 & isfinite(lr);
DTm = repmat(dt_min, size(f1,1), 1);
inc(use) = (f2(use)-f1(use)).*DTm(use)./lr(use);

Dcum = [zeros(size(f1,1),1), cumsum(inc,2)];
D    = Dcum(:,end);

if numel(sz) > 2
    D    = reshape(D,    sz(1:end-1));
    Dcum = reshape(Dcum, sz);
else
    D    = reshape(D, [sz(1) 1]);
    Dcum = reshape(Dcum, sz);
end
if isscalar(D), Dcum = Dcum(:).'; end
end
