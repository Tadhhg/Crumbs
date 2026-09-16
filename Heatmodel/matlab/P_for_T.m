function [P, S] = P_for_T(Tlimit, texp, varargin)
%P_FOR_T  The generator power at which the gel just reaches a temperature
%         limit after a given exposure.  The inverse of T_gel.
%
%   P = P_FOR_T(Tlimit, texp)          CW generator power [W] at which the
%                                      PEAK temperature in the block reaches
%                                      Tlimit [degC] after texp seconds
%   P = P_FOR_T(Tlimit, Inf)           the maximum power the phantom can take
%                                      indefinitely
%   P = P_FOR_T(Tlimit, texp, x, z)    limit set at a specified point instead
%   P = P_FOR_T(Tlimit, texp, x, y, z)
%   [P, S] = P_FOR_T(...)              plus P delivered, the hot-spot location
%                                      and the resulting thermal dose
%   P = P_FOR_T(..., 'delivered')      return the DELIVERED power
%
%   This is exact, not a search: because
%       T = T_blood + (eta P / P_ref) Theta(r,texp)
%   is linear in P,
%       P = P_ref (Tlimit - T_blood) / (eta Theta(r,texp)) .
%   texp and Tlimit may be vectors; they are broadcast, which is what makes
%   the power-versus-exposure-time design chart a one-liner.
%
%   Example -- the safe-operating chart:
%       te = logspace(-1, 2, 60);
%       plot(te, P_for_T(100, te), te, P_for_T(60, te));
%
%   See also T_GEL, TIME_TO_T.

delivered = false;
args = varargin;
while ~isempty(args) && (ischar(args{end}) || isstring(args{end}))
    if strcmpi(char(args{end}),'delivered'), delivered = true;
    else, error('P_for_T:flag','unknown option ''%s''', char(args{end}));
    end
    args(end) = [];
end

[Tl, te] = ndgrid_pair(Tlimit, texp);
th = zeros(size(te));
for q = 1:numel(te)
    [~,Sq] = T_gel(1, te(q), args{:}, 'delivered');   % Theta per delivered watt
    th(q)  = Sq.theta;
end
% T_gel(1 W delivered) gives theta = Theta/P_ref, so Theta = P_ref*theta
Pdel = (Tl - blood_temp())./th;
P    = Pdel;
if ~delivered
    P = Pdel/eta_used();
end
P(th <= 0) = Inf;

if nargout > 1
    S = struct('P_generator', Pdel/eta_used(), 'P_delivered', Pdel, ...
               'Tlimit', Tl, 't_exposure', te, ...
               'theta_per_delivered_watt', th, 'eta', eta_used());
    S.note = ['theta_per_delivered_watt is the temperature rise in K per ' ...
              'watt actually delivered to the gel'];
end
end

function [A,B] = ndgrid_pair(a,b)
a = a(:); b = b(:).';
if isscalar(a),      A = a*ones(size(b)); B = b;
elseif isscalar(b),  A = a;               B = b*ones(size(a));
else,                [A,B] = ndgrid(a,b);
end
end

function T0 = blood_temp()
persistent v
if isempty(v), [~,S] = T_gel(0, 0); v = S.T_blood; end
T0 = v;
end

function e = eta_used()
persistent v
if isempty(v), [~,S] = T_gel(0, 0); v = S.eta; end
e = v;
end
