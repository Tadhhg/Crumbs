function [t, S] = time_to_T(P, Tlimit, varargin)
%TIME_TO_T  How long the generator can run before the gel reaches a
%           temperature limit.  The other inverse of T_gel.
%
%   t = TIME_TO_T(P, Tlimit)           seconds until the PEAK temperature in
%                                      the block first reaches Tlimit [degC]
%                                      at CW generator power P [W]
%   t = TIME_TO_T(P, Tlimit, x, z)     limit set at a specified point
%   t = TIME_TO_T(P, Tlimit, x, y, z)
%   t = TIME_TO_T(..., 'delivered')    P is the delivered power
%
%   Returns Inf if the steady state never reaches the limit -- that power can
%   be run indefinitely.  Returns 0 if the limit is already exceeded at t = 0
%   (impossible unless Tlimit < 37 degC).
%
%   P and Tlimit are broadcast against each other.
%
%   METHOD.  Theta(r,t) rises monotonically for a step of power, so the
%   crossing is found by monotone interpolation of the stored history,
%   log-spaced from 1 ms to 600 s.  Below 1 ms the rise is exactly adiabatic,
%   theta = q t /(rho c), and that limit is used.
%
%   See also T_GEL, P_FOR_T.

delivered = {};
args = varargin;
while ~isempty(args) && (ischar(args{end}) || isstring(args{end}))
    if strcmpi(char(args{end}),'delivered'), delivered = {'delivered'};
    else, error('time_to_T:flag','unknown option ''%s''', char(args{end}));
    end
    args(end) = [];
end

M = model_times();
[Pg, Tl] = ndgrid_pair(P, Tlimit);
t = zeros(size(Pg));
for q = 1:numel(Pg)
    th = T_gel(Pg(q), M.t, args{:}, delivered{:}) ;
    th = th(:).';
    target = Tl(q);
    if th(1) >= target, t(q) = 0; continue; end
    i = find(th >= target, 1);
    if isempty(i)
        % never reached inside the history: is it reached at all?
        Tss = T_gel(Pg(q), Inf, args{:}, delivered{:});
        if Tss < target, t(q) = Inf;
        else,            t(q) = NaN;   % between t_end and steady state
        end
        continue
    end
    w = (target - th(i-1))/(th(i) - th(i-1));
    t(q) = M.t(i-1) + w*(M.t(i) - M.t(i-1));
end

if nargout > 1
    S = struct('t',t,'P_generator',Pg,'Tlimit',Tl, ...
               't_history_end',M.t(end));
    S.note = ['NaN means the limit is reached later than the end of the ' ...
              'computed history but before the steady state'];
end
end

function M = model_times()
persistent v
if isempty(v)
    here = fileparts(mfilename('fullpath'));
    f = fullfile(fileparts(here),'results','heat_model.mat');
    w = load(f,'t');  v = struct('t', w.t);
end
M = v;
end

function [A,B] = ndgrid_pair(a,b)
a = a(:); b = b(:).';
if isscalar(a),      A = a*ones(size(b)); B = b;
elseif isscalar(b),  A = a;               B = b*ones(size(a));
else,                [A,B] = ndgrid(a,b);
end
end
