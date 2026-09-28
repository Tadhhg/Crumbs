function E = v5_envelope(P, t)
%V5_ENVELOPE  Power envelope p(t) in [0, 1] (fraction of the set power on).
%
%   E = v5_envelope(P, t)      t = the stored frame times (uniform, from 0)
%
%   Treatment runs 0 <= t < t_on.  Within it, power comes in bursts of
%   length t_burst separated by rests of t_rest (t_rest = 0: one continuous
%   burst).  Within each burst the power is pulsed with period T_pulse and
%   duty cycle 'duty' (1 = CW); the pulse train restarts with every burst.
%
%   E.pbar(m)   mean of p over frame interval [t(m), t(m+1)], m = 1..nt-1:
%               the heating model is linear, so pulses shorter than a frame
%               enter through their time average (exact for the energy;
%               within a frame the diffusion length sqrt(a dt) is ~0.5 mm at 2 s)
%   E.p_now(n)  p at the frame time t(n) (instantaneous)
%   E.ts, E.ps  exact stair-step envelope (transition times) for plotting
%   E.mean_on   mean of p over the treatment, E.n_pulses, E.fine = sub-frame pulses?

dt = t(2) - t(1);  nt = numel(t);
duty = min(max(P.duty, 0), 1);
Tp = max(P.T_pulse, 1e-6);
t_on = max(P.t_on, 0);
burst = P.t_rest > 0 && isfinite(P.t_burst) && P.t_burst > 0;

% fine sampling, aligned to the frames: at least 20 samples per frame and
% 20 per pulse high/low phase, capped at 2e6 samples in all
dmin = dt/20;
if duty > 0 && duty < 1, dmin = min(dmin, min(duty, 1 - duty)*Tp/20); end
nsub = max(1, ceil(dt/dmin));
nsub = min(nsub, max(1, floor(2e6/(nt - 1))));
dtf = dt/nsub;
tc = (0.5:(nt - 1)*nsub)*dtf;                     % cell centres
p = pfun(tc);
E.pbar = mean(reshape(p, nsub, nt - 1), 1);
E.p_now = pfun(t);

% stair-step for plotting from the sampled transitions
k = find(diff(p) ~= 0);
te = (k*dtf);                                      % transition times (cell edges)
E.ts = [0, reshape([te; te], 1, []), t(end)];
v = [p(1), p(k + 1)];
E.ps = reshape([v; v], 1, []);
E.mean_on = mean(p(tc < t_on));  if isnan(E.mean_on), E.mean_on = 0; end
E.n_pulses = nnz(diff([0 p]) > 0);
E.fine = duty < 1 && Tp < 2*dt;                    % pulses unresolved by the frames
E.dt = dt;

    function y = pfun(tt)
        tb = tt;
        if burst
            tb = mod(tt, P.t_burst + P.t_rest);
            inb = tb < P.t_burst;
        else
            inb = true(size(tt));
        end
        if duty >= 1
            inp = true(size(tt));
        else
            inp = mod(tb, Tp) < duty*Tp;
        end
        y = double(tt < t_on & inb & inp);
    end
end
