function v5_selftest()
%V5_SELFTEST  Checks of the V5 solver against closed forms.  Prints PASS/FAIL.
%   1. E-field table = ../Insulated_monopole_Erms/insulated_monopole_Erms.m (photo preset, 15 W)
%   2. U with the surface pinned and a wide footprint = erfc(|z|/(2 sqrt(a t)))  (docx s4)
%   3. D at early time = Q t/(rho c)   (adiabatic rise, docx s1)
%   4. Z with every face adiabatic stays at 1 K; Wr with h_w only = 0 at t = 0
%   5. v5_eval: linearity in P_del and the switch-off superposition
%   6. 2-D and 3-D agree for a y-uniform problem (U with pinned surface)

ok = true;

% ---- 1 -------------------------------------------------------------------
P = v5_defaults('photo');
x = -5e-3:0.1e-3:P.L + 5e-3;  depth = [0 1 2 3 5 10]*1e-3;
E = v5_erms(P, x, P.d + depth)*sqrt(50*(1 - 0.7));     % 15 W, as the original script
[pk, ix] = max(E, [], 2);
% printed by the original script (source segments 0.02 mm; here 0.05 mm)
ref_pk = [5907 2341 1513 1217 872 439];  ref_x = [-0.8 -0.7 0.1 2.9 3.4 4.0]*1e-3;
fprintf('1. E-field (photo antenna, 15 W):\n');
for i = 1:numel(depth)
    fprintf('   depth %4.1f mm: peak %5.0f V/m at x = %4.1f mm   (original %5.0f at %4.1f)\n', ...
            depth(i)*1e3, pk(i), x(ix(i))*1e3, ref_pk(i), ref_x(i)*1e3);
end
ok = check(ok, 'E-field = insulated_monopole_Erms.m (peaks within 0.5 %, same x)', ...
           max(abs(pk(:).'./ref_pk - 1)) < 5e-3 && max(abs(x(ix) - ref_x)) < 0.051e-3);

% ---- 2 -------------------------------------------------------------------
P = v5_defaults('table');
P.h_s = Inf;  P.h_w = 0;  P.L_c = 1;  P.dh = 0.25e-3;  P.t_end = 120;  P.dt_frame = 10;
P.gel_H = 20e-3;  P.dim = 2;
R = v5_solve(P);
kap = P.k/(P.rho*P.c);
[~, ixm] = min(abs(R.x - 13.5e-3));
err = 0;
for n = [2 7 13]
    th = erfc(abs(R.z)/(2*sqrt(kap*R.t(n))));
    err = max(err, max(abs(double(R.U(:, ixm, n)).' - th)));
end
ok = check(ok, sprintf('U = erfc (max error %.4f K per K)', err), err < 0.01);

% ---- 3 -------------------------------------------------------------------
% early time: T = t/(rho c) (Q + kap t/2 lap Q) + O(t^3); lap Q from the grid
n = 2;                                   % t = 10 s: diffusion length 1.2 mm
Dnum = double(R.D(:, :, n));  dh = R.x(2) - R.x(1);
lapQ = 4*del2(R.Q, dh);                  % del2 = lap/4 in 2-D
Dad = R.t(n)/(P.rho*P.c)*(R.Q + kap*R.t(n)/2*lapQ);
[~, iz] = min(abs(R.z + 6e-3));          % 6 mm deep: away from the pinned surface
cols = R.x > 0 & R.x < P.L;
rel = max(abs(Dnum(iz, cols) - Dad(iz, cols)))/max(Dad(iz, cols));
ok = check(ok, sprintf('D = (Q + a t/2 lap Q) t/(rho c) at 6 mm, t = 10 s (rel. error %.4f)', rel), rel < 0.01);

% ---- 4 -------------------------------------------------------------------
P4 = P;  P4.h_s = 0;  P4.h_w = 0;  P4.dh = 1e-3;  P4.dim = 3;
R4 = v5_solve(P4);
dev = max(abs(double(R4.Z(:)) - 1));
ok = check(ok, sprintf('Z stays at 1 K when adiabatic (max dev %.2e)', dev), dev < 1e-6);
P4.h_w = 8;  R4 = v5_solve(P4);
ok = check(ok, 'Wr = 0 at t = 0', max(abs(double(R4.Wr(:,:,1))), [], 'all') < 1e-9);
ok = check(ok, 'Wr in [0, 1]', min(R4.Wr(:)) > -1e-6 && max(R4.Wr(:)) < 1 + 1e-6);

% ---- 5 -------------------------------------------------------------------
P5 = v5_defaults('table');  P5.dh = 1e-3;  P5.dim = 3;
R5 = v5_solve(P5);
P5.gamma = 0;  P5.T_0 = P5.T_b;  P5.T_room = P5.T_b;
A = v5_eval(R5, P5);  P5b = P5;  P5b.P_set = 2*P5.P_set;  B = v5_eval(R5, P5b);
lin = max(abs((B.T(:) - P5.T_b) - 2*(A.T(:) - P5.T_b)));
ok = check(ok, sprintf('linear in power (max dev %.2e K)', lin), lin < 1e-4);                 % single-precision storage
P5c = P5;  P5c.t_on = 1e9;  C = v5_eval(R5, P5c);
n_on = find(R5.t <= P5.t_on, 1, 'last');
dd = max(abs(A.T(:,:,1:n_on) - C.T(:,:,1:n_on)), [], 'all');
ok = check(ok, 'switch-off changes nothing before t_on', dd < 1e-9);
ok = check(ok, 'gel cools after switch-off', max(A.T(:,:,end), [], 'all') < max(A.T(:,:,n_on), [], 'all'));

% ---- 6 -------------------------------------------------------------------
P6 = P;  P6.dh = 0.5e-3;  P6.dim = 3;  P6.gel_W = 20e-3;
R6 = v5_solve(P6);  P6.dim = 2;  R6b = v5_solve(P6);
d6 = max(abs(double(R6.U(:)) - double(R6b.U(:))));
ok = check(ok, sprintf('2-D = 3-D for a y-uniform surface step (max dev %.2e)', d6), d6 < 1e-6);

% ---- 7 -------------------------------------------------------------------
% power envelope: fast pulses at 50 % duty heat like CW at half power (linear,
% period << diffusion time); rests are off; catheter = the analytic lumped solution
P7 = P5;  P7.t_on = 400;
A = v5_eval(R5, P7);
P7b = P7;  P7b.duty = 0.5;  P7b.T_pulse = 1;  B = v5_eval(R5, P7b);
n300 = find(R5.t == 300);
rA = A.T(:,:,n300) - P7.T_b;  rB = B.T(:,:,n300) - P7.T_b;
dev = max(abs(rB(:) - 0.5*rA(:)))/max(rA(:));
ok = check(ok, sprintf('1 s pulses at 50 %% duty = half-power CW at 300 s (dev %.4f of the peak rise)', dev), dev < 0.005);
P7c = P7;  P7c.t_burst = 30;  P7c.t_rest = 20;  Cc = v5_eval(R5, P7c);
tm = R5.t(1:end-1) + 1;                         % frame midpoints
inRest = mod(tm, 50) >= 30 & tm < P7c.t_on;
ok = check(ok, 'bursts: no power during the rests, full power in the bursts', ...
           all(Cc.env.pbar(inRest) == 0) && all(Cc.env.pbar(~inRest & tm < 400 & mod(tm, 50) < 30) == 1));
ok = check(ok, 'mean power over the treatment = 30/50', abs(Cc.env.mean_on - 0.6) < 1e-3);
P7d = P7;  P7d.gamma = 0.043;  Dd = v5_eval(R5, P7d);
Tc = P7d.T_b + Dd.dTc*(1 - exp(-R5.t/Dd.tau_c));
off = R5.t > P7d.t_on;  Ton = P7d.T_b + Dd.dTc*(1 - exp(-P7d.t_on/Dd.tau_c));
Tc(off) = P7d.T_b + (Ton - P7d.T_b)*exp(-(R5.t(off) - P7d.t_on)/Dd.tau_c);
ok = check(ok, sprintf('catheter tip = analytic lumped solution (max dev %.1e K)', max(abs(Dd.T_cath - Tc))), ...
           max(abs(Dd.T_cath - Tc)) < 1e-9);

if ok, fprintf('\nALL PASS\n'); else, fprintf('\nSOME CHECKS FAILED\n'); end
end

function ok = check(ok, msg, cond)
if cond, fprintf('PASS  %s\n', msg); else, fprintf('FAIL  %s\n', msg); ok = false; end
end
