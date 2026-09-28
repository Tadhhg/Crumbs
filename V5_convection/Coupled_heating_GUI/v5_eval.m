function O = v5_eval(R, P)
%V5_EVAL  Temperature on the antenna plane for every stored frame, from the
%         four responses of v5_solve and the instant settings in P.
%
%   O = v5_eval(R, P)
%     O.T      nz x nx x nt temperature (degC, single)
%     O.T_cath 1 x nt catheter temperature (degC); O.Theta = alpha (T_cath - T_b)
%     O.env    the power envelope (v5_envelope: duty cycle, bursts, rests)
%     O.P_del, O.P_c, O.alpha, O.dTc, O.tau_c   the derived instant quantities
%
%   Convection model.docx s3 and s5, with the generator power P_del p(t):
%     theta_d(t_n) = P_del sum_m pbar_m [D(t_n - t_m) - D(t_n - t_m+1)]
%                    (p piecewise constant per frame; CW = the docx's D(t) - D(t - t_on))
%     theta_s      = sum_k [Theta(t_k) - Theta(t_k-1)] U(t - t_k + dt/2)   (Duhamel, midpoint)
%     T = T_b + theta_d + theta_s + (T_room - T_b) Wr + (T_0 - T_b) Z
%   Fast: a few matrix products, no re-solve.

t = R.t;  nt = numel(t);  dt = t(2) - t(1);
[nz, nx, ~] = size(R.D);
O.env = v5_envelope(P, t);
pbar = O.env.pbar;

% ---------------------------------------------------------- POWER -------
P_fwd = max(P.P_set - P.P_refl, 0);
O.P_del = P_fwd*(1 - P.P_dropoff);                 % reaches the antenna while on
O.P_c   = P.gamma*P_fwd*P.P_dropoff;               % heats the catheter tip while on
Dm = reshape(R.D, [], nt);
dD = diff(Dm, 1, 2);                                % D(lag j) - D(lag j-1), j = 1..nt-1
thd = single(O.P_del)*(dD*single(toep(pbar, nt)));

% ------------------------------------------------- CATHETER + SURFACE ---
O.dTc   = O.P_c/P.k_cat;
O.tau_c = P.m_cat*P.c_cat/P.k_cat;
O.T_cath = zeros(1, nt);  O.T_cath(1) = P.T_b;     % lumped tip, exact for p constant per frame
g = exp(-dt/O.tau_c);
for n = 2:nt
    O.T_cath(n) = P.T_b + (O.T_cath(n-1) - P.T_b)*g + O.dTc*pbar(n-1)*(1 - g);
end
O.alpha = min(1, P.beta/(R.P.d*sqrt(P.flow)));     % d is the solved geometry
O.Theta = O.alpha*(O.T_cath - P.T_b);
dTh = diff([0 O.Theta]);                            % dTh(k) = Theta(t_k) - Theta(t_k-1); dTh(1) = 0
Um = reshape(R.U, [], nt);
Ub = 0.5*(Um(:, 1:end-1) + Um(:, 2:end));           % U at lag m*dt + dt/2, m = 0..nt-2
ths = Ub*single(toep(dTh(2:end), nt));

% ------------------------------------------------------- TOTAL ----------
T = P.T_b + thd + ths + (P.T_room - P.T_b)*reshape(R.Wr, [], nt) ...
          + (P.T_0 - P.T_b)*reshape(R.Z, [], nt);        % single: ~4e-6 K resolution
O.T = reshape(T, nz, nx, nt);
O.theta_d = reshape(thd, nz, nx, nt);
O.theta_s = reshape(ths, nz, nx, nt);
O.t = t;
end

function M = toep(a, nt)
% M(j, n) = a(n - j) for 1 <= j <= n-1, else 0   ((nt-1) x nt): out(:, n) = sum_j X(:, j) a(n - j)
M = zeros(nt - 1, nt);
for n = 2:nt
    M(1:n-1, n) = a(n-1:-1:1);
end
end
