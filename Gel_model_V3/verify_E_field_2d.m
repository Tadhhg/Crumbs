%VERIFY_E_FIELD_2D  Checks on E_field_2d.  Run from this folder.
%   1  the wire field equals the sum of ../matlab/hertzian_dipole_field.m
%      (the kernel that validated the layered model) over the same current
%   2  energy conservation: the outward Poynting flux through a closed
%      cylinder round the wire equals the power dissipated outside it
%   3  the power integral that sets the scale is converged
%   4  the sampled waveform reproduces the closed-form |E(t)| envelope
%   5  the numbers quoted in E_field_2d_notes.md

addpath('../matlab');
P_set = 100;                                     % 20 W delivered
mu0 = 1.25663706212e-6;

% ---- 1: kernel ------------------------------------------------------------
xp = [3 10 20 -5]*1e-3;   zp = [-2 -5 -10 -1]*1e-3;
[~, S] = E_field_2d(xp, zp, 0, P_set);
Ex = zeros(size(xp));  Ez = Ex;
for s = 1:numel(S.xs)
    [ex, ~, ez] = hertzian_dipole_field(S.params.f, S.params.epsr, S.Im*S.Idl(s), ...
                                        xp - S.xs(s), 0*xp, zp - S.params.h);
    Ex = Ex + ex;  Ez = Ez + ez;
end
e1 = max(abs([Ex - S.Ex, Ez - S.Ez]))/max(abs([S.Ex S.Ez]));
fprintf('1  kernel vs hertzian_dipole_field: rel err %.1e\n', e1);

% ---- 2: Poynting flux through a closed cylinder vs dissipation outside ----
% H of an x-directed element at (xs,0,h) seen at (x,0,z):  H = Idl grad G x xhat
% = -Idl (jk + 1/R) G (Z/R) yhat  [Balanis AT3 Eq. (4-8a), rotated to xhat].
k = S.kc;  w = 2*pi*S.params.f;  h = S.params.h;
EH = @(x, z) local_EH(x, z, S.xs, S.Idl, k, w, h);
rc = 2e-3;  x0 = -S.params.L_sheath - 2e-3;  x1 = S.params.L + 2e-3;
xs_ = linspace(x0, x1, 2001);                              % side, rho = rc
[ex, ez, hy] = EH(xs_, h - rc);
P_side = trapz(xs_, -0.5*real(ex.*conj(hy)))*2*pi*rc;      % n = -zhat = rho_hat
rr = linspace(1e-5, rc, 801);                               % end discs
[ex, ez, hy] = EH(x1 + 0*rr, h - rr);  P_cap1 = trapz(rr, -0.5*real(ez.*conj(hy)).*rr)*2*pi;
[ex, ez, hy] = EH(x0 + 0*rr, h - rr);  P_cap0 = trapz(rr,  0.5*real(ez.*conj(hy)).*rr)*2*pi;
P_flux = P_side + P_cap1 + P_cap0;
% dissipation outside the cylinder, for I_m = 1
pad = 3/S.alpha;  xg = (x0 - pad):0.25e-3:(x1 + pad);
rg = [logspace(-5, log10(rc), 120), linspace(rc*1.01, rc + pad, 300)];
[XG, RG] = ndgrid(xg, rg);
[ex, ez] = EH(XG, h - RG);
qv = 0.5*S.sig_eff*(abs(ex).^2 + abs(ez).^2);
qv(XG >= x0 & XG <= x1 & RG < rc) = 0;                      % inside the cylinder
P_vol = trapz(xg, trapz(rg, qv.*RG*2*pi, 2), 1);
fprintf('2  Poynting flux out of cylinder %.4f W, dissipation outside %.4f W, diff %.2f %%\n', ...
        P_flux, P_vol, 100*(P_flux - P_vol)/P_vol);

% ---- 3: convergence of the scale --------------------------------------------
[~, Sa] = E_field_2d(5e-3, -5e-3, 0, P_set);
[~, Sb] = E_field_2d(5e-3, -5e-3, 0, P_set, 'ds', 0.05e-3);
[~, Sc] = E_field_2d(5e-3, -5e-3, 0, P_set, 'dxg', 0.125e-3, 'Nr', 480);
[~, Sd] = E_field_2d(5e-3, -5e-3, 0, P_set, 'npad', 5);
fprintf('3  P1 baseline %.4f W;  ds/2 %+.2f %%;  grid x2 %+.2f %%;  extent 5/alpha %+.2f %%\n', ...
        Sa.P1, 100*(Sb.P1/Sa.P1-1), 100*(Sc.P1/Sa.P1-1), 100*(Sd.P1/Sa.P1-1));

% ---- 4: waveform vs closed-form envelope -----------------------------------
ph = linspace(0, 2*pi, 721);
[Et, S4] = E_field_2d(5e-3, -5e-3, ph, P_set);
Emag = hypot(Et.x, Et.z);
Eform = sqrt(0.5*(2*S4.E_rms^2 + abs(S4.Q)*cos(2*ph + angle(S4.Q))));
fprintf('4  |E(t)| sampled vs closed form: max rel err %.1e;  E_peak %.0f, E_min %.0f V/m\n', ...
        max(abs(Emag - Eform))/max(Emag), S4.E_peak, S4.E_min);

% ---- 5: the numbers in the notes -------------------------------------------
fprintf('5  kc = %.1f %+.1fj 1/m, lambda = %.2f mm, 1/alpha = %.2f mm, sig_eff = %.3f S/m\n', ...
        real(S.kc), imag(S.kc), S.lambda*1e3, 1e3/S.alpha, S.sig_eff);
fprintf('   I_m = %.3f A, |I_feed| = %.3f A at P_del = %.0f W\n', S.Im, abs(S.I_feed), S.P_del);
pts = [5 -2; 5 -5; 5 -10; 13.5 -5; 25 -5; -5 -5]*1e-3;
fprintf('   %7s %7s %9s %9s %9s %9s %9s %12s\n', 'x mm','z mm','|Ex|','|Ez|','ph_x deg','ph_z deg','E_rms','q W/m^3');
for i = 1:size(pts,1)
    [~, s5] = E_field_2d(pts(i,1), pts(i,2), 0, P_set);
    fprintf('   %7.1f %7.1f %9.0f %9.0f %9.1f %9.1f %9.0f %12.3e\n', pts(i,1)*1e3, pts(i,2)*1e3, ...
        s5.Ex_amp, s5.Ez_amp, rad2deg(s5.Ex_phase), rad2deg(s5.Ez_phase), s5.E_rms, s5.q);
end
Et5 = E_field_2d(5e-3, -5e-3, deg2rad(0:45:315), P_set);
fprintf('   phase deg :'); fprintf(' %8d', 0:45:315); fprintf('\n');
fprintf('   Ex V/m    :'); fprintf(' %8.0f', Et5.x); fprintf('\n');
fprintf('   Ez V/m    :'); fprintf(' %8.0f', Et5.z); fprintf('\n');
[~, s6] = E_field_2d(5e-3, -5e-3, 0, P_set, 'ke', pi/(2*27e-3));
[~, s7] = E_field_2d(5e-3, -5e-3, 0, P_set, 'L_sheath', 0);
fprintf('   sensitivity at (5,-5) mm: E_rms %.0f baseline | %.0f quarter-wave ke | %.0f no braid arm\n', ...
        S4.E_rms, s6.E_rms, s7.E_rms);

function [Ex, Ez, Hy] = local_EH(x, z, xs, Idl, k, w, h)
mu0 = 1.25663706212e-6;
Ex = zeros(size(x));  Ez = Ex;  Hy = Ex;  Z = z - h;
for s = 1:numel(xs)
    X = x - xs(s);  R = sqrt(X.^2 + Z.^2);  G = exp(-1j*k*R)./(4*pi*R);
    C2 = -k^2./R.^2 + 3j*k./R.^3 + 3./R.^4;   C1 = -(1j*k./R + 1./R.^2);
    pre = -1j*w*mu0*Idl(s);
    Ex = Ex + pre*(G + G.*(C2.*X.*X + C1)/k^2);
    Ez = Ez + pre*(    G.*(C2.*Z.*X)/k^2);
    Hy = Hy - Idl(s)*(1j*k + 1./R).*G.*Z./R;
end
end
