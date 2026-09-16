function P = fit_closed_form()
%FIT_CLOSED_FORM  Reduce the numerical field map to a compact analytic formula
%   that can be evaluated without running the spectral solver.
%
%   Three nested ansatz are fitted and compared.  All use the EXACT complex
%   gel wavenumber k1 = beta - j*alpha, so attenuation is never a free knob.
%
%   MODEL 1 (point source, single ray)
%       E = K (sin psi)^q R^-p exp(-alpha_fit R)
%
%   MODEL 2 (point source + floor image)
%       E = K | g(R_d,psi_d) e^{-j k1 R_d} + Gamma g(R_r,psi_r) e^{-j k1 R_r} |
%       g(R,psi) = (sin psi)^q R^-p
%
%   MODEL 3 (line source + floor image)  <-- recommended
%       E = K | Int_0^L I(s) [ g(R_d,psi_d) e^{-j k1 R_d}
%                            + Gamma g(R_r,psi_r) e^{-j k1 R_r} ] ds |
%       with the SAME current I(s) = sin(k_e (L - s)) used by the full solver,
%       R_d  = |r - (s, h)|          (direct ray from the element at x = s)
%       R_r  = |r - (s, -2d - h)|    (ray from its image in the gel floor)
%   Model 3 is the physically correct structure: the radiator is a 12 mm line,
%   not a point, and the current vanishes at the tip.

here = fileparts(mfilename('fullpath'));  root = fileparts(here);
M = load(fullfile(root,'results','field_model.mat'));

xc = M.C.L/2;  zc = M.C.h;  d = M.C.d;  L = M.C.L;  ke = M.C.ke;
k1 = M.beta - 1j*M.alpha;

[XX,ZZ] = meshgrid(M.x, M.z);
Rd = sqrt((XX-xc).^2 + (ZZ-zc).^2);        spd = abs(ZZ-zc)./Rd;
Rr = sqrt((XX-xc).^2 + (ZZ+2*d+zc).^2);    spr = abs(ZZ+2*d+zc)./Rr;

Rmin = 4e-3;  Rmax = 28e-3;
msk = Rd >= Rmin & Rd <= Rmax & ZZ <= -1e-3 & abs(XX) <= 45e-3 & spd > 0.05;
dat = M.Erms(msk);  y = log(dat);  n = numel(y);

%% ---------------- MODEL 1 ------------------------------------------------
A = [ones(n,1), log(spd(msk)), -log(Rd(msk)), -Rd(msk)];
c = A\y;
K1 = exp(c(1)); q1 = c(2); p1 = c(3); af1 = c(4);
pred1 = K1 * spd.^q1 .* Rd.^(-p1) .* exp(-af1*Rd);
e1 = 20*log10(pred1(msk)./dat);

%% ---------------- MODEL 2 ------------------------------------------------
    function E = f2(v, Rdv, spdv, Rrv, sprv)
        q = v(1); p = v(2); G = v(3)+1j*v(4); K = exp(v(5));
        E = K*abs(spdv.^q .* Rdv.^(-p) .* exp(-1j*k1*Rdv) ...
                + G*(sprv.^q .* Rrv.^(-p) .* exp(-1j*k1*Rrv)));
    end
best2 = []; J2 = inf;
for G0 = [-0.9 -0.4 0.4 0.9]
    v = fminsearch(@(v) rmsdB(f2(v,Rd(msk),spd(msk),Rr(msk),spr(msk)),dat), ...
        [q1 p1 G0 0 log(K1)], optimset('MaxFunEvals',3e4,'MaxIter',3e4, ...
        'TolX',1e-9,'TolFun',1e-9,'Display','off'));
    J = rmsdB(f2(v,Rd(msk),spd(msk),Rr(msk),spr(msk)),dat);
    if J < J2, J2 = J; best2 = v; end
end
pred2 = f2(best2, Rd, spd, Rr, spr);
e2 = 20*log10(pred2(msk)./dat);

%% ---------------- MODEL 3 (line source) ---------------------------------
Nq = 41;                                       % Simpson nodes along the wire
s  = linspace(0, L, Nq);  hs = s(2)-s(1);
wq = ones(1,Nq); wq(2:2:end-1) = 4; wq(3:2:end-2) = 2; wq = wq*hs/3;
Is = sin(ke*(L - s)) .* wq;                    % current x quadrature weight

    function [Ef, Eg] = f3(v, Xv, Zv)
        q = v(1); p = v(2); G = v(3)+1j*v(4); K = exp(v(5));
        acc = zeros(size(Xv));
        for jj = 1:Nq
            dx2 = (Xv - s(jj)).^2;
            Rdj = sqrt(dx2 + (Zv-zc).^2);       spdj = abs(Zv-zc)./Rdj;
            Rrj = sqrt(dx2 + (Zv+2*d+zc).^2);   sprj = abs(Zv+2*d+zc)./Rrj;
            acc = acc + Is(jj)*( spdj.^q .* Rdj.^(-p) .* exp(-1j*k1*Rdj) ...
                             + G*(sprj.^q .* Rrj.^(-p) .* exp(-1j*k1*Rrj)) );
        end
        Ef = K*abs(acc);  Eg = acc;
    end
Xm = XX(msk); Zm = ZZ(msk);
best3 = []; J3 = inf;
for G0 = [-0.9 -0.4 0.4 0.9]
    v = fminsearch(@(v) rmsdB(f3(v,Xm,Zm),dat), [1 1 G0 0 log(K1)], ...
        optimset('MaxFunEvals',6e3,'MaxIter',6e3,'TolX',1e-8,'TolFun',1e-8, ...
                 'Display','off'));
    J = rmsdB(f3(v,Xm,Zm),dat);
    if J < J3, J3 = J; best3 = v; end
end
pred3 = f3(best3, XX, ZZ);
e3 = 20*log10(pred3(msk)./dat);

%% ---------------- report ------------------------------------------------
fprintf('\n=========== CLOSED-FORM REDUCTION ===========\n');
fprintf('fit region: %.0f <= R <= %.0f mm, depth >= 1 mm, |x| <= 45 mm  (%d pts)\n', ...
        Rmin*1e3, Rmax*1e3, n);
show('MODEL 1 point source, no image', e1);
fprintf('        K=%.4g  q=%.3f  p=%.3f  alpha_fit=%.2f Np/m (exact %.2f)\n', ...
        K1,q1,p1,af1,M.alpha);
show('MODEL 2 point source + image ', e2);
fprintf('        K=%.4g  q=%.3f  p=%.3f  Gamma=%.3f%+.3fj\n', ...
        exp(best2(5)),best2(1),best2(2),best2(3),best2(4));
show('MODEL 3 line source + image  ', e3);
fprintf('        K=%.4g  q=%.3f  p=%.3f  Gamma=%.3f%+.3fj  (|Gamma|=%.3f)\n', ...
        exp(best3(5)),best3(1),best3(2),best3(3),best3(4), abs(best3(3)+1j*best3(4)));
fprintf('        k1 = %.2f - j%.2f rad/m (exact, NOT fitted)\n', real(k1), -imag(k1));

P = struct('model1',pack(K1,q1,p1,[],af1,e1), ...
           'model2',pack(exp(best2(5)),best2(1),best2(2),best2(3)+1j*best2(4),[],e2), ...
           'model3',pack(exp(best3(5)),best3(1),best3(2),best3(3)+1j*best3(4),[],e3), ...
           'k1',k1,'L',L,'ke',ke,'h',zc,'d',d,'x_c',xc, ...
           'Nq',Nq,'Rmin',Rmin,'Rmax',Rmax,'Pdel',M.C.Pdel, ...
           'alpha_exact',M.alpha,'beta',M.beta);
save(fullfile(root,'results','closed_form.mat'),'-struct','P');

%% ---------------- figure -----------------------------------------------
fd = fullfile(root,'figures');
sel = abs(M.x*1e3) <= 45;  mx = max(M.Erms(:));
fh = figure('Visible','off','Position',[100 100 1250 700]);
dat4 = {M.Erms, pred1, pred2, pred3};
ttl  = {'numerical','Model 1 point','Model 2 point+image','Model 3 line+image'};
for q = 1:4
    subplot(2,4,q);
    imagesc(M.x(sel)*1e3, M.z*1e3, 20*log10(dat4{q}(:,sel)/mx));
    set(gca,'YDir','normal'); clim([-45 0]); colormap(turbo);
    title(ttl{q}); xlabel('x (mm)'); if q==1, ylabel('z (mm)'); end
end
for q = 1:3
    subplot(2,4,4+q);
    er = 20*log10(dat4{q+1}./M.Erms); er(~msk) = NaN;
    imagesc(M.x(sel)*1e3, M.z*1e3, er(:,sel)); set(gca,'YDir','normal');
    clim([-5 5]); colormap(gca,'jet'); colorbar;
    title(sprintf('%s err (dB)', ttl{q+1})); xlabel('x (mm)');
    if q==1, ylabel('z (mm)'); end
end
subplot(2,4,8);
[~,ix] = min(abs(M.x - 6e-3));
semilogy(-M.z*1e3, M.Erms(:,ix),'k-','LineWidth',2); hold on;
semilogy(-M.z*1e3, pred1(:,ix),'r--'); semilogy(-M.z*1e3, pred2(:,ix),'g:','LineWidth',1.4);
semilogy(-M.z*1e3, pred3(:,ix),'b-.','LineWidth',1.5);
grid on; xlabel('depth (mm)'); ylabel('|E|_{rms} (V/m)');
legend('numerical','M1','M2','M3','Location','southwest'); title('cut at x = 6 mm');
exportgraphics(fh, fullfile(fd,'fig9_closed_form_fit.png'),'Resolution',150);
close(fh);
fprintf('\nsaved results/closed_form.mat and figures/fig9_closed_form_fit.png\n');
end

% -------------------------------------------------------------------------
function J = rmsdB(pred, dat)
J = sqrt(mean((20*log10(pred./dat)).^2));
if ~isfinite(J), J = 1e6; end
end

function show(nm, e)
fprintf('  %s : rms %5.2f dB | max %5.2f dB | <1dB %3.0f%% | <2dB %3.0f%%\n', ...
        nm, sqrt(mean(e.^2)), max(abs(e)), 100*mean(abs(e)<1), 100*mean(abs(e)<2));
end

function s = pack(K,q,p,G,af,e)
s = struct('K',K,'q',q,'p',p,'Gamma',G,'alpha_fit',af, ...
    'rms_dB',sqrt(mean(e.^2)),'max_dB',max(abs(e)), ...
    'within1dB',mean(abs(e)<1),'within2dB',mean(abs(e)<2));
end
