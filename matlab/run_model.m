function R = run_model(varargin)
%RUN_MODEL  Baseline 2.45 GHz field model for a 12 mm monopole lying on a
%           30 mm deep Phytagel slab.  Produces the phasor field on the
%           y = 0 plane, the derived time-domain quantities, figures, and
%           results/field_model.mat used by E_gel.m.

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
if ~exist(fullfile(root,'results'),'dir'), mkdir(fullfile(root,'results')); end
if ~exist(fullfile(root,'figures'),'dir'), mkdir(fullfile(root,'figures')); end

C = model_config(varargin{:});
fprintf('\n=========== BASELINE MODEL ===========\n');
print_config(C);

cfg = struct('f',C.f,'epsr',C.epsr,'d',C.d,'h',C.h,'L',C.L,'ke',C.ke, ...
             'Idist',C.Idist,'Nfft',C.Nfft,'dx',C.dx,'z',C.z,'Pdel',C.Pdel);
t0 = tic;  o = layered_field(cfg);  fprintf('solver time %.1f s\n', toc(t0));

% ---------------- derived time-domain descriptors -------------------------
% Instantaneous field  E(r,t) = Re{ Etil(r) exp(j w t) }, so
%   |E(r,t)|^2 = 1/2 [ |Etil|^2 + Re{ Q(r) exp(j 2 w t) } ],  Q = Etil . Etil
% (un-conjugated dot product).  The magnitude therefore oscillates at 2w
% between the semi-axes of the polarisation ellipse.
E2  = abs(o.Ex).^2 + abs(o.Ey).^2 + abs(o.Ez).^2;     % |Etil|^2
Q   = o.Ex.^2 + o.Ey.^2 + o.Ez.^2;                    % complex
Erms  = sqrt(E2/2);
Epeak = sqrt((E2 + abs(Q))/2);
Emin  = sqrt(max(E2 - abs(Q),0)/2);
AR    = Emin./Epeak;                                   % axial ratio (0 = linear)

rho   = 1000;                     % [kg/m^3] gel density ~ water
SAR   = 0.5*o.sigma1*E2/rho;      % [W/kg]
Qv    = 0.5*o.sigma1*E2;          % [W/m^3] volumetric power deposition

% ---------------- scalar diagnostics --------------------------------------
xg = o.x(:).';  zg = o.z(:);
Ifeed_unit = feed_current(C);
Rin = 2*o.Ptot/abs(Ifeed_unit)^2;          % input resistance at the feed
Gam = abs((Rin-50)/(Rin+50));              % |Gamma| if X_in were 0

k1 = o.k(2);  alpha = -imag(k1);  beta = real(k1);
lam1 = 2*pi/beta;

fprintf('\n--- medium ---\n');
fprintf('  eps_r(gel)      = %.3f - j%.3f  (tan d = %.4f)\n', ...
        real(C.epsr_gel), -imag(C.epsr_gel), -imag(C.epsr_gel)/real(C.epsr_gel));
fprintf('  sigma_eff       = %.3f S/m\n', o.sigma1);
fprintf('  alpha           = %.2f Np/m  -> field skin depth  %.2f mm\n', alpha, 1e3/alpha);
fprintf('                                   power  1/e depth  %.2f mm\n', 1e3/(2*alpha));
fprintf('  beta            = %.2f rad/m -> wavelength in gel %.2f mm\n', beta, lam1*1e3);
fprintf('  L / lambda_gel  = %.3f\n', C.L/lam1);

fprintf('\n--- power budget (delivered %.0f W) ---\n', C.Pdel);
fprintf('  absorbed in gel        %8.2f W  (%5.1f %%)\n', o.Pgel*o.scale^2, 100*o.Pgel/o.Ptot);
fprintf('  radiated up into air   %8.2f W  (%5.1f %%)\n', o.Pup *o.scale^2, 100*o.Pup /o.Ptot);
fprintf('  out through the bottom %8.2f W  (%5.1f %%)\n', o.Pdn *o.scale^2, 100*o.Pdn /o.Ptot);
fprintf('  antenna current amplitude I_m = %.3f A, feed current %.3f A\n', ...
        o.scale, o.scale*abs(Ifeed_unit));
fprintf('  implied input RESISTANCE R_in = %.1f ohm  -> |Gamma| >= %.3f (X_in=0)\n', Rin, Gam);
fprintf('  (generator spec P_ref/P_fwd = %.2f -> |Gamma| = %.2f)\n', ...
        C.Pref/C.Pfwd, sqrt(C.Pref/C.Pfwd));

% field magnitudes at reference points
fprintf('\n--- E-field (RMS) at reference points, y = 0 ---\n');
probes = [ 6  0.5 ; 6 2 ; 6 5 ; 6 10 ; 6 20 ; 6 29 ; 0 5 ; 12 5 ; 25 5 ; 50 5];
fprintf('     x[mm] depth[mm]   E_rms[V/m]  E_peak[V/m]  SAR[W/kg]   axial ratio\n');
for q = 1:size(probes,1)
    [ix,iz] = nearest_idx(xg, zg, probes(q,1)*1e-3, -probes(q,2)*1e-3);
    fprintf('   %6.1f %8.1f  %11.1f %12.1f %11.1f %13.3f\n', probes(q,1), probes(q,2), ...
        Erms(iz,ix), Epeak(iz,ix), SAR(iz,ix), AR(iz,ix));
end

% penetration metrics along the vertical through the middle of the antenna
[ixm,~] = nearest_idx(xg, zg, C.L/2, -1e-3);
prof = Erms(:,ixm);
z_1e  = interp_cross(zg, prof, prof(1)/exp(1));
z_10  = interp_cross(zg, prof, prof(1)/10);
fprintf('\n--- penetration on the vertical x = %.1f mm ---\n', 1e3*xg(ixm));
fprintf('  E_rms at 0.125 mm depth = %.0f V/m\n', prof(1));
fprintf('  depth where E falls to 1/e  of surface value : %.2f mm\n', -1e3*z_1e);
fprintf('  depth where E falls to 1/10 of surface value : %.2f mm\n', -1e3*z_10);
fprintf('  E_rms at the gel floor (z = -%.0f mm)         : %.1f V/m (%.2f %% of surface)\n', ...
        1e3*C.d, prof(end), 100*prof(end)/prof(1));

% fraction of the deposited power in the top 10 mm
dz = abs(zg(2)-zg(1));
Pz = sum(Qv,2)*C.dx*dz;                     % per unit y-width, arbitrary scale
cum = cumsum(Pz)/sum(Pz);
i10 = find(zg <= -10e-3, 1);
fprintf('  share of y=0-plane power deposition above 10 mm depth: %.1f %%\n', 100*cum(i10));

% ---------------- save ----------------------------------------------------
R = struct('C',C,'x',xg,'z',zg,'Ex',o.Ex,'Ey',o.Ey,'Ez',o.Ez, ...
    'E2',E2,'Q',Q,'Erms',Erms,'Epeak',Epeak,'Emin',Emin,'AR',AR, ...
    'SAR',SAR,'Qv',Qv,'sigma',o.sigma1,'k',o.k,'alpha',alpha,'beta',beta, ...
    'scale',o.scale,'Pgel',o.Pgel*o.scale^2,'Pup',o.Pup*o.scale^2, ...
    'Pdn',o.Pdn*o.scale^2,'Rin',Rin,'omega',2*pi*C.f,'slice',o.slice);
save(fullfile(root,'results','field_model.mat'),'-struct','R');
fprintf('\nsaved results/field_model.mat\n');

make_figures(R, root);
end

% =========================================================================
function print_config(C)
fprintf('  f = %.3f GHz | L = %.1f mm | gel depth = %.1f mm | h = %.2f mm\n', ...
    C.f/1e9, C.L*1e3, C.d*1e3, C.h*1e3);
fprintf('  P_fwd = %.0f W, P_ref(max) = %.0f W -> P_delivered = %.0f W\n', ...
    C.Pfwd, C.Pref, C.Pdel);
fprintf('  current: %s, k_e case ''%s'' -> k_e = %.1f rad/m (L = %.3f lambda_e)\n', ...
    C.Idist, C.ke_case, C.ke, C.L*C.ke/(2*pi));
fprintf('  grid: N = %d, dx = %.3f mm (window +-%.0f mm), dz = %.3f mm\n', ...
    C.Nfft, C.dx*1e3, C.Nfft*C.dx/2*1e3, C.dz*1e3);
end

function I0 = feed_current(C)
switch lower(C.Idist)
    case 'sinusoid',   I0 = sin(C.ke*C.L);
    case 'uniform',    I0 = 1;
    case 'triangular', I0 = 1;
    case 'hertzian',   I0 = 1;
end
if abs(I0) < 1e-6, I0 = 1e-6; end
end

function [ix,iz] = nearest_idx(xg, zg, xq, zq)
[~,ix] = min(abs(xg - xq));
[~,iz] = min(abs(zg - zq));
end

function zc = interp_cross(z, v, target)
i = find(v <= target, 1);
if isempty(i), zc = z(end); return; end
if i == 1, zc = z(1); return; end
w = (v(i-1)-target)/(v(i-1)-v(i));
zc = z(i-1) + w*(z(i)-z(i-1));
end

% =========================================================================
function make_figures(R, root)
fd = fullfile(root,'figures');
xm = R.x*1e3;  zm = R.z*1e3;   L = R.C.L*1e3;
sel = abs(xm) <= 45;
X = xm(sel);

% ---- Fig 1: RMS |E| map, dB ---------------------------------------------
f1 = figure('Visible','off','Position',[100 100 900 420]);
V  = 20*log10(R.Erms(:,sel)/max(R.Erms(:)));
imagesc(X, zm, V); set(gca,'YDir','normal'); clim([-45 0]);
colormap(turbo); cb = colorbar; cb.Label.String = '|E|_{rms} (dB re max)';
hold on; plot([0 L],[0 0],'w-','LineWidth',4); plot([0 L],[0 0],'k-','LineWidth',1.5);
yline(-R.C.d*1e3,'w--','LineWidth',1.5);
xlabel('x  (mm)  -- along the antenna axis'); ylabel('z  (mm)  -- depth into gel');
title(sprintf(['|E|_{rms} on the y = 0 plane, %.2f GHz, %.0f W delivered ' ...
    '(white bar = 12 mm monopole)'], R.C.f/1e9, R.C.Pdel));
exportgraphics(f1, fullfile(fd,'fig1_Erms_dB.png'), 'Resolution',160); close(f1);

% ---- Fig 2: linear |E| with contours ------------------------------------
f2 = figure('Visible','off','Position',[100 100 900 420]);
imagesc(X, zm, R.Erms(:,sel)); set(gca,'YDir','normal');
clim([0 prctile(R.Erms(:),99.5)]);
colormap(parula); cb = colorbar; cb.Label.String = '|E|_{rms}  (V/m)';
hold on;
lev = [200 500 1000 2000 5000 10000];
[c1,h1] = contour(X, zm, R.Erms(:,sel), lev, 'w-','LineWidth',0.8);
clabel(c1,h1,'Color','w','FontSize',8);
plot([0 L],[0 0],'r-','LineWidth',4);
xlabel('x  (mm)'); ylabel('z  (mm)');
title('|E|_{rms} (V/m), y = 0 plane, contours labelled in V/m');
exportgraphics(f2, fullfile(fd,'fig2_Erms_linear.png'), 'Resolution',160); close(f2);

% ---- Fig 3: SAR ----------------------------------------------------------
f3 = figure('Visible','off','Position',[100 100 900 420]);
imagesc(X, zm, log10(max(R.SAR(:,sel),1e-2))); set(gca,'YDir','normal');
colormap(hot); cb = colorbar; cb.Label.String = 'log_{10} SAR  (W/kg)';
hold on; plot([0 L],[0 0],'c-','LineWidth',4);
xlabel('x  (mm)'); ylabel('z  (mm)');
title(sprintf('Specific absorption rate, sigma = %.2f S/m, rho = 1000 kg/m^3', R.sigma));
exportgraphics(f3, fullfile(fd,'fig3_SAR.png'), 'Resolution',160); close(f3);

% ---- Fig 4: components ---------------------------------------------------
f4 = figure('Visible','off','Position',[100 100 1300 350]);
comp = {abs(R.Ex(:,sel))/sqrt(2), abs(R.Ey(:,sel))/sqrt(2), abs(R.Ez(:,sel))/sqrt(2)};
nm = {'|E_x|_{rms}','|E_y|_{rms}','|E_z|_{rms}'};
mx = max(cellfun(@(v) max(v(:)), comp));
for q = 1:3
    subplot(1,3,q);
    imagesc(X, zm, 20*log10(comp{q}/mx)); set(gca,'YDir','normal'); clim([-45 0]);
    colormap(turbo); hold on; plot([0 L],[0 0],'w-','LineWidth',3);
    title(nm{q}); xlabel('x (mm)'); if q==1, ylabel('z (mm)'); end
end
cb = colorbar; cb.Label.String = 'dB re max component';
exportgraphics(f4, fullfile(fd,'fig4_components.png'), 'Resolution',160); close(f4);

% ---- Fig 5: depth profiles ----------------------------------------------
f5 = figure('Visible','off','Position',[100 100 780 460]);
xs = [0 3 6 9 12 18 25 40];
co = lines(numel(xs)); lg = cell(1,numel(xs));
for q = 1:numel(xs)
    [~,ix] = min(abs(xm - xs(q)));
    semilogy(-zm, R.Erms(:,ix), 'LineWidth',1.6, 'Color',co(q,:)); hold on;
    lg{q} = sprintf('x = %d mm', xs(q));
end
grid on; xlabel('depth into gel  (mm)'); ylabel('|E|_{rms}  (V/m)');
title(sprintf('Depth profiles, y = 0, %.0f W delivered', R.C.Pdel));
legend(lg,'Location','northeast'); xlim([0 R.C.d*1e3]);
exportgraphics(f5, fullfile(fd,'fig5_depth_profiles.png'), 'Resolution',160); close(f5);

% ---- Fig 6: instantaneous field snapshots -------------------------------
f6 = figure('Visible','off','Position',[100 100 1100 620]);
ph = [0 pi/4 pi/2 3*pi/4];
for q = 1:4
    subplot(2,2,q);
    ex = real(R.Ex(:,sel)*exp(1j*ph(q)));
    ez = real(R.Ez(:,sel)*exp(1j*ph(q)));
    ey = real(R.Ey(:,sel)*exp(1j*ph(q)));
    Ei = sqrt(ex.^2+ey.^2+ez.^2);
    imagesc(X, zm, Ei); set(gca,'YDir','normal'); clim([0 prctile(R.Epeak(:),99.5)]);
    colormap(parula); hold on; plot([0 L],[0 0],'r-','LineWidth',3);
    title(sprintf('|E(r,t)| at \\omega t = %.2f rad', ph(q)));
    if q>2, xlabel('x (mm)'); end; if mod(q,2)==1, ylabel('z (mm)'); end
end
sgtitle(['Instantaneous |E| at \omega t = 0 .. 3\pi/4  (drive period 408 ps; ' ...
         '|E| itself repeats every 204 ps)']);
exportgraphics(f6, fullfile(fd,'fig6_time_snapshots.png'), 'Resolution',160); close(f6);

% ---- Fig 7: polarisation / axial ratio ----------------------------------
f7 = figure('Visible','off','Position',[100 100 900 420]);
imagesc(X, zm, R.AR(:,sel)); set(gca,'YDir','normal'); clim([0 1]);
colormap(parula); cb = colorbar;
cb.Label.String = 'axial ratio  E_{min}/E_{peak}  (0 = linear, 1 = circular)';
hold on; plot([0 L],[0 0],'r-','LineWidth',4);
xlabel('x (mm)'); ylabel('z (mm)');
title('Polarisation ellipticity of E(r,t) -- how much |E| ripples in time');
exportgraphics(f7, fullfile(fd,'fig7_axial_ratio.png'), 'Resolution',160); close(f7);

% ---- Fig 8: |E|(t) at three points --------------------------------------
f8 = figure('Visible','off','Position',[100 100 780 420]);
T = 1/R.C.f;  t = linspace(0, 2*T, 600);
pts = [6 2; 6 10; 6 25];  co = lines(3); lg = {};
for q = 1:3
    [~,ix] = min(abs(xm - pts(q,1)));  [~,iz] = min(abs(-zm - pts(q,2)));
    E2p = R.E2(iz,ix); Qp = R.Q(iz,ix);
    Et = sqrt(max(0.5*(E2p + real(Qp*exp(2j*R.omega*t))),0));
    plot(t*1e12, Et, 'LineWidth',1.6,'Color',co(q,:)); hold on;
    lg{q} = sprintf('x=%g mm, depth=%g mm', pts(q,1), pts(q,2)); %#ok<AGROW>
end
grid on; xlabel('t  (ps)'); ylabel('|E(r,t)|  (V/m)'); set(gca,'YScale','log');
title('Instantaneous field magnitude vs time (ripples at 2f = 4.9 GHz)');
legend(lg,'Location','best');
exportgraphics(f8, fullfile(fd,'fig8_E_vs_time.png'), 'Resolution',160); close(f8);

fprintf('wrote 8 figures to figures/\n');
end
