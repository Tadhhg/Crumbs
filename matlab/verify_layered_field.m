function verify_layered_field()
%VERIFY_LAYERED_FIELD  Validation suite for the spectral layered-media solver.
%
%   TEST 1  Homogeneous limit.  With all three layers set to the same lossy
%           medium the numerical solver must reproduce the closed-form
%           Hertzian-dipole field of HERTZIAN_DIPOLE_FIELD exactly.
%   TEST 2  Grid convergence (spatial sample dx).
%   TEST 3  Window convergence (FFT window size).
%   TEST 4  Asymptotic depth decay must approach the plane-wave attenuation
%           constant alpha of the lossy medium.
%   TEST 5  Power-budget consistency: the three power channels must sum to the
%           normalisation target, and the split must be insensitive to the grid.

fprintf('\n================ VALIDATION OF layered_field.m ================\n');
f = 2.45e9;
epsg = saline_stogryn(f, 25, 0.154);      % 0.9%% saline at 25 degC

%% ---------------- TEST 1: homogeneous limit ---------------------------
fprintf('\n--- TEST 1: homogeneous limit vs analytic Hertzian dipole ---\n');
cfg = base_cfg(f);
cfg.epsr  = [epsg epsg epsg];             % NO interfaces
cfg.Idist = 'hertzian';
cfg.L     = 1e-4;                         % dipole moment I*dl, I = 1 A
cfg.h     = 1.0e-3;
cfg.Nfft  = 1024;  cfg.dx = 0.25e-3;
cfg.z     = -(0.5e-3:0.25e-3:20e-3);
cfg.Pdel  = 1;

o = layered_field(cfg);
o = unscale(o);                           % remove the sqrt(Pdel/Ptot) factor

% (a) vertical line x = 0, y = 0  -- pure E_x (broadside), tests TE+TM sum
ix0 = cfg.Nfft/2 + 1;
Zc  = o.z - cfg.h;                        % source sits at z = +h
[ax,ay,az] = hertzian_dipole_field(f, epsg, cfg.L, 0*Zc, 0*Zc, Zc);
report('vertical line x=y=0', [o.Ex(:,ix0) o.Ey(:,ix0) o.Ez(:,ix0)], [ax(:) ay(:) az(:)]);

% (b) horizontal line along x at z = -5 mm -- endfire cut, strong E_z (TM)
iz  = find(abs(o.z + 5e-3) == min(abs(o.z + 5e-3)), 1);
sel = abs(o.x) <= 40e-3;
xs  = o.x(sel).';
zq  = (o.z(iz)-cfg.h) + 0*xs;
[ax,ay,az] = hertzian_dipole_field(f, epsg, cfg.L, xs, 0*xs, zq);
report('line along x, z=-5 mm', ...
       [o.Ex(iz,sel).' o.Ey(iz,sel).' o.Ez(iz,sel).'], [ax ay az]);

% (c) line along y at x = 0, mid-depth -- tests the TE branch specifically
izs  = find(abs(o.z - o.slice.z) < 1e-12, 1);
ys   = xs;
zq   = (o.z(izs)-cfg.h) + 0*ys;
[ax,ay,az] = hertzian_dipole_field(f, epsg, cfg.L, 0*ys, ys, zq);
report('line along y, x=0', [o.slice.Ex(ix0,sel).' o.slice.Ey(ix0,sel).' ...
                             o.slice.Ez(ix0,sel).'], [ax ay az]);

% (d) off-axis diagonal point set -- tests the full TE/TM mixing
[ax,ay,az] = hertzian_dipole_field(f, epsg, cfg.L, xs, xs, zq);
dg = zeros(numel(xs),3);
for q = 1:numel(xs)
    iX = find(abs(o.x - xs(q)) == min(abs(o.x - xs(q))),1);
    iY = find(abs(o.y - xs(q)) == min(abs(o.y - xs(q))),1);
    dg(q,:) = [o.slice.Ex(iX,iY) o.slice.Ey(iX,iY) o.slice.Ez(iX,iY)];
end
report('diagonal x=y line', dg, [ax ay az]);

%% ---------------- TEST 2/3: grid and window convergence ---------------
fprintf('\n--- TEST 2/3: grid and window convergence (full 3-layer problem) ---\n');
ref = [];
combos = { {1024, 0.25e-3}, {2048, 0.25e-3}, {2048, 0.125e-3}, {1536, 0.20e-3} };
for c = 1:numel(combos)
    cfg = base_cfg(f);
    cfg.epsr = [1 epsg 1];
    cfg.Nfft = combos{c}{1};  cfg.dx = combos{c}{2};
    o = layered_field(cfg);
    ix = find(abs(o.x - 6e-3) == min(abs(o.x - 6e-3)),1);
    Em = sqrt(abs(o.Ex(:,ix)).^2 + abs(o.Ey(:,ix)).^2 + abs(o.Ez(:,ix)).^2);
    probe = [Em(4) Em(20) Em(60)];
    if isempty(ref), ref = probe; end
    fprintf(['  N=%4d dx=%.3f mm  win=+-%5.1f mm | Pgel=%6.2f Pup=%6.2f ' ...
             'Pdn=%5.2f W | |E| probes %8.1f %8.1f %8.1f V/m (dev %5.2f%%)\n'], ...
        cfg.Nfft, cfg.dx*1e3, cfg.Nfft*cfg.dx/2*1e3, o.Pgel, o.Pup, o.Pdn, ...
        probe, 100*max(abs(probe-ref)./abs(ref)));
end

%% ---------------- TEST 4: asymptotic attenuation ----------------------
fprintf('\n--- TEST 4: deep-field decay vs plane-wave attenuation ---\n');
cfg = base_cfg(f);  cfg.epsr = [1 epsg 1];
cfg.Nfft = 2048;  cfg.dx = 0.25e-3;
o = layered_field(cfg);
k1 = o.k(2);  alpha = -imag(k1);
ix = find(abs(o.x - 6e-3) == min(abs(o.x - 6e-3)),1);
Em = sqrt(abs(o.Ex(:,ix)).^2 + abs(o.Ey(:,ix)).^2 + abs(o.Ez(:,ix)).^2);
zc = o.z;
i1 = find(zc <= -12e-3, 1);  i2 = find(zc <= -22e-3, 1);
alpha_num = -log(Em(i2)/Em(i1))/abs(zc(i2)-zc(i1));
fprintf('  plane-wave alpha       = %7.2f Np/m  (skin depth %5.2f mm)\n', alpha, 1e3/alpha);
fprintf('  numerical decay 12-22mm= %7.2f Np/m  (includes 1/r spreading)\n', alpha_num);
fprintf('  ratio                  = %7.3f  (expect >1: geometric spreading adds)\n', alpha_num/alpha);

%% ---------------- TEST 5: power budget --------------------------------
fprintf('\n--- TEST 5: power budget closure ---\n');
fprintf('  Pgel+Pup+Pdn = %.6f W ; normalisation target Pdel = %.6f W\n', ...
        o.Ptot*o.scale^2, cfg.Pdel);
fprintf('  fractions: gel %.1f%%  up(air) %.1f%%  through-bottom %.1f%%\n', ...
        100*o.Pgel/o.Ptot, 100*o.Pup/o.Ptot, 100*o.Pdn/o.Ptot);
fprintf('\n================ END VALIDATION ================\n\n');
end

% -------------------------------------------------------------------------
function cfg = base_cfg(f)
epsg = saline_stogryn(f, 25, 0.154);
cfg = struct('f',f,'epsr',[1 epsg 1],'d',30e-3,'h',0.5e-3,'L',12e-3, ...
    'ke',2*pi*f/299792458*sqrt((1+real(epsg))/2),'Idist','sinusoid', ...
    'Nfft',1024,'dx',0.25e-3,'z',-(0.125e-3:0.25e-3:29.875e-3),'Pdel',150);
end

function o = unscale(o)
o.Ex = o.Ex/o.scale; o.Ey = o.Ey/o.scale; o.Ez = o.Ez/o.scale;
if isfield(o,'slice')
    o.slice.Ex = o.slice.Ex/o.scale; o.slice.Ey = o.slice.Ey/o.scale;
    o.slice.Ez = o.slice.Ez/o.scale;
end
end

function report(name, num, ana)
err = abs(num - ana);
den = max(abs(ana(:)));
fprintf('  %-26s  max rel err = %9.3e   (rms %9.3e)\n', name, ...
    max(err(:))/den, sqrt(mean(err(:).^2))/den);
end
