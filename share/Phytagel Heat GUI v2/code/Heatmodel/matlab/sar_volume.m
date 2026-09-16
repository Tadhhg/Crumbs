function out = sar_volume(cfg)
%SAR_VOLUME  Volumetric microwave power deposition q(x,y,z) [W/m^3] on a
%            3-D block of the gel, for the 12 mm monopole of the field model.
%
%   This is the SAME spectral (Sommerfeld / Weyl) three-layer solver as
%   ../../matlab/layered_field.m -- same Weyl expansion, same TE/TM
%   transmission-line recursion, same branch handling, same FFT
%   normalisation, same delivered-power normalisation.  The ONLY differences
%   are in what it keeps:
%
%     layered_field.m  keeps the complex phasor field on the y = 0 plane
%                      (plus one horizontal slice) and discards the rest;
%     sar_volume.m     keeps  q = 0.5*sigma_eff*|Etilde|^2  on a 3-D block
%                      (the thermal region of interest) and discards the
%                      phasors outside the y = 0 plane.
%
%   Duplicating the solver rather than re-running the old one is deliberate:
%   the 3-D block is what the heat equation needs, and holding the full 3-D
%   complex phasor field at 0.25 mm over the whole 256 mm FFT window would be
%   ~200 GB.  The duplication is CHECKED, not assumed: verify_heat_model.m
%   compares the y = 0 plane of this solver against E_gel.m point by point
%   and against the stored delivered-power normalisation of the field model.
%
%   WHY q AND NOT SAR.  SAR [W/kg] = q/rho is the dosimetry quantity; the heat
%   equation needs q [W/m^3].  They differ only by rho, kept at 1000 kg/m^3
%   in both models.
%
%   TIME AVERAGING.  The heat capacity of the gel is enormous compared with
%   one microwave period (204 ps for |E|, cf. a thermal time constant of
%   seconds), so only the CYCLE-AVERAGED deposition matters:
%       q = 0.5*sigma_eff*|Etilde|^2 = sigma_eff*E_rms^2 .
%   The 4.9 GHz ripple in |E(t)| derived in the field model averages out
%   exactly and cannot appear in the temperature.
%
%   cfg FIELDS  (all those of layered_field.m, plus)
%     zn     vector of z NODE positions to evaluate, 0 >= zn >= -d.  Unlike
%            layered_field.m these are nodes, not cell centres, so that they
%            coincide with the thermal grid.  Must be uniformly spaced with
%            an EVEN number of intervals (Simpson).
%     xlim   [xmin xmax] of the block  [m]
%     ylim   [ymin ymax] of the block  [m]
%
%   OUTPUT struct
%     q      [nx x ny x nz] volumetric power deposition [W/m^3] at Pdel
%     x,y,z  the block node coordinates [m]
%     Ex,Ey,Ez  complex phasors on the y = 0 plane, [nz x Nfft], for
%            cross-validation against E_gel.m
%     Pgel, Pup, Pdn, Ptot, scale, sigma1, k

eps0 = 8.8541878128e-12;          % [F/m]  CODATA 2018
mu0  = 1.25663706212e-6;          % [H/m]  CODATA 2018
w    = 2*pi*cfg.f;

N  = cfg.Nfft;  dx = cfg.dx;  dy = dx;
x  = (-N/2:N/2-1)*dx;   y = x;

% ---- spectral grid (rad/m), FFT ordering ---------------------------------
kv = 2*pi/(N*dx) * ifftshift(-N/2:N/2-1);
[KX,KY] = ndgrid(kv,kv);
KR2 = KX.^2 + KY.^2;   KR = sqrt(KR2);

% ---- layer wavenumbers and vertical wavenumbers --------------------------
k = w*sqrt(mu0*eps0*cfg.epsr(:).');             % [k0 k1 k2]
k(imag(k)>0) = conj(k(imag(k)>0));              % enforce Im(k) <= 0
kz = cell(1,3);
for n = 1:3
    tt = sqrt(k(n)^2 - KR2);
    tt(imag(tt)>0) = -tt(imag(tt)>0);           % enforce Im(kz) <= 0
    kz{n} = tt;
end
kz0 = kz{1}; kz1 = kz{2};

% ---- modal wave impedances -----------------------------------------------
epsC = eps0*cfg.epsr(:).';
Z0te = w*mu0./kz{1}; Z1te = w*mu0./kz{2}; Z2te = w*mu0./kz{3};
Z0tm = kz{1}/(w*epsC(1)); Z1tm = kz{2}/(w*epsC(2)); Z2tm = kz{3}/(w*epsC(3));

% ---- interface coefficients ----------------------------------------------
r01a = (Z1te - Z0te)./(Z1te + Z0te);   r12a = (Z2te - Z1te)./(Z2te + Z1te);
r01b = (Z0tm - Z1tm)./(Z0tm + Z1tm);   r12b = (Z1tm - Z2tm)./(Z1tm + Z2tm);
t01a = 1 + r01a;                       t01b = 1 + r01b;

d    = cfg.d;
tau2 = exp(-2j*kz1*d);
Da_den = 1 + r01a.*r12a.*tau2;
Db_den = 1 + r01b.*r12b.*tau2;
Ra = (r01a + r12a.*tau2)./Da_den;
Rb = (r01b + r12b.*tau2)./Db_den;

% ---- source spectrum ------------------------------------------------------
It = current_spectrum(kv, cfg);
IT = repmat(It(:), 1, N);
% Kept for VERIFY_V2 test 6b.  A jump of dI at the feed puts a delta in the
% line charge, and integrating by parts shows that its whole signature is
% kx*Itilde(kx) -> dI as kx -> infinity.  A continuous I gives kx*Itilde -> 0
% (Riemann-Lebesgue).  So the spectral tail decides whether the point charge
% is there, exactly, with no grid and no threshold to argue about -- which is
% better than inferring it from a peak that refuses to converge.
if isfield(cfg,'debug_current') && cfg.debug_current
    dbg_It = It;  dbg_kv = kv;
else
    dbg_It = [];  dbg_kv = [];
end

KRg = KR; KRg(KRg==0) = eps;

A0d = mu0*IT.*exp(-1j*kz0*cfg.h)./(2j*kz0);
a0d = -1j*w*A0d.*(KY./KRg);
b0d = -1j*A0d.*KX.*kz0./(mu0*KRg);
A0u = mu0*IT.*exp(+1j*kz0*cfg.h)./(2j*kz0);
a0u = -1j*w*A0u.*(KY./KRg);
b0u = +1j*A0u.*KX.*kz0./(mu0*KRg);

Da = t01a.*a0d./Da_den;
Db = t01b.*b0d./Db_den;

norm2d = 1/(N*dx*N*dy);

% ---- power absorbed in the bath, and the flux that escapes above it ------
% v2.  With AIR above (v1) the upper half-space is lossless: every plane above
% the source carries the same real flux, so "power radiated upward" was the
% whole story and the plane it was measured on did not matter.
%
% With a LIQUID above it is not.  The bath absorbs, so the flux through a
% plane depends on the plane, and the absorbed part must be counted or the
% delivered-power normalisation below divides by too small a total and
% inflates the field in the gel -- which is precisely the error v2 exists to
% remove.
%
% So integrate the dissipation over a column of bath 0 <= z <= Zw exactly as
% Pgel is integrated over the slab, and take Pup as the RESIDUAL flux leaving
% the top of that column.  Nothing is double counted, and the budget closes
% for any Zw: Zw only decides how much of the loss is called "absorbed in the
% bath" rather than "escaped upward".
%
% The integral is done in the SPECTRAL domain.  By Parseval, with this FFT
% normalisation,  sum|e|^2 dx dy == norm2d * sum|E_spectral|^2 , so no inverse
% transform is needed for a quantity that is only ever integrated.  (That
% identity is checked against the gel loop in verify_v2.m, test 3.)
sigma0 = w*eps0*(-imag(cfg.epsr(1)));           % [S/m] bath conductivity
if ~isfield(cfg,'Zw'),  cfg.Zw  = 40e-3; end    % bath column integrated [m]
if ~isfield(cfg,'dzw'), cfg.dzw =  1e-3; end
ztop = cfg.Zw;
aU = (a0u + Ra.*a0d).*exp(-1j*kz0*ztop);
bU = (b0u + Rb.*b0d).*exp(-1j*kz0*ztop);
Pup = 0.5*real(sum(sum(abs(aU).^2./conj(Z0te) + Z0tm.*abs(bU).^2)))*norm2d;

Pwat = 0;
if sigma0 > 0
    zw = (0:cfg.dzw:cfg.Zw).';
    if mod(numel(zw),2) ~= 1, zw = zw(1:end-1); end     % odd -> Simpson
    wW = ones(1,numel(zw)); wW(2:2:end-1) = 4; wW(3:2:end-2) = 2;
    wW = wW*(zw(2)-zw(1))/3;
    for iw = 1:numel(zw)
        zc = zw(iw);
        [Ex0,Ey0,Ez0] = bath_spectral_E(zc, cfg.h, a0d, b0d, a0u, b0u, ...
                                        Ra, Rb, kz0, KX, KY, KRg, Z0tm, ...
                                        w, epsC(1));
        Pwat = Pwat + 0.5*sigma0*sum(sum(abs(Ex0).^2 + abs(Ey0).^2 + ...
                                         abs(Ez0).^2))*norm2d*wW(iw);
    end
end
% v2 self-check.  Tangential E must be CONTINUOUS across z = 0.  That is the
% sharpest available test of the up-going TM sign in BATH_SPECTRAL_E, so it is
% exposed here rather than argued for in a comment; verify_v2.m test 1 runs it.
if isfield(cfg,'debug_interface') && cfg.debug_interface
    [bx, by] = bath_spectral_E(0, cfg.h, a0d, b0d, a0u, b0u, Ra, Rb, kz0, ...
                               KX, KY, KRg, Z0tm, w, epsC(1));
    adz = Da;  auz = Da.*r12a.*exp(-2j*kz1*d);
    bdz = Db;  buz = Db.*r12b.*exp(-2j*kz1*d);
    Esz = adz + auz;   Bdz = bdz - buz;
    gx =  Esz.*(KY./KRg) + Z1tm.*Bdz.*(KX./KRg);
    gy = -Esz.*(KX./KRg) + Z1tm.*Bdz.*(KY./KRg);
    dbg = struct('Ex_bath',bx,'Ey_bath',by,'Ex_gel',gx,'Ey_gel',gy);
else
    dbg = [];
end
clear aU bU Ex0 Ey0 Ez0 a0u b0u a0d b0d A0d A0u IT Ra Rb Z0te Z0tm kz0 KR KR2

% ---- power transmitted downward out of the slab --------------------------
a2 = (1+r12a).*Da.*exp(-1j*kz1*d);
b2 = (1+r12b).*Db.*exp(-1j*kz1*d);
Pdn = 0.5*real(sum(sum(abs(a2).^2./conj(Z2te) + Z2tm.*abs(b2).^2)))*norm2d;
clear a2 b2 Z2te Z2tm Z1te r01a r01b t01a t01b Da_den Db_den tau2

% ================= evaluate on the node depths ============================
zn = cfg.zn(:);   nz = numel(zn);
if mod(nz,2) ~= 1
    error('sar_volume:simpson', ...
        'numel(zn) must be odd (an even number of intervals) for Simpson');
end
dz = abs(zn(2)-zn(1));

% Simpson weights for the depth integral of the dissipated power
wS = ones(1,nz); wS(2:2:end-1) = 4; wS(3:2:end-2) = 2; wS = wS*dz/3;

% block index ranges on the FFT grid
ix = find(x >= cfg.xlim(1)-1e-12 & x <= cfg.xlim(2)+1e-12);
iy = find(y >= cfg.ylim(1)-1e-12 & y <= cfg.ylim(2)+1e-12);
nx = numel(ix);  ny = numel(iy);

sigma1 = w*eps0*(-imag(cfg.epsr(2)));           % effective conductivity [S/m]

E2blk = zeros(nx,ny,nz);                        % |Etilde|^2 on the block
Ex = zeros(nz,N); Ey = zeros(nz,N); Ez = zeros(nz,N);
iy0 = N/2 + 1;                                  % index of y = 0
Pgel = 0;

for iz = 1:nz
    zc = zn(iz);
    ad = Da.*exp( 1j*kz1*zc);
    au = Da.*r12a.*exp(-1j*kz1*(2*d + zc));     % overflow-safe form
    bd = Db.*exp( 1j*kz1*zc);
    bu = Db.*r12b.*exp(-1j*kz1*(2*d + zc));

    Es = ad + au;  Bs = bd + bu;  Bd = bd - bu;

    Exs =  Es.*(KY./KRg) + Z1tm.*Bd.*(KX./KRg);
    Eys = -Es.*(KX./KRg) + Z1tm.*Bd.*(KY./KRg);
    Ezs =  Bs.*KRg/(w*epsC(2));

    ex = fftshift(fft2(Exs))*norm2d;
    ey = fftshift(fft2(Eys))*norm2d;
    ez = fftshift(fft2(Ezs))*norm2d;

    E2full = abs(ex).^2 + abs(ey).^2 + abs(ez).^2;

    % dissipated power: exact in x,y (Parseval / full grid), Simpson in z
    Pgel = Pgel + 0.5*sigma1*sum(E2full(:))*dx*dy*wS(iz);

    E2blk(:,:,iz) = E2full(ix,iy);
    Ex(iz,:) = ex(:,iy0).';  Ey(iz,:) = ey(:,iy0).';  Ez(iz,:) = ez(:,iy0).';
end

Ptot  = Pgel + Pwat + Pup + Pdn;                % v2: the bath is in the budget
scale = sqrt(cfg.Pdel/Ptot);                    % amplitude scale factor

out = struct();
out.x = x(ix);  out.y = y(iy);  out.z = zn;
out.q = 0.5*sigma1*E2blk*scale^2;               % [W/m^3] at Pdel
out.Ex = Ex*scale;  out.Ey = Ey*scale;  out.Ez = Ez*scale;
out.x_full = x;
out.scale = scale;  out.sigma1 = sigma1;  out.k = k;
out.Pgel = Pgel*scale^2;  out.Pup = Pup*scale^2;  out.Pdn = Pdn*scale^2;
out.Pwat = Pwat*scale^2;                        % v2: absorbed in the bath
out.frac_gel = Pgel/Ptot;                       % v2: was ~0.996 with air above
out.Ptot = Ptot;  out.cfg = cfg;  out.dbg = dbg;
out.It = dbg_It;  out.kv = dbg_kv;

% power captured inside the block (trapezoid in x,y; Simpson in z)
wx = ones(nx,1); wx([1 end]) = 0.5;
wy = ones(1,ny); wy([1 end]) = 0.5;
Wxy = (wx*wy)*dx*dy;
Pblk = 0;
for iz = 1:nz
    Pblk = Pblk + sum(sum(out.q(:,:,iz).*Wxy))*wS(iz);
end
out.Pblock = Pblk;
out.block_fraction = Pblk/out.Pgel;
end



% =========================================================================
function It = current_spectrum(kx, cfg)
%CURRENT_SPECTRUM  Itilde(kx) = Int I(x) exp(+j kx x) dx   [A*m]
%
%   Feed at x = 0, open tip at x = L, unit reference current amplitude.
%
%   THE BRAID RETURN ARM (cfg.L_sheath, see MODEL_CONFIG).  The element's
%   current is carried on 0 <= x <= L exactly as before.  If L_sheath > 0 a
%   second segment on -L_sheath <= x <= 0 carries the return current on the
%   outside of the coax braid,
%       I(x) = I(0) * cos(pi x / (2 L_sheath)) ,
%   continuous at the feed and zero at the far end.  Truncating the current at
%   x = 0 instead is what put a point charge Q = I(0)/w there and made the
%   deposition singular.
%
%   The two segments are integrated with SEPARATE Simpson rules meeting at
%   x = 0.  That matters: I(x) has a kink there for every ke_case except the
%   quarter-wave baseline, and a single rule straddling a kink would lose the
%   fourth-order accuracy the element quadrature is relying on.  It also keeps
%   the element segment's samples IDENTICAL to the truncated code, so
%   L_sheath = 0 reproduces v1/v2 bit for bit.
%
%   Identical to the private function of the same name in layered_field.m.
if strcmpi(cfg.Idist,'hertzian')
    It = cfg.L*ones(size(kx));
    return
end
% ---- the radiating element, 0 <= x <= L (unchanged from v1/v2) ----------
M  = 4001;                                      % odd -> Simpson
xs = linspace(0, cfg.L, M);
switch lower(cfg.Idist)
    case 'sinusoid'          % Balanis, Antenna Theory 3rd ed., Eq. (4-56)
        Ix = sin(cfg.ke*(cfg.L - xs));
    case 'uniform'
        Ix = ones(size(xs));
    case 'triangular'
        Ix = 1 - xs/cfg.L;
    otherwise
        error('sar_volume:Idist','unknown current distribution');
end
It = seg_ft(xs, Ix, kx);

% ---- the braid return arm, -L_sheath <= x <= 0 --------------------------
% Ix(1) is I(0+), so the arm starts from exactly the value the element ends on
% and I(x) is continuous -- which is the whole point: a continuous I has no
% delta in dI/dx and therefore no point charge at the feed.
if isfield(cfg,'L_sheath') && cfg.L_sheath > 0
    Ms = 1001;                                  % odd -> Simpson
    xf = linspace(-cfg.L_sheath, 0, Ms);
    If = Ix(1)*cos(pi*xf/(2*cfg.L_sheath));
    It = It + seg_ft(xf, If, kx);
end
end

% =========================================================================
function It = seg_ft(xs, Ix, kx)
%SEG_FT  Simpson quadrature of Int Ix(x) exp(+j kx x) dx over one segment.
%   xs uniformly spaced, an odd number of samples.
M  = numel(xs);
hs = xs(2) - xs(1);
Wt = ones(1,M); Wt(2:2:end-1) = 4; Wt(3:2:end-2) = 2; Wt = Wt*hs/3;
It = (Ix(:).'.*Wt) * exp(1j*(xs(:))*kx(:).');
end
