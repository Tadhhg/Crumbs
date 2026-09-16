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

% ---- power radiated upward into air --------------------------------------
ztop = cfg.h + 5e-3;
aU = (a0u + Ra.*a0d).*exp(-1j*kz0*ztop);
bU = (b0u + Rb.*b0d).*exp(-1j*kz0*ztop);
Pup = 0.5*real(sum(sum(abs(aU).^2./conj(Z0te) + Z0tm.*abs(bU).^2)))*norm2d;
clear aU bU a0u b0u a0d b0d A0d A0u IT Ra Rb Z0te Z0tm kz0 KR KR2

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

Ptot  = Pgel + Pup + Pdn;
scale = sqrt(cfg.Pdel/Ptot);                    % amplitude scale factor

out = struct();
out.x = x(ix);  out.y = y(iy);  out.z = zn;
out.q = 0.5*sigma1*E2blk*scale^2;               % [W/m^3] at Pdel
out.Ex = Ex*scale;  out.Ey = Ey*scale;  out.Ez = Ez*scale;
out.x_full = x;
out.scale = scale;  out.sigma1 = sigma1;  out.k = k;
out.Pgel = Pgel*scale^2;  out.Pup = Pup*scale^2;  out.Pdn = Pdn*scale^2;
out.Ptot = Ptot;  out.cfg = cfg;

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
%   Identical to the private function of the same name in layered_field.m.
if strcmpi(cfg.Idist,'hertzian')
    It = cfg.L*ones(size(kx));
    return
end
M  = 4001;
xs = linspace(0, cfg.L, M);
hs = xs(2)-xs(1);
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
Wt = ones(1,M); Wt(2:2:end-1) = 4; Wt(3:2:end-2) = 2; Wt = Wt*hs/3;
It = (Ix.*Wt) * exp(1j*(xs(:))*kx(:).');
end
