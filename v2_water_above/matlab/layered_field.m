function out = layered_field(cfg)
%LAYERED_FIELD  Exact time-harmonic E-field of a horizontal wire antenna lying
%               above a lossy dielectric slab, by plane-wave (Sommerfeld)
%               spectral synthesis.
%
%   GEOMETRY  (z up, gel surface at z = 0)
%       layer 0 :  z > 0        air        (source: x-directed wire at z = h)
%       layer 1 : -d < z < 0    lossy gel  (thickness d)
%       layer 2 :  z < -d       air (or backing)
%
%   METHOD
%   The x-directed current filament is written as a magnetic vector potential
%   A = x_hat*A_x with the Weyl plane-wave (Sommerfeld) representation
%       exp(-j k R)/(4 pi R) = 1/(8 pi^2) Int Int exp(-j(kx x + ky y + kz|z|))
%                                                 /(j kz) dkx dky
%   [W. C. Chew, "Waves and Fields in Inhomogeneous Media", IEEE Press 1995,
%    Sec. 2.2, Eq. (2.2.27); originally Weyl 1919 / Sommerfeld 1909].
%   Each spectral plane wave is split into TE (s) and TM (p) polarisations and
%   propagated through the three-layer stack with the transmission-line
%   analogy of
%   [K. A. Michalski and J. R. Mosig, "Multilayered media Green's functions in
%    integral equation formulations", IEEE Trans. Antennas Propagat., vol. 45,
%    no. 3, pp. 508-519, 1997, doi:10.1109/8.558666],
%   using interface (Fresnel) coefficients written in terms of the modal wave
%   impedances  Z^TE = w*mu0/kz  and  Z^TM = kz/(w*eps).
%   The spatial field is recovered by a 2-D inverse Fourier transform (FFT).
%
%   TIME CONVENTION  exp(+j w t);  lossy media have epsr = eps' - j*eps''.
%
%   cfg FIELDS
%     f      frequency [Hz]
%     epsr   [epsr0 epsr1 epsr2] complex relative permittivities
%     d      gel thickness [m]
%     h      height of wire axis above the gel surface [m], h >= 0
%     L      antenna length [m]
%     ke     effective wavenumber of the current wave on the antenna [rad/m]
%     Idist  'sinusoid' | 'uniform' | 'triangular' | 'hertzian'
%     Nfft   FFT points per side
%     dx     spatial sample [m]
%     z      vector of depths to evaluate (negative, inside the gel) [m]
%     Pdel   delivered (net) power [W] used to scale the result
%
%   OUTPUT struct: complex phasor field on the y = 0 plane
%   (Ex,Ey,Ez : numel(z) x Nfft), the power budget, and the scale factor.

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
    t = sqrt(k(n)^2 - KR2);
    t(imag(t)>0) = -t(imag(t)>0);               % enforce Im(kz) <= 0
    kz{n} = t;
end
kz0 = kz{1}; kz1 = kz{2};

% ---- modal wave impedances -----------------------------------------------
epsC = eps0*cfg.epsr(:).';                      % complex permittivities
Z0te = w*mu0./kz{1}; Z1te = w*mu0./kz{2}; Z2te = w*mu0./kz{3};
Z0tm = kz{1}/(w*epsC(1)); Z1tm = kz{2}/(w*epsC(2)); Z2tm = kz{3}/(w*epsC(3));

% ---- interface coefficients ----------------------------------------------
% TE amplitude a = E along s_hat  (transmission-line VOLTAGE analogue)
r01a = (Z1te - Z0te)./(Z1te + Z0te);   r12a = (Z2te - Z1te)./(Z2te + Z1te);
% TM amplitude b = H along s_hat  (transmission-line CURRENT analogue)
r01b = (Z0tm - Z1tm)./(Z0tm + Z1tm);   r12b = (Z1tm - Z2tm)./(Z1tm + Z2tm);
t01a = 1 + r01a;                       t01b = 1 + r01b;

d    = cfg.d;
tau2 = exp(-2j*kz1*d);                          % two-way phase across slab
Da_den = 1 + r01a.*r12a.*tau2;
Db_den = 1 + r01b.*r12b.*tau2;
Ra = (r01a + r12a.*tau2)./Da_den;               % generalised slab reflection
Rb = (r01b + r12b.*tau2)./Db_den;

% ---- source spectrum ------------------------------------------------------
It = current_spectrum(kv, cfg);                 % Itilde(kx) [A*m], 1 x N
IT = repmat(It(:), 1, N);                       % independent of ky

KRg = KR; KRg(KRg==0) = eps;                    % guard removable kx=ky=0 point

A0d = mu0*IT.*exp(-1j*kz0*cfg.h)./(2j*kz0);     % primary down-going A_x at z=0
a0d = -1j*w*A0d.*(KY./KRg);
b0d = -1j*A0d.*KX.*kz0./(mu0*KRg);
A0u = mu0*IT.*exp(+1j*kz0*cfg.h)./(2j*kz0);     % primary up-going A_x at z=0
a0u = -1j*w*A0u.*(KY./KRg);
b0u = +1j*A0u.*KX.*kz0./(mu0*KRg);

% ---- amplitudes inside the gel (layer 1), referenced to z = 0 ------------
Da = t01a.*a0d./Da_den;
Db = t01b.*b0d./Db_den;

norm2d = 1/(N*dx*N*dy);        % ifft & Parseval normalisation (identical)

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
clear aU bU Ex0 Ey0 Ez0 a0u b0u a0d b0d A0d A0u IT Ra Rb Z0te Z0tm kz0 KR KR2

% ---- power transmitted downward out of the slab --------------------------
a2 = (1+r12a).*Da.*exp(-1j*kz1*d);
b2 = (1+r12b).*Db.*exp(-1j*kz1*d);
Pdn = 0.5*real(sum(sum(abs(a2).^2./conj(Z2te) + Z2tm.*abs(b2).^2)))*norm2d;
clear a2 b2 Z2te Z2tm Z1te r01a r01b t01a t01b Da_den Db_den tau2

% ================= evaluate on y = 0 and accumulate dissipated power =======
zz = cfg.z(:);   nz = numel(zz);
Ex = zeros(nz,N); Ey = zeros(nz,N); Ez = zeros(nz,N);
sigma1 = w*eps0*(-imag(cfg.epsr(2)));           % effective conductivity [S/m]
dz = abs(zz(2)-zz(1));
Pgel = 0;
iy0 = N/2 + 1;                                  % index of y = 0
izmid = round(nz/2);
slice = [];

for iz = 1:nz
    zc = zz(iz);
    ad = Da.*exp( 1j*kz1*zc);                   % down-going at depth zc
    au = Da.*r12a.*exp(-1j*kz1*(2*d + zc));     % up-going (overflow-safe form)
    bd = Db.*exp( 1j*kz1*zc);
    bu = Db.*r12b.*exp(-1j*kz1*(2*d + zc));

    Es = ad + au;                 % TE E-amplitude along s_hat
    Bs = bd + bu;                 % TM H-amplitude along s_hat
    Bd = bd - bu;

    Exs =  Es.*(KY./KRg) + Z1tm.*Bd.*(KX./KRg);
    Eys = -Es.*(KX./KRg) + Z1tm.*Bd.*(KY./KRg);
    Ezs =  Bs.*KRg/(w*epsC(2));

    ex = fftshift(fft2(Exs))*norm2d;
    ey = fftshift(fft2(Eys))*norm2d;
    ez = fftshift(fft2(Ezs))*norm2d;

    Pgel = Pgel + 0.5*sigma1*sum(abs(ex(:)).^2 + abs(ey(:)).^2 + ...
                                 abs(ez(:)).^2)*dx*dy*dz;

    Ex(iz,:) = ex(:,iy0).';  Ey(iz,:) = ey(:,iy0).';  Ez(iz,:) = ez(:,iy0).';
    if iz == izmid, slice = struct('Ex',ex,'Ey',ey,'Ez',ez,'z',zc); end
end

Ptot  = Pgel + Pwat + Pup + Pdn;                % v2: the bath is in the budget
scale = sqrt(cfg.Pdel/Ptot);                    % amplitude scale factor

out = struct('x',x,'y',y,'z',zz,'Ex',Ex*scale,'Ey',Ey*scale,'Ez',Ez*scale, ...
    'scale',scale,'Pgel',Pgel,'Pup',Pup,'Pdn',Pdn,'Pwat',Pwat, ...
    'frac_gel',Pgel/Ptot,'Ptot',Ptot, ...
    'k',k,'sigma1',sigma1,'cfg',cfg);
if ~isempty(slice)
    out.slice = struct('Ex',slice.Ex*scale,'Ey',slice.Ey*scale, ...
                       'Ez',slice.Ez*scale,'z',slice.z);
end
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
%   quarter-wave baseline, and one rule straddling a kink would lose the
%   fourth-order accuracy the element quadrature relies on.  It also leaves the
%   element segment's samples IDENTICAL to the truncated code, so L_sheath = 0
%   reproduces v1/v2 bit for bit.
if strcmpi(cfg.Idist,'hertzian')
    It = cfg.L*ones(size(kx));                  % point dipole, moment I*L
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
        error('layered_field:Idist','unknown current distribution');
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
