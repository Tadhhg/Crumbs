function [Et, S] = E_field_2d(x, z, phase, P_set, varargin)
%E_FIELD_2D  Electric-field vector at a point of the gel plane at a given
%            phase of the 2.45 GHz cycle.  The 27 mm monopole lies parallel to
%            the gel surface; the gel is treated as ONE unbounded lossy
%            dielectric (Helmholtz equation with the complex permittivity of
%            Experimental_setup_table.xlsx).
%
%   Et = E_FIELD_2D(x, z, phase, P_set)
%     x, z    [m]   point(s) in the plane y = 0 that contains the antenna.
%                   Antenna frame: feed at x = 0, tip at x = L = 27 mm, wire
%                   axis at z = +h; gel surface z = 0, gel is z <= 0.
%     phase   [rad] omega*t, scalar or vector
%     P_set   [W]   generator setting; delivered power = eta*P_set (eta = 0.2)
%     Et.x, Et.z    [size(x) x numel(phase)]  E_x and E_z in V/m, i.e.
%                   E(x,z,t) = Re{ Etil(x,z) exp(j*phase) }.  E_y = 0 on y = 0.
%   [Et, S] = ...   S carries the phasor (S.Ex, S.Ez), amplitudes and phase
%                   angles, E_rms, E_peak, E_min (semi-axes of the polarisation
%                   ellipse), the heating rate S.q [W/m^3], and the constants.
%
%   Name/value overrides:  'eta' 0.2 | 'ke' (current wavenumber; default = the
%   medium's k_c) | 'L_sheath' 27e-3 (0 = truncated monopole) | 'epsr' |
%   'h' 0.5e-3 | 'a' 0.45e-3 | 'ds' 0.1e-3.
%
%   EQUATIONS (see E_field_2d_notes.md for sources and justification)
%   (1) Helmholtz, lossy dielectric:  grad^2 E + kc^2 E = 0,
%       kc = w sqrt(mu0 eps0 (eps' - j eps''))          [Ellingson Vol.2 s3.3-3.4]
%   (2) monopole current:  I(x') = Im sin(ke (L - x')), 0 <= x' <= L
%                                                       [Balanis AT3 Eq. 4-56]
%       return current on the braid, I(0) cos(pi x'/(2 L_sheath)), -L_sheath <= x' < 0
%                                                       [v2_water_above/README.md s1]
%   (3) field = superposition of x-directed current elements I(x') ds in the
%       unbounded medium                                [Balanis AT3 Eqs. 4-8a..4-10c,
%                                                        ../matlab/hertzian_dipole_field.m]
%   (4) scale:  P_del = int 1/2 w eps0 eps'' |E|^2 dV over the whole medium
%       (nothing is reflected, everything radiated is absorbed)
%                                                       [E_field_model_writeup.md s6.4]
%   (5) E(t) = Re{ Etil e^{j w t} }; |E(t)|^2 = 1/2 [ |Etil|^2 + |Q| cos(2wt + arg Q) ],
%       Q = Ex^2 + Ez^2                                 [E_field_model_writeup.md s10.1]
%
%   Time convention exp(+j w t); epsr = eps' - j eps'', Im(kc) <= 0.

% ------------------------------------------------------------ PARAMETERS --
p = struct( ...
    'f',        2.45e9, ...             % [Hz]  setup table
    'epsr',     43.98 - 1j*13.5, ...    % [-]   setup table (eps'' is the total
    ...                                 %       measured loss factor; 1.7 S/m
    ...                                 %       would be eps'' = 12.5, -8 %)
    'L',        27e-3, ...              % [m]   setup table, monopole length
    'eta',      0.2, ...                % [-]   delivered / set power (user:
    ...                                 %       "runs at 20 % of the set power")
    'h',        0.5e-3, ...             % [m]   wire axis above the gel surface
    ...                                 %       (../matlab/model_config.m)
    'a',        0.45e-3, ...            % [m]   conductor radius (model_config)
    'L_sheath', 27e-3, ...              % [m]   braid return arm, = L as in v2.1
    'ke',       [], ...                 % [1/m] current wavenumber; [] -> kc
    'ds',       0.1e-3, ...             % [m]   quadrature segment along the wire
    'dxg',      0.25e-3, ...            % [m]   x step of the power integral
    'Nr',       240, ...                % [-]   log-spaced rho points of the same
    'npad',     3);                     % [-]   its extent, in field 1/e lengths
for q1 = 1:2:numel(varargin), p.(varargin{q1}) = varargin{q1+1}; end

eps0 = 8.8541878128e-12;  mu0 = 1.25663706212e-6;      % CODATA 2018
w  = 2*pi*p.f;
kc = w*sqrt(mu0*eps0*p.epsr);   if imag(kc) > 0, kc = conj(kc); end   % (1)
if isempty(p.ke), p.ke = kc; end
sig_eff = w*eps0*(-imag(p.epsr));                     % [S/m] = w eps0 eps''
P_del   = p.eta*P_set;

% ------------------------------------------------- CURRENT ON THE WIRE --
xs  = (-p.L_sheath + p.ds/2) : p.ds : (p.L - p.ds/2);   % segment midpoints
Idl = wire_current(xs, p.ke, p.L, p.L_sheath)*p.ds;     % I_m = 1 A     (2)

% --------------------------------------- SCALE FROM DELIVERED POWER (4) --
P1 = dissipated_power(p, kc, w, sig_eff, xs, Idl);      % [W] at I_m = 1
Im = sqrt(P_del/P1);

% ------------------------------------------------- FIELD AT THE POINT(S) --
[Ex, Ez] = wire_field(x, z, xs, Im*Idl, kc, w, p.h);    % phasors       (3)

% ---------------------------------------------------- WAVEFORM IN TIME --
ph  = reshape(phase, 1, []);                                            % (5)
Et  = struct('x', real(Ex(:).*exp(1j*ph)), 'z', real(Ez(:).*exp(1j*ph)));
if ~isscalar(Ex)
    Et.x = reshape(Et.x, [size(Ex) numel(ph)]);
    Et.z = reshape(Et.z, [size(Ez) numel(ph)]);
end

if nargout > 1
    E2 = abs(Ex).^2 + abs(Ez).^2;   Q = Ex.^2 + Ez.^2;
    S = struct('Ex', Ex, 'Ez', Ez, ...
               'Ex_amp', abs(Ex), 'Ez_amp', abs(Ez), ...
               'Ex_phase', angle(Ex), 'Ez_phase', angle(Ez), ...
               'E_rms',  sqrt(E2/2), ...
               'E_peak', sqrt((E2 + abs(Q))/2), ...
               'E_min',  sqrt((E2 - abs(Q))/2), ...
               'Q', Q, 'q', 0.5*sig_eff*E2, ...            % [W/m^3] heating
               'kc', kc, 'ke', p.ke, 'lambda', 2*pi/real(kc), ...
               'alpha', -imag(kc), 'sig_eff', sig_eff, ...
               'Im', Im, 'I_feed', Im*sin(p.ke*p.L), 'P_del', P_del, ...
               'P1', P1, 'xs', xs, 'Idl', Idl, 'params', p);
end
end

% =========================================================================
function I = wire_current(xs, ke, L, Ls)
% Balanis AT3 Eq. (4-56) on the element; v2.1 cosine arm on the braid.
I  = zeros(size(xs));
on = xs >= 0 & xs <= L;
I(on) = sin(ke*(L - xs(on)));
if Ls > 0
    arm = xs < 0 & xs >= -Ls;
    I(arm) = sin(ke*L)*cos(pi*xs(arm)/(2*Ls));
end
end

function [Ex, Ez] = wire_field(x, z, xs, Idl, k, w, h)
% Sum of x-directed current elements Idl at (xs, 0, h), observed at (x, 0, z).
% Kernel is ../matlab/hertzian_dipole_field.m verbatim: exact field of a
% Hertzian dipole in an unbounded medium, E = -j w mu0 Idl [I + grad grad/k^2] G,
% G = exp(-jkR)/(4 pi R)  [Balanis AT3 Sec. 4.2 Eqs. (4-8a)-(4-10c)].
mu0 = 1.25663706212e-6;
Ex = zeros(size(x));  Ez = zeros(size(x));
Z  = z - h;
for s = 1:numel(xs)
    X  = x - xs(s);
    R  = sqrt(X.^2 + Z.^2);
    G  = exp(-1j*k*R)./(4*pi*R);
    C2 = -k^2./R.^2 + 3j*k./R.^3 + 3./R.^4;
    C1 = -(1j*k./R + 1./R.^2);
    pre = -1j*w*mu0*Idl(s);
    Ex = Ex + pre*(G + G.*(C2.*X.*X + C1)/k^2);
    Ez = Ez + pre*(    G.*(C2.*Z.*X)/k^2);
end
end

function P1 = dissipated_power(p, kc, w, sig_eff, xs, Idl)
% P = int 1/2 sig_eff |E|^2 dV over the whole medium outside the conductor,
% for I_m = 1 A.  The field is axisymmetric about the wire, so dV = 2 pi rho
% drho dx and the plane y = 0 below the wire (rho = h - z) samples it.
% rho is log-spaced from the conductor radius a because |E| ~ 1/rho near a
% line charge; the outer limits are ~6 power e-folding lengths.
persistent cache
key = sprintf('%.9g,', [real(kc) imag(kc) real(p.ke) imag(p.ke) ...
                        p.L p.L_sheath p.a p.h p.ds p.dxg p.Nr p.npad]);
if ~isempty(cache) && isKey(cache, key), P1 = cache(key); return; end
if isempty(cache), cache = containers.Map(); end

pad  = p.npad/(-imag(kc));   Rmax = p.a + pad;           % npad field 1/e lengths
xg   = (-p.L_sheath - pad) : p.dxg : (p.L + pad);
rg   = logspace(log10(p.a), log10(Rmax), p.Nr);
[XG, RG] = ndgrid(xg, rg);
[Ex, Ez] = wire_field(XG, p.h - RG, xs, Idl, kc, w, p.h);
qv = 0.5*sig_eff*(abs(Ex).^2 + abs(Ez).^2);              % [W/m^3]
P1 = trapz(xg, trapz(rg, qv.*RG*2*pi, 2), 1);
cache(key) = P1;
end
