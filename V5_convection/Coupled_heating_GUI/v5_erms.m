function [E, info] = v5_erms(P, x, rho)
%V5_ERMS  E_rms (V/m) of the insulated monopole at 1 W delivered, on the
%         axisymmetric grid (x along the element, rho = distance from its axis).
%
%   [E, info] = v5_erms(P, x, rho)     x, rho vectors (m); E is numel(rho) x numel(x)
%
%   Same physics as ../Insulated_monopole_Erms/insulated_monopole_Erms.m:
%   King's insulated-antenna wavenumber k_L, a standing wave on the element
%   (zero at the tip), the return charge spread over the exposed feed block,
%   and the thin-wire field E = -jwA - grad(phi) with the gel Green's function
%   and the charge fields scaled by 2 eps_g/(eps_w + eps_g).  The field only
%   depends on (x, rho), so the 3-D solver interpolates this table.
%   E scales as sqrt(P_delivered).

c0 = 299792458;  mu0 = 4e-7*pi;  eps0 = 8.8541878128e-12;
w  = 2*pi*P.f;   k0 = w/c0;
eps_w = P.eps_w1 - 1j*P.eps_w2;
eps_g = P.eps_g1 - 1j*P.eps_g2;
k_ins = k0*sqrt(P.eps_ins);  k_w = k0*sqrt(eps_w);  k_g = k0*sqrt(eps_g);
a = P.a;  b = P.b;  L = P.L;  l_f = P.l_f;

% 1. current wavenumber on the insulated element (King, Trembly & Strohbehn 1983)
k_L = k_ins*sqrt(1 + besselh(0,2,k_w*b)/(k_w*b*besselh(1,2,k_w*b)*log(b/a)));

% 2. feed current for 1 W delivered (open-ended insulated line)
C_ins = 2*pi*eps0*P.eps_ins/log(b/a);
Z_in  = -1j*(k_L/(w*C_ins))*cot(k_L*L);
I0    = sqrt(2*1/real(Z_in));

% 3. current and line charge on the element and the feed block
dxs = min(0.05e-3, L/200);
nb  = max(2, round(l_f/dxs));  ne = max(20, round(L/dxs));
xb  = ((1:nb) - 0.5)*l_f/nb - l_f;   xe = ((1:ne) - 0.5)*L/ne;
xs  = [xb xe];  wx = [repmat(l_f/nb, 1, nb) repmat(L/ne, 1, ne)];
I = [I0*(1 + xb/l_f),  I0*sin(k_L*(L - xe))/sin(k_L*L)];
q = [repmat((1j/w)*I0/l_f, 1, nb),  (1j/w)*(-I0*k_L*cos(k_L*(L - xe))/sin(k_L*L))];

% 4. field on the (x, rho) grid
tau = 2*eps_g/(eps_w + eps_g);
[X, R0] = meshgrid(x(:).', rho(:));
Ex = zeros(size(X));  Er = Ex;
cA = -1j*w*mu0/(4*pi);  cQ = -tau/(4*pi*eps0*eps_g);
for n = 1:numel(xs)
    DX = X - xs(n);  R = sqrt(DX.^2 + R0.^2);
    G  = exp(-1j*k_g*R)./R;
    dGR = -(1 + 1j*k_g*R).*G./R.^2;           % (dG/dR)/R
    Ex = Ex + (cA*I(n))*wx(n)*G + (cQ*q(n)*wx(n))*dGR.*DX;
    Er = Er + (cQ*q(n)*wx(n))*dGR.*R0;
end
E = sqrt((abs(Ex).^2 + abs(Er).^2)/2);

info = struct('k_L', k_L, 'Z_in', Z_in, 'I0', I0, 'lambda_L', 2*pi/real(k_L), ...
              'sigma', w*eps0*P.eps_g2);
end
