function [Ex,Ey,Ez] = hertzian_dipole_field(f, epsr, Idl, X, Y, Z)
%HERTZIAN_DIPOLE_FIELD  Exact E-field of an x-directed Hertzian dipole at the
%   origin of an UNBOUNDED homogeneous (possibly lossy) medium.
%
%   Used only as an independent analytical reference to validate LAYERED_FIELD.
%
%   E(r) = -j*w*mu * I*dl * [ I_dyad + grad grad / k^2 ] G(R) . x_hat
%   with the scalar Green function G(R) = exp(-j k R)/(4 pi R).
%   [C. A. Balanis, "Advanced Engineering Electromagnetics", 2nd ed., Wiley
%    2012, Ch. 6 (dyadic Green functions); equivalently Balanis, "Antenna
%    Theory", 3rd ed., Sec. 4.2, Eqs. (4-8a)-(4-10c) rotated to x_hat.]
%
%   Time convention exp(+j w t); epsr = eps' - j*eps''.

eps0 = 8.8541878128e-12;
mu0  = 1.25663706212e-6;
w    = 2*pi*f;
k    = w*sqrt(mu0*eps0*epsr);
if imag(k) > 0, k = conj(k); end

R  = sqrt(X.^2 + Y.^2 + Z.^2);
G  = exp(-1j*k*R)./(4*pi*R);

% second derivatives of G:
%   d2G/dxi dxj = G[ (-k^2/R^2 + 3jk/R^3 + 3/R^4) xi xj - (jk/R + 1/R^2) dij ]
C2 = -k^2./R.^2 + 3j*k./R.^3 + 3./R.^4;
C1 = -(1j*k./R + 1./R.^2);

Gxx = G.*(C2.*X.*X + C1);
Gyx = G.*(C2.*Y.*X);
Gzx = G.*(C2.*Z.*X);

pre = -1j*w*mu0*Idl;
Ex = pre*(G + Gxx/k^2);
Ey = pre*(    Gyx/k^2);
Ez = pre*(    Gzx/k^2);
end
