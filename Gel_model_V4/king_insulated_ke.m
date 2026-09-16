function [ke, info] = king_insulated_ke(a, b, epsr_ins, epsr_medium, f)
%KING_INSULATED_KE  Wavenumber of the current on an insulated wire antenna in
%                   a dense (high-permittivity, lossy) medium.
%
%   ke = KING_INSULATED_KE(a, b, epsr_ins, epsr_medium, f)
%     a            [m]  conductor radius
%     b            [m]  outer radius of the insulation
%     epsr_ins     [-]  relative permittivity of the insulation (PTFE 2.1)
%     epsr_medium  [-]  complex relative permittivity of the surrounding
%                       medium, eps' - j eps'' (time convention exp(+j w t))
%     f            [Hz] frequency
%   ke is returned in the exp(+j w t) convention, Im(ke) <= 0, ready for
%   E_field_2d(..., 'ke', ke).  info carries lambda, the 1/e length and the
%   inputs.
%
%   King, Trembly & Strohbehn, "The electromagnetic field of an insulated
%   antenna in a conducting or dielectric medium", IEEE Trans. MTT-31 (1983)
%   574-583, Eq. (7); King & Smith, Antennas in Matter (1981), ch. 8:
%
%       k_L = k_2 [ 1 + H0(k_3 b) / ( k_3 b H1(k_3 b) ln(b/a) ) ]^(1/2)
%
%   k_2 = insulation wavenumber, k_3 = medium wavenumber, H = Hankel functions
%   of the first kind (their exp(-i w t) convention, Im k_3 >= 0).  Valid when
%   |k_3| >> |k_2| and |k_L b| << 1 -- both hold here (|k_3|/k_2 ~ 6,
%   |k_L b| ~ 0.06).  The formula is for an unbounded medium; with the wire
%   0.7 mm above a gel of lower permittivity the true value lies between the
%   water and gel answers -- which differ by only 4 % (127-30j vs 132-28j
%   for RG-178), so the choice of ambient hardly matters.

eps0 = 8.8541878128e-12;  mu0 = 1.25663706212e-6;
k0 = 2*pi*f*sqrt(mu0*eps0);
k2 = k0*sqrt(epsr_ins);
k3 = k0*sqrt(conj(epsr_medium));                 % exp(-i w t): eps' + i eps''
if imag(k3) < 0, k3 = conj(k3); end
x  = k3*b;
kL = k2*sqrt(1 + besselh(0, 1, x)./(x.*besselh(1, 1, x)*log(b/a)));
if imag(kL) < 0, kL = -kL; end                   % decaying wave, Im >= 0
ke = conj(kL);                                   % back to exp(+j w t)
info = struct('ke', ke, 'lambda', 2*pi/real(ke), 'L_1e', -1/imag(ke), ...
              'k2', k2, 'k3', conj(k3), 'a', a, 'b', b, 'epsr_ins', epsr_ins, ...
              'epsr_medium', epsr_medium, 'f', f);
end
