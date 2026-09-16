function [epsr, sig_eff] = tissue_cole_cole(f, tissue)
%TISSUE_COLE_COLE  Complex relative permittivity of a body tissue from the
%   four-pole Cole-Cole parametric model.
%
%       eps_r(w) = eps_inf + SUM_n  d_eps_n / (1 + (j w tau_n)^(1-alpha_n))
%                          + sigma_i / (j w eps0)
%
%   MODEL / SOURCE
%     S. Gabriel, R. W. Lau and C. Gabriel, "The dielectric properties of
%     biological tissues: III. Parametric models for the dielectric spectrum
%     of tissues", Phys. Med. Biol., vol. 41, no. 11, pp. 2271-2293, 1996,
%     doi:10.1088/0031-9155/41/11/003; parameter tables in
%     C. Gabriel, "Compilation of the Dielectric Properties of Body Tissues at
%     RF and Microwave Frequencies", Report AL/OE-TR-1996-0037, Brooks AFB,
%     1996.  The MUSCLE parameter column used below is transcribed from
%     Table 1 of the open-access reproduction in
%     "An Analytical Framework for Frequency-Dependent Electromagnetic Power
%     Absorption in Biological Tissues", arXiv:2604.06225, Table 1
%     (column "Muscle"), which cites Gabriel et al. (1996) as its source.
%
%   VERIFICATION
%     At 2.45 GHz this returns eps' = 52.73 and sigma_eff = 1.74 S/m, which is
%     the value the IT'IS Foundation tissue database serves for muscle -- an
%     independent check that the parameters have been transcribed correctly.
%
%   Time convention exp(+j w t): epsr = eps' - j*eps''.

if nargin < 2, tissue = 'muscle'; end
eps0 = 8.8541878128e-12;
w = 2*pi*f;

switch lower(tissue)
    case 'muscle'
        einf = 4.0;   sig = 2.0e-1;
        de  = [50.0    7000      1.2e6    2.5e7  ];
        tau = [7.23e-12 353.68e-9 318.31e-6 2.274e-3];
        alf = [0.10    0.10      0.10     0.00   ];
    case 'liver'
        einf = 4.0;   sig = 2.0e-2;
        de  = [39.0    6000      5.0e4    3.0e7  ];
        tau = [8.84e-12 530.52e-9 22.74e-6 15.915e-3];
        alf = [0.10    0.20      0.20     0.05   ];
    case 'fat'          % fat, not infiltrated
        einf = 2.5;   sig = 1.0e-2;
        de  = [3.0     15        3.3e4    1.0e7  ];
        tau = [7.96e-12 15.92e-9 159.15e-6 7.958e-3];
        alf = [0.20    0.10      0.05     0.01   ];
    otherwise
        error('tissue_cole_cole:tissue','no parameters stored for "%s"', tissue);
end

epsr = einf + sig./(1j*w*eps0);
for n = 1:4
    epsr = epsr + de(n)./(1 + (1j*w*tau(n)).^(1-alf(n)));
end
sig_eff = -imag(epsr).*w*eps0;
end
