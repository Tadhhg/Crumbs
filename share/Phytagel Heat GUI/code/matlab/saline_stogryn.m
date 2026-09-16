function [epsr, info] = saline_stogryn(f, T, N)
%SALINE_STOGRYN  Complex relative permittivity of an aqueous NaCl solution.
%
%   [epsr, info] = SALINE_STOGRYN(f, T, N) returns the complex relative
%   permittivity  epsr = eps' - 1j*eps''  of a NaCl solution.
%
%   INPUTS
%     f  : frequency [Hz]
%     T  : temperature [degC]           (valid 0-40 degC)
%     N  : normality [mol/L]            (valid 0-3 N)
%
%   MODEL / SOURCE
%     Single-Debye relaxation plus an ionic-conduction loss term, with the
%     regression coefficients of:
%       A. Stogryn, "Equations for Calculating the Dielectric Constant of
%       Saline Water", IEEE Trans. Microwave Theory Tech., vol. MTT-19,
%       no. 8, pp. 733-736, 1971.  doi:10.1109/TMTT.1971.1127617
%     The coefficients are transcribed from the open-access reproduction in
%       D. H. Gadani, V. A. Rana, S. P. Bhatnagar, A. N. Prajapati, A. D. Vyas,
%       "Effect of salinity on the dielectric properties of water",
%       Indian J. Pure & Applied Physics, vol. 50, pp. 405-410, June 2012.
%       Equations (1)-(10) of that paper.  (Open access: nopr.niscpr.res.in)
%     The high-frequency limit eps_inf = 4.9 is the value adopted by
%       L. A. Klein and C. T. Swift, "An Improved Model for the Dielectric
%       Constant of Sea Water at Microwave Frequencies", IEEE Trans. Antennas
%       Propagat., vol. AP-25, no. 1, pp. 104-111, Jan. 1977
%       ("best fit to available high frequency data is achieved if eps_inf is
%        assumed a constant value of 4.9").
%
%   SIGN CONVENTION
%     exp(+1j*omega*t) time convention, so a lossy medium has
%     epsr = eps' - 1j*eps''  with eps'' >= 0.

eps0 = 8.8541878128e-12;         % [F/m] CODATA 2018 vacuum permittivity
eps_inf = 4.9;                   % Klein & Swift 1977

% --- static permittivity of pure water, Gadani eq.(6) [= Stogryn 1971] ------
es0 = 87.74 - 0.40008*T + 9.398e-4*T.^2 + 1.410e-6*T.^3;

% --- salt correction to static permittivity, Gadani eq.(4) -----------------
aN  = 1.000 - 0.2551*N + 5.151e-2*N.^2 - 6.889e-3*N.^3;

% --- static permittivity of the solution, Gadani eq.(2) --------------------
es  = es0 .* aN;

% --- 2*pi*tau for pure water [s], Gadani eq.(7) ----------------------------
twopitau0 = 1.1109e-10 - 3.824e-12*T + 6.938e-14*T.^2 - 5.096e-16*T.^3;

% --- salt correction to relaxation time, Gadani eq.(5) ---------------------
bNT = 1.000 + 0.1463e-2*N.*T - 0.04896*N - 0.02967*N.^2 + 5.644e-3*N.^3;

% --- relaxation time of the solution, Gadani eq.(3) ------------------------
twopitau = twopitau0 .* bNT;
tau = twopitau/(2*pi);

% --- ionic conductivity at 25 degC, Gadani eq.(10) -------------------------
sig25 = N .* (10.394 - 2.3776*N + 0.68258*N.^2 - 0.13538*N.^3 + 1.0086e-2*N.^4);

% --- temperature correction of conductivity, Gadani eq.(9) -----------------
D   = 25 - T;
sig = sig25 .* (1.000 - 1.962e-2*D + 8.08e-5*D.^2 ...
      - D.*N.*(3.020e-5 + 3.922e-5*D + N.*(1.721e-5 - 6.584e-6*D)));

% --- assemble complex permittivity, Gadani eq.(1) --------------------------
w = 2*pi*f;
epsp  = eps_inf + (es - eps_inf)./(1 + (w.*tau).^2);
epspp = (es - eps_inf).*(w.*tau)./(1 + (w.*tau).^2) + sig./(w*eps0);
epsr  = epsp - 1j*epspp;

info = struct('eps_s',es,'eps_inf',eps_inf,'tau_s',tau,'sigma_ionic_Spm',sig, ...
              'eps_real',epsp,'eps_imag',epspp, ...
              'sigma_eff_Spm', w*eps0.*epspp, ...
              'loss_tangent', epspp./epsp);
end
