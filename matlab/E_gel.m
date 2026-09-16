function [Emag, S] = E_gel(x, z, t)
%E_GEL  THE DELIVERABLE.  Electric-field magnitude at a point in the 2-D gel
%       plane as a function of time.
%
%   Emag = E_GEL(x, z, t)   returns |E(x,z,t)| in V/m.
%   [Emag, S] = E_GEL(...)  also returns a struct with the full description.
%
%   INPUTS  (SI units -- metres and seconds)
%     x   along-antenna coordinate.  The 12 mm monopole occupies 0 <= x <= 12 mm
%         at height h above the surface; x = 0 is the feed, x = 12 mm the tip.
%     z   height.  Gel surface z = 0, gel floor z = -30 mm, so z is NEGATIVE
%         inside the gel.
%     t   time.  Scalar or array.  May be [] to get only the envelope.
%   x, z, t are broadcast against each other.
%
%   THE FUNCTION
%   The excitation is time-harmonic at f = 2.45 GHz, so the physical field is
%       E(x,z,t) = Re{ Etil(x,z) exp(j w t) },   w = 2 pi f
%   where Etil = (Ex, Ey, Ez) is the complex phasor computed by LAYERED_FIELD.
%   Writing Etil_i = |Etil_i| exp(j phi_i), each Cartesian component is a pure
%   cosine, and the MAGNITUDE follows in closed form:
%
%       |E(x,z,t)|^2 = SUM_i |Etil_i|^2 cos^2(w t + phi_i)
%                    = 1/2 [ |Etil|^2 + Re{ Q exp(j 2 w t) } ]
%
%   with     |Etil|^2 = Etil . conj(Etil)   (a real number)
%            Q        = Etil . Etil         (NO conjugate -- complex)
%
%   so the deliverable is
%
%       |E(x,z,t)| = sqrt( 1/2 [ |Etil(x,z)|^2 + |Q(x,z)| cos(2 w t + arg Q) ] )
%
%   The magnitude therefore oscillates at TWICE the drive frequency between
%   the semi-axes of the local polarisation ellipse:
%       E_peak = sqrt((|Etil|^2 + |Q|)/2)      (major semi-axis)
%       E_min  = sqrt((|Etil|^2 - |Q|)/2)      (minor semi-axis)
%       E_rms  = |Etil|/sqrt(2)                (time RMS -- the quantity that
%                                               sets heating and SAR)
%   Linear polarisation gives |Q| = |Etil|^2 and E_min = 0; circular gives
%   Q = 0 and a magnitude that is constant in time.
%
%   POWER SCALING
%   The stored solution is normalised to P_delivered = 150 W
%   (200 W forward - 50 W maximum reflected).  Because Maxwell's equations are
%   linear and the medium properties are taken as power-independent,
%       |E| scales as sqrt(P_delivered).
%   For 200 W delivered multiply the result by sqrt(200/150) = 1.155.
%
%   VALIDITY:  -30 mm <= z <= 0 and |x| <= 128 mm, on the plane y = 0.
%   Accuracy is limited by the assumptions listed in the write-up, Sec. 8 --
%   principally the current distribution on the monopole (about +-10 % on |E|)
%   and the gel permittivity (about +-15 %).

persistent M
if isempty(M)
    here = fileparts(mfilename('fullpath'));
    f = fullfile(fileparts(here),'results','field_model.mat');
    if ~exist(f,'file')
        error('E_gel:noModel',['results/field_model.mat not found. ' ...
              'Run  run_model  first.']);
    end
    M = load(f);
end

% -- interpolate the complex phasor components onto the requested points ----
[Xq, Zq] = deal(x, z);
Ex = interp2(M.x, M.z, M.Ex, Xq, Zq, 'linear', 0);
Ey = interp2(M.x, M.z, M.Ey, Xq, Zq, 'linear', 0);
Ez = interp2(M.x, M.z, M.Ez, Xq, Zq, 'linear', 0);

E2 = abs(Ex).^2 + abs(Ey).^2 + abs(Ez).^2;      % |Etil|^2
Q  = Ex.^2 + Ey.^2 + Ez.^2;                     % complex

w = M.omega;
if nargin < 3 || isempty(t)
    Emag = sqrt(E2/2);                          % default: return the RMS value
else
    Emag = sqrt(max(0.5*(E2 + real(Q.*exp(2j*w*t))), 0));
end

if nargout > 1
    S = struct();
    S.Ex_phasor = Ex;  S.Ey_phasor = Ey;  S.Ez_phasor = Ez;
    S.E_rms   = sqrt(E2/2);
    S.E_peak  = sqrt((E2 + abs(Q))/2);
    S.E_min   = sqrt(max(E2 - abs(Q),0)/2);
    S.axial_ratio = S.E_min./S.E_peak;
    S.Q       = Q;
    S.phase_of_Q = angle(Q);
    S.omega   = w;
    S.f       = w/(2*pi);
    S.SAR     = 0.5*M.sigma*E2/1000;             % W/kg, rho = 1000 kg/m^3
    S.Pv      = 0.5*M.sigma*E2;                  % W/m^3
    S.sigma   = M.sigma;
    S.Pdel    = M.C.Pdel;
    S.instantaneous = @(tt) [real(Ex.*exp(1j*w*tt)), real(Ey.*exp(1j*w*tt)), ...
                             real(Ez.*exp(1j*w*tt))];
end
end
