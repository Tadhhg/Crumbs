function [Exs, Eys, Ezs] = bath_spectral_E(zc, h, a0d, b0d, a0u, b0u, ...
                                           Ra, Rb, kz0, KX, KY, KRg, Z0tm, ...
                                           w, epsC0)
%BATH_SPECTRAL_E  Spectral E-field components in the bath (layer 0) at height
%                 zc >= 0, referenced the same way as the gel expressions.
%
%   Below the wire (0 <= zc < h) the field is the primary DOWN-going wave plus
%   the wave the slab reflects back UP.  Above the wire only up-going waves
%   remain -- the bath is a half-space, so nothing comes back down.
%
%   SIGN CONVENTION.  a0u/b0u are referenced at z = 0 with the up-going TM
%   flip already built in (note the +1j in b0u against the -1j in b0d), which
%   is the convention the upward-flux formula uses.  The gel expressions
%   instead keep both directions in one convention and apply the flip when
%   the fields are formed (Bd = bd - bu).  S_TM converts between the two.
%   S_TM = +1 is not a guess.  Continuity of tangential E across z = 0 forces
%   Z0tm*(1 - r01b) == Z1tm*(1 + r01b), and substituting
%   r01b = (Z0tm - Z1tm)/(Z0tm + Z1tm) makes both sides 2*Z0tm*Z1tm/(Z0tm+Z1tm)
%   identically -- true for any impedances, so it is the convention, not a fit.
%   The same choice simultaneously makes normal D continuous (Bs matches on
%   both sides, and Ez carries the 1/eps).  Both are checked numerically in
%   verify_v2.m test 1; S_TM = -1 fails them by ~1 %.
S_TM = +1;

if zc <= h
    ad =         a0d.*exp( 1j*kz0*zc);
    bd =         b0d.*exp( 1j*kz0*zc);
    au =     Ra.*a0d.*exp(-1j*kz0*zc);
    bu = S_TM*Rb.*b0d.*exp(-1j*kz0*zc);
else
    ad = 0;  bd = 0;
    au =      (a0u + Ra.*a0d).*exp(-1j*kz0*zc);
    bu = S_TM*(b0u + Rb.*b0d).*exp(-1j*kz0*zc);
end

Es = ad + au;   Bs = bd + bu;   Bd = bd - bu;
Exs =  Es.*(KY./KRg) + Z0tm.*Bd.*(KX./KRg);
Eys = -Es.*(KX./KRg) + Z0tm.*Bd.*(KY./KRg);
Ezs =  Bs.*KRg/(w*epsC0);
end
