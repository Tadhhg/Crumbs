function [epsr, info] = bath_medium(f, T_C, NaCl_pct)
%BATH_MEDIUM  Complex relative permittivity of the liquid the antenna is
%             immersed in -- the bath lying ON TOP of the gel.
%
%   [epsr, info] = BATH_MEDIUM(f, T_C, NaCl_pct)
%
%   v1 of this project put AIR above the gel (model_config: epsr_air = 1) and
%   the antenna radiated into it.  The bench photographs show the antenna
%   fully immersed in the bath instead, so layer 0 is a liquid with a
%   permittivity close to the gel's own -- which is why v2 exists.
%
%   NaCl_pct is w/v %.  Two cases bracket the bench:
%     0     deionised / tap water   eps'' small, the bath barely absorbs
%     0.9   physiological saline    eps'' ~ the gel's, the bath absorbs as
%                                   hard as the gel does
%   The salinity of the actual bath is NOT known (see v2 README), so both are
%   run and reported as a bracket rather than one being presented as the
%   answer.
%
%   Uses the same Stogryn (1971) saline model as the gel, so the bath and the
%   gel are described by one consistent model rather than two.
%
%   See also SALINE_STOGRYN, MODEL_CONFIG.

if nargin < 3, NaCl_pct = 0; end
N = NaCl_pct*10/58.443;                 % [mol/L] normality, M(NaCl)=58.443
[epsr, gi] = saline_stogryn(f, T_C, N);
info = struct('NaCl_pct',NaCl_pct, 'N_molL',N, 'T_C',T_C, ...
              'sigma_eff_Spm',gi.sigma_eff_Spm, 'eps_real',real(epsr), ...
              'source','Stogryn (1971) saline model, as used for the gel');
end
