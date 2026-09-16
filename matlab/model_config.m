function C = model_config(varargin)
%MODEL_CONFIG  Central definition of every physical constant and geometric
%              parameter used by the 2.45 GHz Phytagel field model.
%
%   Every entry carries the source it was taken from.  Override any field with
%   name/value pairs, e.g.  model_config('Pdel',200,'ke_case','bare_iface').
%
%   COORDINATES
%     origin = feed point of the monopole, projected onto the gel surface
%     x      = along the antenna axis; the 12 mm element spans 0 <= x <= L
%     y      = transverse horizontal;  the modelled plane is y = 0
%     z      = height; gel surface at z = 0, gel occupies -d <= z <= 0

c0 = 299792458;                 % [m/s] exact, SI definition of the metre

C = struct();

% ---------------------------------------------------------------- SOURCE ---
C.f     = 2.45e9;   % [Hz]  ISM band; stated by the user
C.Pfwd  = 200;      % [W]   generator maximum forward power (user)
C.Pref  = 50;       % [W]   generator maximum reflected power (user)
C.Pdel  = C.Pfwd - C.Pref;   % [W] net delivered power, worst case = 150 W

% ------------------------------------------------------------- GEOMETRY ---
C.L  = 12e-3;       % [m]   monopole length (user)
C.d  = 30e-3;       % [m]   gel depth (user)
C.h  = 0.5e-3;      % [m]   height of the wire AXIS above the gel surface.
                    %       Default = radius of a typical 0.86/2.98 mm coax
                    %       inner conductor + catheter wall, i.e. the wire is
                    %       resting on the surface.  Swept in run_sensitivity.
C.a  = 0.45e-3;     % [m]   conductor radius (assumed; only used for reporting)

% ------------------------------------------------- MEDIUM (baseline gel) ---
% Baseline: 1 %% w/v Phytagel in 0.9 %% w/v NaCl saline at 25 degC.
% At 1 %% w/v the gellan network occupies <1 %% by volume, so the complex
% permittivity is that of the saline to within the uncertainty of the model
% (see writeup Sec. 4).  Saline permittivity from the Stogryn (1971) model.
C.T_C   = 25;       % [degC] assumed bench temperature
C.NaCl_pct = 0.9;   % [%% w/v]
C.N     = C.NaCl_pct*10/58.443;   % [mol/L] normality; M(NaCl)=58.443 g/mol
                                  % (IUPAC standard atomic weights 2021)

C.epsr_air = 1 - 0j;              % dry air at 1 atm (eps_r = 1.00059, taken
                                  % as 1; error < 0.06 %%)
C.epsr_gel = saline_stogryn(C.f, C.T_C, C.N);
C.backing  = 'air';               % what lies below the gel: 'air' | 'pec'

% ------------------------------------------------- ANTENNA CURRENT WAVE ---
% Effective wavenumber k_e of the current wave on the monopole.  This is the
% single largest modelling uncertainty; three defensible hypotheses are kept
% and compared in run_sensitivity.m (writeup Sec. 6).
C.k0 = 2*pi*C.f/c0;                                    % [rad/m] free space
% H1 'quarterwave' : the 12 mm element is a resonant quarter wave, L = lam_e/4
% H2 'bare_iface'  : bare wire on the air/gel interface, eps_eff = (1+eps')/2
% H3 'bare_air'    : the current wave sees air only
C.ke_cases = struct('quarterwave', pi/(2*C.L), ...
                    'bare_iface',  C.k0*sqrt((1+real(C.epsr_gel))/2), ...
                    'bare_air',    C.k0);
C.ke_case = 'quarterwave';                             % baseline hypothesis
C.ke      = C.ke_cases.(C.ke_case);
C.Idist   = 'sinusoid';                                % Balanis Eq. (4-56)

% ------------------------------------------------------ NUMERICAL GRID ----
C.Nfft = 1024;              % FFT points per transverse direction
C.dx   = 0.25e-3;           % [m] spatial sample -> window +-128 mm
C.dz   = 0.25e-3;           % [m] depth sample
C.z    = -(C.dz/2 : C.dz : C.d - C.dz/2);   % cell-centred depths in the gel

% ------------------------------------------------------------ OVERRIDES ---
for q = 1:2:numel(varargin)
    C.(varargin{q}) = varargin{q+1};
end
if any(strcmp(varargin(1:2:end),'ke_case'))
    C.ke = C.ke_cases.(C.ke_case);
end
if any(strcmp(varargin(1:2:end),'T_C')) || any(strcmp(varargin(1:2:end),'NaCl_pct'))
    C.N = C.NaCl_pct*10/58.443;
    C.epsr_gel = saline_stogryn(C.f, C.T_C, C.N);
end
if strcmpi(C.backing,'pec')
    C.epsr_bot = -1e7j;      % near-PEC backing (large loss -> |r12| -> 1)
else
    C.epsr_bot = C.epsr_air;
end
C.epsr = [C.epsr_air, C.epsr_gel, C.epsr_bot];
end
