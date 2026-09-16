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

% ------------------------------------ THE RETURN CURRENT ON THE BRAID -----
% v1 and v2 prescribed the current on 0 <= x <= L and then simply STOPPED it.
% Under the quarter-wave hypothesis the current is at its MAXIMUM where it
% stops, so charge continuity dI/dx + j w rho = 0 puts a point charge
% Q = I_m/w at the feed and the deposition -- hence the temperature -- diverges
% there.  Both versions lived with it by BLANKING a 2 mm sphere about the feed
% and quoting the peak outside.  That is the grey semicircle in the GUI maps,
% and the photographs of the 50 W run (../50W CW 6 min/) show it is wrong: the
% gel heats smoothly there, the hottest zone sits a couple of millimetres along
% the element and a couple of millimetres down, and the blob reaches BACK past
% the junction under the braid.  No hole, no spike, and heating behind the
% feed that the old model could not produce at all.
%
% WHY THE POINT CHARGE IS AN ARTEFACT, AND WHAT REPLACES IT.  The current
% arriving at the junction does not stop and it does not pile up as charge: it
% returns.  There is no ground plane and no balun in this bench setup, so the
% return path is the OUTSIDE of the bare, immersed coax braid.  The structure
% is not a monopole over a ground plane at all -- it is an asymmetric dipole
% whose second arm is the braid.  Truncating the current at x = 0 is exactly
% the statement "the second arm does not exist", and the point charge is the
% price of saying it.  v2's own notes list "current on the coax braid, which is
% bare and immersed -- not modelled" as an unfixed cause of the remaining
% error; this is that item.
%
% So the braid carries an arm of length L_sheath behind the feed:
%
%     I(x) = sin(ke (L - x))                     0 <= x <= L    (unchanged)
%     I(x) = I(0) * cos(pi x / (2 L_sheath))
%                                       -L_sheath <= x <  0     (new)
%
% I is continuous at the feed by construction, so there is no delta in dI/dx
% and no point charge.  With L_sheath = L the arm is the mirror of the element
% and dI/dx vanishes at x = 0 from BOTH sides, i.e. the charge density is zero
% at the feed and largest at the two far ends -- which is what a dipole does.
% The cosine form is used rather than sin(ke(L_sheath + x)) because it stays
% bounded for every L_sheath; the two are identical when L_sheath = lam_e/4.
%
% CHOOSING IT.  L_sheath = L makes the structure symmetric and is the natural
% default with nothing measured.  Shorter means the sheath wave dies faster and
% concentrates charge at the junction (max |dI/dx| = I(0) pi/(2 L_sheath),
% which is 6x the element's own peak at 2 mm and equal to it at 12 mm); longer
% spreads it further back.  It is NOT derived -- the sheath wave is not solved
% for -- so sweep it.  The GUI exposes it.
%
% L_sheath = 0 reproduces v1/v2 EXACTLY, bit for bit: the arm segment is then
% empty and the element quadrature is untouched.  VERIFY_V2 test 2 relies on
% that, and the GUI keeps it as a preset so the two can be seen side by side.
C.L_sheath = 12e-3; % [m]   braid arm behind the feed.  0 = v1/v2 truncated.

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

% ------------------------------------------- BATH ABOVE THE GEL (v2) ------
% v1 put AIR above the gel and had the antenna radiate into it (epsr_air = 1
% against the gel's 74.3 -- an enormous contrast, which is what drove 99.6 %%
% of the delivered power into the gel).  The bench photographs of the 50 W
% run show the antenna FULLY IMMERSED, so layer 0 is a liquid whose
% permittivity is close to the gel's own, the interface nearly disappears
% electrically, and a large share of the power goes up into the bath instead.
% That is the change v2 exists to make.
C.upper         = 'water';   % 'water' | 'air'   ('air' reproduces v1 exactly)
C.bath_NaCl_pct = 0;         % [%% w/v] 0 = deionised/tap, 0.9 = saline.
                             %       NOT MEASURED -- both ends are run as a
                             %       bracket; see v2 README.
C.T_bath_C      = 25;        % [degC] bath temperature
C.Zw            = 40e-3;     % [m] depth of bath integrated for absorption
C.dzw           = 1e-3;      % [m] step for that integral

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

% v2: resolve the medium above the gel AFTER the overrides, so that both
% model_config('upper','air') and model_config('bath_NaCl_pct',0.9) work.
switch lower(C.upper)
    case 'air'
        C.epsr_upper = C.epsr_air;
    case 'water'
        C.epsr_upper = bath_medium(C.f, C.T_bath_C, C.bath_NaCl_pct);
    otherwise
        error('model_config:upper', ...
              'unknown upper medium ''%s'' (want ''water'' or ''air'')', C.upper);
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
