function P = model_defaults()
%MODEL_DEFAULTS  Every constant of the V4 gel-heating model, in one struct,
%                set to the bench run that was photographed (50 W CW, 6 min,
%                ../50W CW 6 min/, 2026-08-25).  The GUI's "Constants &
%                boundaries" tab edits these fields; gel_heating_solve(P)
%                caches one solution per distinct P.
%
%   Sources: [T] Experimental_setup_table.xlsx (../Gel_model_V3)
%            [Ph] measured on the photographs -- ../Gel_model_V3/gel_heating_notes.md s5
%            [v2] ../v2_water_above/README.md (conditions on the day)
%            [I&D] Incropera & DeWitt, Fundamentals of Heat and Mass Transfer
%            [K] King, Trembly & Strohbehn 1983 / King & Smith 1981 (king_insulated_ke)
%   Units are SI throughout; temperatures in degC.

P = struct();

% ------------------------------------------------------------- ANTENNA --
P.f        = 2.45e9;      % [Hz]  microwave frequency                        [T]
P.L        = 12e-3;       % [m]   element length: block-width scale of the
                          %       photographs gives 12.1 mm (12.5-13 with
                          %       perspective).  The table's 27 mm is another
                          %       antenna.                                    [Ph]
P.a        = 0.155e-3;    % [m]   conductor radius.  RG-178-class cable
                          %       ASSUMED from the braid OD (~1.6 mm) and the
                          %       stranded look of the element               [Ph]
P.b_ins    = 0.42e-3;     % [m]   outer radius of the PTFE left on the inner
                          %       conductor (RG-178 dielectric 0.84 mm).  The
                          %       lossy medium starts here.  = a for a bare wire
P.epsr_ins = 2.1;         % [-]   PTFE
P.h_axis   = 1.1e-3;      % [m]   wire axis above the gel surface             [Ph]
P.L_sheath = 4e-3;        % [m]   return-current arm on the braid.  The
                          %       photographs show only a faint lobe under the
                          %       braid; 4 mm reproduces "faint", 27 mm made the
                          %       braid the hottest thing in the tank (V3 notes s6.1)
P.ke       = [];          % [1/m] current wavenumber on the element.  [] = King's
                          %       insulated-antenna value from a, b_ins,
                          %       epsr_ins and the water (127-30j, lambda 49 mm:
                          %       the element is a quarter wave).  A number
                          %       overrides it (V3 used the gel's 344-52j).   [K]
P.eta      = 0.2;         % [-]   delivered / generator setting ("runs at 20 %
                          %       of the set power").  Not measured.
P.P_ref    = 100;         % [W]   generator setting of the stored unit response

% ----------------------------------------------------------------- GEL --
P.rho      = 1000;        % [kg/m^3]  1 % Phytagel in saline, >98 % water
P.c        = 4180;        % [J/(kg K)]
P.k        = 0.60;        % [W/(m K)]
P.epsr_gel = 43.98;       % [-]   real permittivity                            [T]
P.epsi_gel = 13.5;        % [-]   loss factor (total: ionic + dipolar); the
                          %       table's 1.7 S/m would be 12.5                [T]

% --------------------------------------------- WATER ABOVE THE GEL --
P.T_bath0  = 25;          % [degC] water temperature the gel has equilibrated
                          %        with before the run (everything at 25 degC
                          %        on the photographed day)                    [v2]
P.Q_circ   = 5/60000;     % [m^3/s] circulation flow, 5 L/min                 [T]
P.A_flow   = 30e-3*30e-3; % [m^2]  flow cross-section over the gel, ASSUMED =
                          %        the gel face
P.nu_w     = 8.9e-7;      % [m^2/s] water at 25 degC (6.98e-7 at 37)         [I&D A.6]
P.k_w      = 0.607;       % [W/(m K)]
P.Pr_w     = 6.1;         % [-]
P.epsr_w   = 73 - 8.8j;   % [-]   water permittivity, for King's ke           [T]
P.h_top    = [];          % [W/(m^2 K)] [] = laminar flat plate
                          %       Nu = 0.664 Re^1/2 Pr^1/3 (~1400); a number
                          %       overrides; Inf pins the face at the water
                          %       temperature                                  [I&D ch. 7]

% ------------------------------------------ TANK WALLS: SIDES, BOTTOM --
P.T_room   = 25;          % [degC] room air
P.h_air    = 10;          % [W/(m^2 K)] free convection + radiation outside   [I&D ch. 9]
P.t_wall   = 5e-3;        % [m]   acrylic, ASSUMED from the photographs
P.k_wall   = 0.19;        % [W/(m K)] acrylic
P.h_wall   = [];          % [W/(m^2 K)] [] = h_air in series with the wall (~7.9);
                          %       a number overrides; 0 = adiabatic

% ------------------------------------------------------ BLOCK AND GRID --
P.x_min    = -6.5e-3;     % [m]   block edge behind the feed (feed at x = 0)  [Ph]
P.x_max    = 23.5e-3;     % [m]   block edge beyond the tip: 30 mm block      [Ph, T]
P.W_gel    = 30e-3;       % [m]   width across the antenna (y)                [T]
P.H_gel    = 30e-3;       % [m]   depth below the surface.  The photographs
                          %       show ~20 mm; irrelevant within 6 min         [T]
P.dh       = 0.5e-3;      % [m]   grid spacing.  Cost ~ dh^-5 (nodes x steps):
                          %       0.5 mm ~10 s, 0.25 mm ~5 min, 1 mm < 1 s
P.t_end    = 720;         % [s]   length of the stored response
P.dt_frame = 2;           % [s]   stored frame interval (interpolated between)
end
