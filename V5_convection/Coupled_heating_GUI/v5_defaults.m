function P = v5_defaults(preset)
%V5_DEFAULTS  Every constant of the V5 coupled heating model, in one struct.
%
%   P = v5_defaults()          % 'table': ../Experimental_setup_table.xlsx
%   P = v5_defaults('photo')   % the photographed run, ../../50W CW 6 min
%
%   The equations are in ../Catheter_heating/Convection model.docx (sections 1-6):
%       rho c dT/dt = k lap(T) + Q_d,   Q_d = sigma E_rms^2 (P_del / P_ref)
%       -k dT/dz = h_s (T - T_diffusion)  on the gel surface,
%       T_diffusion = T_b + alpha(x) (T_cath(t) - T_b),  alpha = min(1, beta/(d sqrt f)) on 0 <= x <= L_c
%       T_cath: lumped, m c dT/dt = P_c - k_c (T - T_b),  P_c = gamma P_set P_dropoff
%   E_rms is the insulated-monopole field (../Insulated_monopole_Erms).
%
%   Fields are SI; temperatures in degC.  P.group says, for the GUI, which
%   fields change the stored solution ('solve') and which only rescale it
%   ('instant').  Sources: [T] setup table, [Ph] the photographs
%   (../../Gel_model_V3/gel_heating_notes.md s5), [Doc] Convection model.docx
%   s6, [V4] ../../Gel_model_V4/model_defaults.m, [K] King et al. 1983.

if nargin < 1, preset = 'table'; end
P = struct();
P.preset = preset;

% ------------------------------------------------ STARTING CONDITIONS -----
P.P_set     = 60;         % [W]    generator input power                      [T]
P.P_refl    = 0;          % [W]    reflected power                            [T]
P.P_dropoff = 0.7;        % [-]    fraction lost before the antenna (0.7 -> 30 % reaches it)  [T]
P.t_on      = 360;        % [s]    treatment time: power ends here            [Doc: 6 min run]

% ------------------------------------------- POWER ENVELOPE (instant) -----
P.duty      = 1;          % [-]    duty cycle within a burst (1 = CW)
P.T_pulse   = 10;         % [s]    pulse period (on for duty*T_pulse, then off)
P.t_burst   = 60;         % [s]    burst length: pulsing for t_burst, then ...
P.t_rest    = 0;          % [s]    ... a rest of t_rest (0 = no rests: one burst)
P.T_b       = 37;         % [degC] circulating water temperature              [T]
P.T_0       = 37;         % [degC] gel temperature at t = 0 (equilibrated)    [T]
P.T_room    = 25;         % [degC] room air outside the tank walls            [Doc]
P.flow      = 5;          % [L/min] circulation flow, f in alpha              [T]

% --------------------------------------------------- ANTENNA (solve) ------
P.f         = 2.45e9;     % [Hz]  microwave frequency                         [T]
P.L         = 27e-3;      % [m]   element length, feed to tip                 [T]
P.b         = 0.9e-3;     % [m]   outer radius of the insulated element =
                          %       the table's 1.8 mm diameter (read as the
                          %       surface the field leaves from, r0)          [T]
P.a         = 0.27e-3;    % [m]   conductor radius.  ASSUMED: a 50 ohm PTFE
                          %       line has b/a = exp(50 sqrt(2.1)/60) = 3.35
P.eps_ins   = 2.1;        % [-]   PTFE
P.l_f       = 1.5e-3;     % [m]   exposed metal feed block behind the feed    [Ph]
P.r_cable   = 1.5e-3;     % [m]   coax outer radius (drawing only)            [T: 3 mm]
P.d         = 2e-3;       % [m]   element axis to gel surface (also d in alpha)  [T]

% ------------------------------------------------------ MEDIA (solve) -----
P.eps_g1    = 43.98;      % [-]   gel permittivity, real                      [T]
P.eps_g2    = 13.5;       % [-]   gel permittivity, loss factor -> sigma = w eps0 eps'' = 1.84 S/m
                          %       (the table quotes 1.7 S/m)                  [T]
P.eps_w1    = 73;         % [-]   water, real                                 [T]
P.eps_w2    = 8.8;        % [-]   water, loss factor                          [T]
P.rho       = 1000;       % [kg/m^3] gel density (1 % Phytagel)               [V4]
P.c         = 4180;       % [J/(kg K)]                                        [V4]
P.k         = 0.60;       % [W/(m K)]                                         [V4]

% ------------------------------------------------- BOUNDARIES (solve) -----
P.h_s       = 1400;       % [W/(m^2 K)] water-to-gel surface; Inf pins it     [Doc, V4]
P.h_w       = 8;          % [W/(m^2 K)] tank walls to room; 0 = adiabatic     [Doc, V4]
P.L_c       = NaN;        % [m]   heated catheter length on the surface;
                          %       NaN = the element length L                  [Doc]

% --------------------------------------------- SURFACE TRANSFER (instant) --
P.beta      = 9e-4;       % [m (L/min)^0.5] alpha = min(1, beta/(d sqrt f)) = 0.2  [Doc guess]

% ---------------------------------------------- CATHETER (instant) --------
P.m_cat     = 2e-3;       % [kg]  heated tip mass                             [Doc guess]
P.c_cat     = 1500;       % [J/(kg K)]                                        [Doc guess]
P.k_cat     = 0.1;        % [W/K] tip-to-water loss                           [Doc guess]
P.gamma     = 0.043;      % [-]   share of the lost power heating the tip     [Doc guess]

% ------------------------------------------------- GEL BLOCK (solve) ------
% Solved window: the table's 30 x 30 mm gel enlarged 20 % (length, width) and
% 30 % deeper at the user's request (2026-09-28), so the heat front and the
% catheter-footprint edges stay clear of the block faces over a 10 min run.
P.gel_L     = 36e-3;      % [m]   length along the antenna                    [T x 1.2]
P.gel_H     = 39e-3;      % [m]   height (depth below the surface)            [T x 1.3]
P.gel_W     = 36e-3;      % [m]   width across the antenna (ASSUMED = length)
P.x_start   = -4.5e-3;    % [m]   gel edge behind the feed (feed at x = 0);
                          %       ASSUMED: 27 mm element centred on the block

% ------------------------------------------------ NUMERICS (solve) --------
P.dim       = 3;          % 2 = (x,z) slice, infinite in y (fast, overestimates); 3 = full block
P.dh        = 0.5e-3;     % [m]   grid spacing
P.t_end     = 600;        % [s]   length of the stored responses
P.dt_frame  = 2;          % [s]   frame interval (exact in time: modal solution)

% ------------------------------------------------ PIGMENT (display) -------
P.T_tlc     = 40.5;       % [degC] thermochromic start (first pink).  No datasheet:
                          %        FITTED to the photographed bloom edges, which sit at
                          %        40-42 degC in this model at 29, 91, 181 and 332 s
                          %        (photo preset, 3-D, 0.35 mm).  Trades off with P_dropoff.
P.bw_tlc    = 23;         % [K]    start -> blue.  FITTED: the core is pink at 29 s, green
                          %        at 91 s, teal with blue hints at 332 s (peak 44/53/62 degC)

switch lower(preset)
    case 'table'
    case 'photo'          % 50 W CW, 6 min, 2026-08-25
        P.P_set   = 50;
        P.t_on    = 360;
        P.T_b     = 25;   P.T_0 = 25;  P.T_room = 25;   % bath and gel at 25 degC on the day [v2 README]
        P.L       = 12e-3;                              % [Ph]
        P.a       = 0.15e-3;  P.b = 0.42e-3;            % RG-178 with the PTFE left on [Ph]
        P.d       = 1.1e-3;                             % [Ph]
        P.r_cable = 0.9e-3;                             % RG-178 braid + jacket [Ph]
        P.x_start = -9.5e-3;                            % photographed block -6.5..23.5 mm, grown 3 mm each side [Ph]
    otherwise
        error('v5_defaults:preset', 'unknown preset "%s" (table | photo)', preset);
end
end
