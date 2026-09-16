function T = thermal_config(varargin)
%THERMAL_CONFIG  Central definition of every thermal constant, boundary
%                condition and grid parameter used by the heating model.
%
%   Companion to ../../matlab/model_config.m, which defines the ELECTRO-
%   MAGNETIC side.  Every entry carries the source it was taken from.
%   Override any field with name/value pairs, e.g.
%       thermal_config('surf_bc','robin','k_gel',0.55)
%
%   COORDINATES  (identical to the field model)
%     x  along the antenna axis; the 12 mm monopole spans 0 <= x <= 12 mm
%     y  transverse horizontal; the antenna lies on y = 0
%     z  height; gel surface z = 0, gel floor z = -30 mm, so z <= 0 in the gel
%
%   THE REGION OF INTEREST
%   The brief restricts attention to a 30 mm x 30 mm region below the antenna.
%   That square is taken to be the vertical plane y = 0, spanning 30 mm along
%   x and the full 30 mm depth of the gel -- i.e. exactly the plane on which
%   the field model reports |E|.  Heat conduction is nevertheless solved in
%   THREE dimensions on the 30 x 30 x 30 mm block containing it, because the
%   thermal diffusion length sqrt(4*alpha*t) reaches 2.4 mm at 10 s and 5.9 mm
%   at 60 s, which is comparable to the transverse (y) width of the deposition
%   lobe.  A 2-D solve in the y = 0 plane would ignore spreading in y and
%   over-predict the temperature substantially.

T = struct();

% ------------------------------------------------------------- COOLANT ----
% "constant heat sink (blood): water at 37 degC, 5 L/min" (user).
T.T_blood  = 37;            % [degC] coolant / far-field temperature
T.T_init   = 37;            % [degC] initial gel temperature (in equilibrium
                            %        with the coolant before the power is on)
T.Qflow    = 5/60000;       % [m^3/s] = 5 L/min (user)
T.coolant  = 'water';       % as specified by the user

% ------------------------------------------------- GEL THERMAL PROPERTIES --
% 1 %% w/v Phytagel in 0.9 %% w/v NaCl saline: >98 %% water by mass, so the
% thermal properties are those of water, slightly reduced by the polymer.
%
%   rho  1000 kg/m^3  -- kept EXACTLY equal to the density already used to
%        convert power density to SAR in ../../matlab/run_model.m, so that the
%        two models are numerically consistent.  Water at 37 degC is 993.
%   c    4180 J/(kg K) -- water at 37 degC is 4178 (Incropera Table A.6).
%        Measured agar tissue phantoms give 4090-4130 J/(kg K)
%        [Choi et al., J. Ther. Ultrasound 5:14 (2017),
%         doi:10.1186/s40349-017-0093-z].
%   k    0.60 W/(m K) -- water at 37 degC is 0.628.  Adding 1-2 %% of a
%        polysaccharide lowers it by a few percent: measured water-agar gels
%        [Zhang et al., J. Chem. Eng. Data 56:859-864 (2011),
%         doi:10.1021/je100570h] and measured agar tissue phantoms
%        (0.52 +- 0.06 and 0.57 +- 0.10 W/(m K), same J. Ther. Ultrasound
%        paper) bracket it.  Swept 0.50-0.63 in run_heat_sensitivity.
T.rho_gel = 1000;           % [kg/m^3]
T.c_gel   = 4180;           % [J/(kg K)]
T.k_gel   = 0.60;           % [W/(m K)]

% Temperature at which the gel PERMITTIVITY is evaluated, passed through to
% model_config.  Kept at 25 degC -- the baseline of the field model -- so the
% two models are numerically identical and the cross-validation in
% VERIFY_HEAT_MODEL test 5 is exact.  The gel actually sits in 37 degC blood;
% the Stogryn model gives sigma_eff = 2.723 S/m at 25 degC and 2.782 at 37,
% and eps' 74.29 vs 70.87, so the difference is a couple of percent.  Case H
% of run_heat_sensitivity quantifies it.
T.T_medium = 25;            % [degC]

% ---------------------------------------------------- POWER NORMALISATION --
% The field model is normalised to DELIVERED power P_del = P_fwd - P_ref.
% The variable the user asked for is the GENERATOR power, so carry the
% efficiency explicitly.  eta = 150/200 = 0.75 is the worst-case baseline of
% the field model (200 W forward, 50 W maximum reflected).
T.Pfwd_ref = 200;                       % [W] reference generator power
T.eta      = 150/200;                   % [-] delivered / forward
T.Pdel_ref = T.Pfwd_ref*T.eta;          % [W] = 150 W, matches the field model

% ----------------------------------------------------- BOUNDARY CONDITIONS --
% TOP (z = 0), the gel/blood interface.
%   'dirichlet' : T = T_blood everywhere on the surface.  THE BASELINE --
%                 this is what was asked for ("model the surface as just a
%                 constant 37C").  It is the perfect-heat-sink limit and
%                 therefore a LOWER BOUND on the gel temperature.
%   'robin'     : -k dT/dz = h (T - T_blood), h from convection_h().  The
%                 physical case; run and reported alongside the baseline.
T.surf_bc = 'dirichlet';
T.h_surf  = [];             % [W/(m^2 K)] filled from convection_h() if empty

% THE CATHETER.  It lies on the surface along 0 <= x <= L at y = 0 and is
% itself washed by the same 37 degC flow, so it is NOT a heat source and NOT
% an insulating blanket -- it is a second cooled surface.  Three treatments:
%   'transparent' : the surface BC is unchanged under the catheter footprint.
%                   Correct in the Dirichlet limit (a cooled body resting on a
%                   surface that is already pinned at 37 degC changes nothing)
%                   and the BASELINE.
%   'cooled'      : the footprint gets its own coefficient h_cath, computed
%                   from cylinder-in-cross-flow convection in series with the
%                   catheter wall conduction (convection_h().cath.h_eff).
%   'insulating'  : the footprint is adiabatic.  The OPPOSITE extreme -- a
%                   catheter that blocks the coolant and conducts nothing.
%                   Included to bracket the answer, not because it is likely.
T.cath_bc    = 'transparent';
T.cath_L     = 12e-3;       % [m] antenna length (= model_config L)
T.cath_x0    = -20e-3;      % [m] shaft extends back past the feed
T.cath_D     = 1.5e-3;      % [m] outer diameter (ASSUMED)
T.h_cath     = [];          % [W/(m^2 K)] filled from convection_h() if empty

% BOTTOM (z = -d).  The baseline of the field model is an air-backed dish, so
% the thermal baseline is adiabatic.  'dirichlet' (water-cooled dish) is swept.
T.bot_bc = 'adiabatic';

% SIDES (the four vertical faces of the ROI).  The instruction given -- "the
% rest can be assumed to be minor temperature rises taken care of by
% homeostasis" -- is exactly a Dirichlet condition at T_blood.  'adiabatic'
% is the opposite extreme (no heat leaves sideways) and is swept.
T.side_bc = 'dirichlet';

% ------------------------------------------------------ REGION OF INTEREST -
T.L_ant = 12e-3;                        % [m] monopole length
T.xc    = T.L_ant/2;                    % ROI centred on the antenna midpoint
T.Wroi  = 30e-3;                        % [m] ROI width in x and in y
T.d     = 30e-3;                        % [m] gel depth (= ROI depth)
T.xlim  = T.xc + [-1 1]*T.Wroi/2;       % [-9, +21] mm
T.ylim  =        [-1 1]*T.Wroi/2;       % [-15, +15] mm

% ------------------------------------------------ THE FEED SINGULARITY ----
% The field model prescribes the antenna current as I(x) = sin(k_e (L-x)) on
% 0 <= x <= L.  That current is MAXIMUM at the feed x = 0 and is then
% truncated, which puts a fictitious point charge there and makes |E| -- and
% so q -- blow up at (x,y,z) = (0,0,0).  It is an artefact of the prescribed
% current, not physics: at the tip x = L the current goes smoothly to zero and
% there is no such spike.  In a real applicator the current continues onto the
% coax outer conductor, which the field model does not include (its Sec. 11,
% limitation 3).  The deposition inside a couple of millimetres of the feed is
% therefore not trustworthy, and the peak temperature is reported EXCLUDING a
% sphere of this radius about the feed.  The unrestricted peak is reported too,
% and the radius is swept in run_heat_sensitivity.
%
% SUPERSEDED BY THE BRAID ARM.  model_config's L_sheath carries the return
% current back along the outside of the coax, so I(x) is continuous at the feed
% and there is no point charge left to hide.  GUI_MODEL therefore passes
% r_feed_excl = 0 whenever L_sheath > 0, and this radius applies only to the
% v1/v2 comparison case L_sheath = 0.  Setting it to 0 disables the mask.
T.r_feed_excl = 2e-3;       % [m]  0 = no exclusion (the braid arm handles it)
T.feed_xyz    = [0 0 0];    % the feed point
T.feed_fix    = 'none';     % 'none' (baseline) | 'cap'  -- see feed_regularise

% -------------------------------------------------------------- GRID ------
T.dh    = 0.25e-3;          % [m] uniform node spacing in x, y and z.
                            %     Chosen equal to dx of the field model so the
                            %     EM solution is evaluated directly on the
                            %     thermal nodes with no interpolation.
T.decim = 2;                % storage decimation for the 3-D history (0.5 mm)

% -------------------------------------------------------------- TIME ------
T.t_end = 600;              % [s] end of the fine transient run
T.Nt    = 241;              % number of output times
T.t_min = 1e-3;             % [s] first output time (log-spaced from here)
T.cfl   = 0.90;             % safety factor on the explicit stability limit

% ----------------------------------------------------- DERIVED / SOURCES ---
T.alpha_gel = T.k_gel/(T.rho_gel*T.c_gel);      % [m^2/s] thermal diffusivity

% ------------------------------------------------------------ OVERRIDES ---
for q = 1:2:numel(varargin), T.(varargin{q}) = varargin{q+1}; end
T.alpha_gel = T.k_gel/(T.rho_gel*T.c_gel);
T.Pdel_ref  = T.Pfwd_ref*T.eta;
T.xlim = T.xc + [-1 1]*T.Wroi/2;
T.ylim =        [-1 1]*T.Wroi/2;

% fill the convection coefficients if they were not supplied
if isempty(T.h_surf) || isempty(T.h_cath)
    H = convection_h('Tw',T.T_blood,'Qflow',T.Qflow,'fluid',T.coolant, ...
                     'Wch',T.Wroi,'Hch',T.Wroi,'Lplate',T.Wroi, ...
                     'Dcath',T.cath_D);
    if isempty(T.h_surf), T.h_surf = H.h_used;        end
    if isempty(T.h_cath), T.h_cath = H.cath.h_eff;    end
    T.conv = H;
end

T.t_out = unique([0, logspace(log10(T.t_min), log10(T.t_end), T.Nt-1)]);
end
