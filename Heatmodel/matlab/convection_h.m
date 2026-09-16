function H = convection_h(varargin)
%CONVECTION_H  Convective heat-transfer coefficient at the gel surface for a
%              5 L/min flow of 37 degC water (the blood analogue).
%
%   H = CONVECTION_H()                        baseline
%   H = CONVECTION_H('Qflow',5/60000,'Wch',30e-3,'Hch',30e-3,'Lplate',30e-3)
%
%   Returns every intermediate quantity plus FOUR correlation estimates,
%   because the answer depends entirely on a flow geometry the brief does not
%   specify.  The spread between them is the honest uncertainty and is swept
%   in RUN_HEAT_SENSITIVITY.
%
%   WHAT IS SPECIFIED AND WHAT IS NOT
%     specified : coolant = water, T = 37 degC, volumetric flow = 5 L/min
%     NOT specified : the cross-section it flows through, hence the velocity;
%                     where the 30 mm gel patch sits relative to the leading
%                     edge of the flow.
%   Both are needed for h.  The baseline assumes the phantom sits in the floor
%   of a channel whose cross-section equals the region of interest,
%   30 mm x 30 mm, with the patch starting at the leading edge.
%
%   THE FOUR ESTIMATES
%   (a) 'plate_lam'  laminar boundary layer growing from the leading edge of a
%       flat plate of length L:      Nu_L = 0.664 Re_L^(1/2) Pr^(1/3)
%       [Incropera & DeWitt, Fundamentals of Heat and Mass Transfer, Ch. 7,
%        average Nusselt number for laminar flow over an isothermal plate;
%        valid Pr >~ 0.6, Re_L < 5e5].  This is the BASELINE, because both
%       entrance lengths below are metres long, so over a 30 mm patch the
%       flow is an undeveloped boundary layer, not duct flow.
%   (b) 'plate_turb' turbulent flat plate, mixed boundary layer:
%       Nu_L = (0.037 Re_L^(4/5) - 871) Pr^(1/3)      [same source]
%   (c) 'duct_dittus' Dittus-Boelter for a duct:
%       Nu_D = 0.023 Re_D^(4/5) Pr^0.4  (heating)     [same source]
%   (d) 'duct_lam'   fully developed laminar duct flow, one wall heated at
%       constant temperature, Nu_D ~ 4.9 for a square duct with three
%       adiabatic walls (order-of-magnitude bound, the LOW end).
%
%   ENTRANCE-LENGTH CHECK (reported): for laminar duct flow
%       L_hydro ~ 0.05 Re_D D_h ,    L_therm ~ 0.05 Re_D Pr D_h .
%   If both greatly exceed the 30 mm patch, the flat-plate estimates apply.
%
%   The catheter is treated separately: it is a small cylinder in cross flow,
%   Churchill-Bernstein,
%       Nu_D = 0.3 + 0.62 Re^(1/2) Pr^(1/3) / [1+(0.4/Pr)^(2/3)]^(1/4)
%              * [1 + (Re/282000)^(5/8)]^(4/5)
%   in series with the conduction resistance of its wall.  See H.cath.

p = struct( ...
    'Tw',      37, ...            % [degC] coolant temperature
    'fluid',   'water', ...       % 'water' (as specified) or 'blood'
    'Qflow',   5/60000, ...       % [m^3/s] 5 L/min
    'Wch',     30e-3, ...         % [m] channel width  (ASSUMED = ROI width)
    'Hch',     30e-3, ...         % [m] channel height (ASSUMED)
    'Lplate',  30e-3, ...         % [m] streamwise length of the heated patch
    'Dcath',   1.5e-3, ...        % [m] catheter outer diameter (ASSUMED)
    'twall',   0.25e-3, ...       % [m] catheter wall thickness (ASSUMED)
    'kwall',   0.25, ...          % [W/mK] PTFE/Pebax catheter wall
    'correlation','plate_lam');
for q = 1:2:numel(varargin), p.(varargin{q}) = varargin{q+1}; end

if strcmpi(p.fluid,'blood'), W = water_properties('blood');
else,                        W = water_properties(p.Tw);  end

A   = p.Wch*p.Hch;                       % [m^2] flow cross-section
U   = p.Qflow/A;                         % [m/s] bulk velocity
Dh  = 4*A/(2*(p.Wch+p.Hch));             % [m] hydraulic diameter
ReD = U*Dh/W.nu;
ReL = U*p.Lplate/W.nu;
Pr  = W.Pr;

% ---- (a) laminar flat plate, average over 0..L ---------------------------
Nu_a = 0.664*sqrt(ReL)*Pr^(1/3);
% ---- (b) turbulent flat plate, mixed boundary layer ----------------------
Nu_b = max(0.037*ReL^0.8 - 871, 0)*Pr^(1/3);
% ---- (c) Dittus-Boelter duct --------------------------------------------
Nu_c = 0.023*ReD^0.8*Pr^0.4;
% ---- (d) fully developed laminar duct, one wall heated -------------------
Nu_d = 4.9;

h_a = Nu_a*W.k/p.Lplate;
h_b = Nu_b*W.k/p.Lplate;
h_c = Nu_c*W.k/Dh;
h_d = Nu_d*W.k/Dh;

switch lower(p.correlation)
    case 'plate_lam',   h = h_a;
    case 'plate_turb',  h = h_b;
    case 'duct_dittus', h = h_c;
    case 'duct_lam',    h = h_d;
    otherwise, error('convection_h:corr','unknown correlation ''%s''',p.correlation);
end

% ---- the catheter --------------------------------------------------------
% Outer surface: cylinder in cross flow (Churchill-Bernstein).
ReC = U*p.Dcath/W.nu;
Nu_c1 = 0.3 + (0.62*sqrt(ReC)*Pr^(1/3))/(1+(0.4/Pr)^(2/3))^0.25 * ...
        (1 + (ReC/282000)^(5/8))^(4/5);
h_cath_out = Nu_c1*W.k/p.Dcath;
% Series with the wall conduction resistance (plane-wall approximation, the
% wall is thin compared with the radius).
R_wall = p.twall/p.kwall;                       % [m^2 K/W]
h_cath = 1/(1/h_cath_out + R_wall);

H = struct();
H.params = p;   H.fluid = W;
H.A = A;  H.U = U;  H.Dh = Dh;  H.ReD = ReD;  H.ReL = ReL;  H.Pr = Pr;
H.L_hydro = 0.05*ReD*Dh;                        % laminar entrance lengths
H.L_therm = 0.05*ReD*Pr*Dh;
H.Nu = struct('plate_lam',Nu_a,'plate_turb',Nu_b,'duct_dittus',Nu_c,'duct_lam',Nu_d);
H.h  = struct('plate_lam',h_a,'plate_turb',h_b,'duct_dittus',h_c,'duct_lam',h_d);
H.h_used = h;
H.correlation = p.correlation;
H.cath = struct('Re',ReC,'Nu',Nu_c1,'h_out',h_cath_out,'R_wall',R_wall, ...
                'h_eff',h_cath);
% "cooling length": the thickness of gel whose conduction resistance equals
% the surface film resistance.  Dirichlet (h -> Inf) is a good approximation
% only if this is small compared with the depth over which power is deposited.
H.note = ['cooling length k_gel/h; compare with the ~1.3 mm 1/e depth of ' ...
          'the volumetric power deposition'];

% ---- bulk temperature rise of the coolant (the user asked us to ignore it;
%      this quantifies that it is safe to do so) ---------------------------
H.mdot     = W.rho*p.Qflow;                     % [kg/s]
H.dT_per_W = 1/(H.mdot*W.cp);                   % [K/W] coolant rise per watt
end
