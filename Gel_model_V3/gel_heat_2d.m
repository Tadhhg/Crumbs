%GEL_HEAT_2D  Transient 2-D heat conduction in the 30 mm x 30 mm gel of the
%             bench setup, with convective boundaries.
%
%   Solves the heat equation
%
%       dT/dt = alpha * (d2T/dx2 + d2T/dy2),      alpha = k/(rho c)
%
%   on 0 <= x,y <= 30 mm.  Every face is a CONVECTIVE (Robin) boundary,
%
%       -k dT/dn = h (T - T_inf)
%
%   because nothing on the bench pins a temperature:
%       top     (y = 30 mm)  37 degC water circulating at 5 L/min over the
%                            gel.  h from the laminar flat-plate correlation,
%                            ~1500 W/(m^2 K), so Bi = hL/k ~ 70 and the face
%                            sits within ~0.1 K of 37 degC -- effectively fixed.
%       sides & bottom       the acrylic tank wall, whose outside is in still
%                            room air at 25 degC.  Free convection + radiation
%                            in series with the wall conduction gives
%                            ~8 W/(m^2 K), Bi ~ 0.4.  This face is NOT fixed at
%                            25 degC: it floats a few degrees below the gel
%                            behind it.  h_air = 0 is the adiabatic limit and
%                            a very large h_air recovers the fixed-25 degC case,
%                            so both brackets are one edit away.
%       initial condition    T = 25 degC everywhere (gel at room temperature)
%
%   Explicit forward-time / centred-space (FTCS) finite differences with a
%   ghost-node layer for the Robin faces.  Everything happens in this one
%   script and it produces ONE figure: the temperature map is redrawn as the
%   simulation advances and ends on the final (steady-state) field.
%
%   Experimental constants are those in Experimental_setup_table.xlsx (same
%   folder).  Gel thermal properties are the ones in
%   ../Heatmodel/matlab/thermal_config.m (1 % Phytagel in saline, >98 % water).

clear; close all;

% ---------------------------------------------- EXPERIMENTAL SETUP TABLE --
% Experimental_setup_table.xlsx, verbatim.  This conduction-only model uses
% the water temperature, the circulation flow rate and the gel size; the
% permittivities, conductivity, power and antenna entries belong to the
% microwave source term, which is not in this script yet.
setup = struct();
setup.epsr_gel    = 43.98;      % [-]     relative real permittivity, gel
setup.epsi_gel    = 13.5;       % [-]     relative imaginary permittivity, gel
setup.sigma_gel   = 1.7;        % [S/m]   conductivity, gel
setup.T_water     = 37;         % [degC]  water temperature
setup.Q_irrig     = 0.03/60000; % [m^3/s] irrigation flow, 0.03 L/min
setup.Q_circ      = 5/60000;    % [m^3/s] circulation flow, 5 L/min
setup.epsr_water  = 73;         % [-]     relative real permittivity, water
setup.epsi_water  = 8.8;        % [-]     relative imaginary permittivity, water
setup.P_reflected = 0;          % [W]     reflected microwave power
setup.P_dropoff   = 0.7;        % [-]     power drop-off
setup.H_gel       = 30e-3;      % [m]     gel height
setup.L_gel       = 30e-3;      % [m]     gel length
setup.f           = 2.45e9;     % [Hz]    microwave frequency
setup.L_antenna   = 27e-3;      % [m]     monopole antenna length

% ------------------------------------------------------- GEL PROPERTIES --
rho   = 1000;                   % [kg/m^3]
c     = 4180;                   % [J/(kg K)]
k     = 0.60;                   % [W/(m K)]
alpha = k/(rho*c);              % [m^2/s] = 1.435e-7
T_init = 25;                    % [degC]  gel starts at room temperature

% ------------------------------------------ TOP FACE: circulating water --
% Laminar boundary layer over a 30 mm plate, Nu_L = 0.664 Re_L^(1/2) Pr^(1/3)
% (Incropera & DeWitt Ch. 7), as in ../Heatmodel/matlab/convection_h.m.  The
% flow cross-section is ASSUMED equal to the gel face, 30 x 30 mm.
nu_w  = 6.98e-7;   k_w = 0.628;   Pr_w = 4.61;    % water at 37 degC, Table A.6
U_w   = setup.Q_circ/(setup.L_gel*setup.H_gel);  % [m/s] bulk velocity, 0.093
Re_L  = U_w*setup.L_gel/nu_w;                    % 3980, laminar
h_top = 0.664*sqrt(Re_L)*Pr_w^(1/3)*k_w/setup.L_gel;   % [W/(m^2 K)] ~1460
T_top = setup.T_water;

% --------------------------------- SIDES AND BOTTOM: tank wall, room air --
T_room = 25;                    % [degC]  room air
h_air  = 10;                    % [W/(m^2 K)] outside of the wall to the room:
                                %   free convection on a vertical wall a few
                                %   cm tall at dT ~ 10 K is 5-7 (Incropera
                                %   Ch. 9, h ~ 1.42 (dT/L)^0.25) plus ~5 of
                                %   radiation for emissivity ~0.9.  0 = adiabatic.
t_wall = 5e-3;                  % [m]     acrylic wall thickness (ASSUMED
                                %         from the photographs)
k_wall = 0.19;                  % [W/(m K)] acrylic
h_wall = 1/(1/h_air + t_wall/k_wall);            % [W/(m^2 K)] ~7.9, in series
Bi_wall = h_wall*setup.H_gel/k;                  % ~0.4: far from fixed-T

% --------------------------------------------------------------- NUMERICS --
dh      = 0.25e-3;              % [m]     grid spacing (121 x 121 nodes)
t_end   = 4*3600;               % [s]     4 h.  With the top pinned and the
                                %         other faces nearly adiabatic the
                                %         slowest mode has tau ~ (2H/pi)^2/alpha
                                %         = 2540 s; 4 h is 5.7 tau.
t_frame = 60;                   % [s]     simulated time between figure refreshes
iso     = 27:2:35;              % [degC]  isotherms drawn on the map

% ------------------------------------------------------------------ GRID --
x = 0:dh:setup.L_gel;   y = 0:dh:setup.H_gel;
nx = numel(x);          ny = numel(y);
[X, Y] = meshgrid(x, y);                    % T(iy,ix): row 1 is the bottom face

% Robin ghost node: T_ghost = T_in - 2 dh (h/k) (T_face - T_inf).  The
% coefficient 2 dh h/k also enters the FTCS stability limit,
%     dt <= dh^2 / ( alpha (4 + 2 dh (h1 + h2)/k) )
% with h1 + h2 the two faces meeting at the worst corner.
g_top  = 2*dh*h_top/k;
g_wall = 2*dh*h_wall/k;
dt     = 0.9*dh^2/(alpha*(4 + g_top + g_wall));
nsteps = ceil(t_end/dt);
dt     = t_end/nsteps;                      % integer number of steps to t_end
r      = alpha*dt/dh^2;

% ----------------------------------------------------- INITIAL CONDITION --
T  = T_init*ones(ny, nx);
Tg = zeros(ny+2, nx+2);                     % T with a one-node ghost border

% ---------------------------------------------------------------- FIGURE --
figure('Name','2-D gel heating','Color','w','Position',[100 100 700 580]);
im = imagesc(x*1e3, y*1e3, T);
set(gca,'YDir','normal');  axis image;  hold on;
colormap(turbo);  clim([T_room T_top]);
cb = colorbar;  cb.Label.String = 'T (\circC)';
warning('off', 'MATLAB:contour:ConstantData');   % the t = 0 field is uniform
[~, hc] = contour(X*1e3, Y*1e3, T, iso, 'w-', 'LineWidth', 0.7, ...
                  'ShowText', 'on', 'LabelSpacing', 300);
xlabel('x (mm)');  ylabel('y (mm)');
hdr = sprintf('h_{top} = %.0f,  h_{wall} = %.1f W m^{-2} K^{-1}', h_top, h_wall);
ttl = title({'30 x 30 mm gel:  37 \circC water above,  tank wall to 25 \circC room air', ...
             [hdr '     t = 0 s']});

% ------------------------------------------------------------- TIME LOOP --
steps_per_frame = max(1, round(t_frame/dt));
for s = 1:nsteps
    % ghost layer from the Robin condition on each face
    Tg(2:end-1, 2:end-1) = T;
    Tg(end, 2:end-1) = T(end-1, :) - g_top *(T(end, :) - T_top);    % top
    Tg(1,   2:end-1) = T(2,     :) - g_wall*(T(1,   :) - T_room);   % bottom
    Tg(2:end-1, 1)   = T(:,     2) - g_wall*(T(:,   1) - T_room);   % left
    Tg(2:end-1, end) = T(:, end-1) - g_wall*(T(:, end) - T_room);   % right

    % 5-point Laplacian on every real node, ghosts supplying the faces
    T = T + r*( Tg(1:end-2, 2:end-1) + Tg(3:end, 2:end-1) ...
              + Tg(2:end-1, 1:end-2) + Tg(2:end-1, 3:end) - 4*T );

    if mod(s, steps_per_frame) == 0 || s == nsteps
        im.CData = T;
        hc.ZData = T;
        if s == nsteps
            ttl.String{2} = sprintf('%s     t = %.0f s  (steady state)', hdr, s*dt);
        else
            ttl.String{2} = sprintf('%s     t = %.0f s', hdr, s*dt);
        end
        drawnow limitrate;
    end
end
