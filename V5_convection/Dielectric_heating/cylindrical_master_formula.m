%% E_rms calculation for monopole catheter antenna - cylindrical master formula
% E_rms in the gel below a 2.45 GHz monopole. Values from V5_convection/Experimental_setup_table.xlsx
clear; clc;

%% Constants
c    = 3e8;              % speed of light (m/s)
eta0 = 377;              % impedance of free space (Ohm)
f    = 2.45e9;           % microwave frequency (Hz)

%% Power
P_in      = 60;          % generator input power (W)
P_refl    = 0;           % reflected power (W)
P_dropoff = 0.7;         % fraction lost between generator and monopole
P = (P_in - P_refl) * (1 - P_dropoff);    % power radiated into the medium (W)

%% Geometry
L  = 27e-3;              % monopole length, feedpoint to tip (m)
r0 = 1.8e-3 / 2;         % monopole radius = coax inner diameter / 2 (m)
h  = 2e-3;               % monopole axis to gel surface (m)

%% Media
eps1_gel = 43.98;  eps2_gel = 13.5;       % gel relative permittivity (real, imaginary)
eps1_w   = 73;     eps2_w   = 8.8;        % water relative permittivity (real, imaginary)

%% Attenuation constant
% alpha = (2*pi*f/c) * sqrt( eps'/2 * ( sqrt(1 + (eps''/eps')^2) - 1 ) )
alpha   = @(e1, e2) (2*pi*f/c) * sqrt(e1/2 * (sqrt(1 + (e2/e1)^2) - 1));
alpha_g = alpha(eps1_gel, eps2_gel);      % gel (Np/m)
alpha_w = alpha(eps1_w,   eps2_w);        % water gap (Np/m)

%% Rectangle in the gel: 10 mm deep, 5 mm past the feedpoint (x = 0) and tip (x = L)
x = -5e-3  : 0.1e-3 : L + 5e-3;           % along the antenna (m)
z = -10e-3 : 0.1e-3 : 0;                  % vertical, z = 0 at gel surface (m)
[X, Z] = meshgrid(x, z);
D = -Z;                                   % depth below gel surface (m)

%% Distance to the nearest point on the antenna, split into water and gel path
x_past = max(0, max(-X, X - L));          % distance past feedpoint/tip, 0 alongside
r   = sqrt((h + D).^2 + x_past.^2);       % radial distance (m)
d_w = r .* h ./ (h + D) - r0;             % path through the water gap (m)
d_g = r .* D ./ (h + D);                  % path through the gel (m)

%% Cylindrical master formula
% E_rms = sqrt( P/(2*pi*r*L) * eta0/sqrt(eps') ) * exp(-alpha*(r - r0))
E_rms = sqrt(P ./ (2*pi*r*L) * eta0/sqrt(eps1_gel)) .* exp(-(alpha_w*d_w + alpha_g*d_g));
% Screenshot form (decay inside the sqrt):
% E_rms = sqrt(P ./ (2*pi*r*L) * eta0/sqrt(eps1_gel) .* exp(-(alpha_w*d_w + alpha_g*d_g)));

%% Result
fprintf('P = %.1f W,  alpha_gel = %.2f Np/m,  alpha_water = %.2f Np/m\n', P, alpha_g, alpha_w);
fprintf('E_rms: max %.0f V/m (surface under the element), min %.0f V/m\n', max(E_rms(:)), min(E_rms(:)));

%% Plot (same style as Gel_model_V3/run_E_rms_2d.m)
fg = figure('Name', 'E_rms', 'Color', 'w', 'Position', [80 80 1200 520]);
ax = axes(fg);  ax.Toolbar.Visible = 'off';
imagesc(ax, x*1e3, z*1e3, log10(E_rms));
set(ax, 'YDir', 'normal');  axis(ax, 'image');  hold(ax, 'on');
colormap(ax, turbo);  clim(ax, log10([300 2000]));
levels = [50 70 100 150 200 300 500 700 1000 1500 2000 3000 5000];
[C, hc] = contour(ax, x*1e3, z*1e3, E_rms, levels, 'k-', 'LineWidth', 0.6);
clabel(C, hc, 'FontSize', 8, 'Color', 'k', 'LabelSpacing', 400);
plot(ax, [0 L]*1e3, [h h]*1e3, 'k-', 'LineWidth', 4);                            % element
plot(ax, [x(1) 0]*1e3, [h h]*1e3, '-', 'Color', [.5 .5 .5], 'LineWidth', 4);     % feed cable
ylim(ax, [-10 h*1e3 + 1]);  xlabel(ax, 'x (mm)');  ylabel(ax, 'z (mm)');
cb = colorbar(ax, 'Ticks', log10([300 500 700 1000 1500 2000]), ...
              'TickLabels', {'300','500','700','1k','1.5k','2k'});
cb.Label.String = 'E_{rms}  (V/m)';
title(ax, {sprintf('Time-RMS E-field in the gel plane, P_{in} = %g W (%g W delivered), 2.45 GHz', P_in, P), ...
           sprintf('cylindrical master formula, \\alpha_{gel} = %.1f Np/m;  peak %.0f V/m at the surface under the element', ...
                   alpha_g, max(E_rms(:)))});

exportgraphics(fg, 'Erms_map.png', 'Resolution', 300);
