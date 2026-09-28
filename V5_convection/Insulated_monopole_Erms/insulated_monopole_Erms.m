%% INSULATED MONOPOLE - E_rms in the gel below the photographed antenna (50W CW 6 min)
% King insulated-antenna current on the element + thin-wire field (current and charge) in the gel.
clear; clc;

%% Constants
c0   = 299792458;        % speed of light (m/s)
mu0  = 4e-7*pi;          % permeability of free space (H/m)
eps0 = 8.854e-12;        % permittivity of free space (F/m)
f    = 2.45e9;           % microwave frequency (Hz)
w    = 2*pi*f;           % angular frequency (rad/s)
k0   = w/c0;             % free-space wavenumber (rad/m)

%% Power
P_set     = 50;          % generator setting in the photographed run (W)
P_dropoff = 0.7;         % fraction lost between generator and monopole
P = P_set*(1 - P_dropoff);                 % power delivered to the antenna (W)

%% Geometry (photographed antenna, measured from Gel_model_V3/mirrored.jpg at 105.8 px/mm)
L   = 12e-3;             % element length, feed to tip (m)
a   = 0.15e-3;           % inner conductor radius, RG-178 (m)
b   = 0.42e-3;           % PTFE insulation outer radius, RG-178 (m)
l_f = 1.5e-3;            % exposed metal feed block (solder/ferrule) behind the feed (m)
h   = 1.1e-3;            % element axis above the gel surface (m)

%% Media (e^{jwt} convention: eps = eps' - j eps'')
eps_ins = 2.1;                             % PTFE
eps_w   = 73 - 8.8j;                       % water around the element
eps_g   = 43.98 - 13.5j;                   % gel
k_ins = k0*sqrt(eps_ins);  k_w = k0*sqrt(eps_w);  k_g = k0*sqrt(eps_g);

%% 1. Wavenumber along the insulated element (King, Trembly & Strohbehn, IEEE MTT-31, 1983)
% k_L = k_ins * sqrt( 1 + H0(k_w b) / ( k_w b H1(k_w b) ln(b/a) ) )
k_L = k_ins*sqrt(1 + besselh(0,2,k_w*b)/(k_w*b*besselh(1,2,k_w*b)*log(b/a)));

%% 2. Feed current from the delivered power (open-ended insulated line)
% Z_in = -j Z_c cot(k_L L),  Z_c = k_L/(w C'),  C' = 2 pi eps0 eps_ins / ln(b/a)
C_ins = 2*pi*eps0*eps_ins/log(b/a);
Z_in  = -1j*(k_L/(w*C_ins))*cot(k_L*L);
I0    = sqrt(2*P/real(Z_in));              % peak feed current (A)

%% 3. Current and charge along the antenna
% element     (0 < x' < L):    I = I0 sin(k_L(L - x'))/sin(k_L L)   standing wave, zero at the tip
% feed block  (-l_f < x' < 0): I = I0 (1 + x'/l_f)                  return charge spread over the block
% line charge (continuity):    q = (j/w) dI/dx'
dx  = 0.02e-3;
xs  = [-l_f+dx/2 : dx : -dx/2,  dx/2 : dx : L-dx/2];
blk = xs < 0;
I = I0*sin(k_L*(L - xs))/sin(k_L*L);             I(blk) = I0*(1 + xs(blk)/l_f);
q = (1j/w)*(-I0*k_L*cos(k_L*(L - xs))/sin(k_L*L)); q(blk) = (1j/w)*I0/l_f;

%% 4. Field in the gel: E = -jwA - grad(phi), Green's function G = exp(-jk_g R)/R
% A = mu0/(4 pi) sum I G dx,  phi = tau/(4 pi eps0 eps_g) sum q G dx
% tau = 2 eps_g/(eps_w + eps_g): the charges sit in the water above the gel (quasi-static image)
tau = 2*eps_g/(eps_w + eps_g);
x = -5e-3  : 0.1e-3 : L + 5e-3;            % along the antenna, feed at x = 0 (m)
z = -10e-3 : 0.1e-3 : 0;                   % vertical, z = 0 at gel surface (m)
[X, Z] = meshgrid(x, z);
RHO = h - Z;                               % distance from the element axis (m)
Ex = 0;  Erho = 0;
for n = 1:numel(xs)
    DX = X - xs(n);  R = sqrt(DX.^2 + RHO.^2);
    G  = exp(-1j*k_g*R)./R;
    dG = -(1 + 1j*k_g*R).*G./R;            % dG/dR
    Ex   = Ex - 1j*w*mu0/(4*pi)*I(n)*G*dx - tau*q(n)*dx/(4*pi*eps0*eps_g)*dG.*DX./R;
    Erho = Erho - tau*q(n)*dx/(4*pi*eps0*eps_g)*dG.*RHO./R;
end
E_rms = sqrt((abs(Ex).^2 + abs(Erho).^2)/2);   % (V/m)

%% Result
[~, ix] = max(E_rms, [], 2);               % x of the peak at each depth
fprintf('k_L = %.1f %+.1fj rad/m (lambda_L/4 = %.1f mm),  Z_in = %.1f %+.1fj ohm,  I0 = %.2f A\n', ...
        real(k_L), imag(k_L), pi/2/real(k_L)*1e3, real(Z_in), imag(Z_in), I0);
for d = [0 1 2 3 5 10]
    iz = find(abs(z + d*1e-3) < 1e-9);
    fprintf('depth %4.1f mm: peak E_rms %5.0f V/m at x = %4.1f mm\n', d, E_rms(iz, ix(iz)), x(ix(iz))*1e3);
end

%% Plot (same style as Gel_model_V3/run_E_rms_2d.m)
fg = figure('Name', 'E_rms', 'Color', 'w', 'Position', [80 80 1000 560]);
ax = axes(fg);  ax.Toolbar.Visible = 'off';
imagesc(ax, x*1e3, z*1e3, log10(E_rms));
set(ax, 'YDir', 'normal');  axis(ax, 'image');  hold(ax, 'on');
colormap(ax, turbo);  clim(ax, log10([100 6000]));
levels = [50 70 100 150 200 300 500 700 1000 1500 2000 3000 5000];
[C, hc] = contour(ax, x*1e3, z*1e3, E_rms, levels, 'k-', 'LineWidth', 0.6);
clabel(C, hc, 'FontSize', 8, 'Color', 'k', 'LabelSpacing', 400);
plot(ax, x(ix)*1e3, z*1e3, 'w--', 'LineWidth', 1.2);                              % peak at each depth
plot(ax, 3.35, -3.35, 'wo', 'MarkerSize', 9, 'LineWidth', 1.5);                   % photographed hot spot, 92 s
plot(ax, [0 L]*1e3, [h h]*1e3, 'k-', 'LineWidth', 4);                             % element
plot(ax, [x(1) 0]*1e3, [h h]*1e3, '-', 'Color', [.5 .5 .5], 'LineWidth', 4);      % feed cable
ylim(ax, [-10 h*1e3 + 1]);  xlabel(ax, 'x (mm)');  ylabel(ax, 'z (mm)');
cb = colorbar(ax, 'Ticks', log10([100 200 300 500 1000 2000 3000 5000]), ...
              'TickLabels', {'100','200','300','500','1k','2k','3k','5k'});
cb.Label.String = 'E_{rms}  (V/m)';
title(ax, {sprintf('Time-RMS E-field in the gel plane, insulated %g mm monopole, P_{set} = %g W (%g W delivered), 2.45 GHz', ...
                   L*1e3, P_set, P), ...
           'dashed: x of the peak at each depth;  o: photographed hot spot at 92 s'});

exportgraphics(fg, 'Erms_insulated_monopole.png', 'Resolution', 300);
