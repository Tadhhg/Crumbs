%RUN_E_FIELD_2D  ONE figure: the E-field vector on the gel plane at eight
%                phases of the 2.45 GHz cycle (omega*t = 0, pi/4, ..., 7pi/4).
%
%   Colour  = instantaneous |E(x,z,t)| (log scale), arrows = direction of E.
%   Uses E_field_2d.m (see E_field_2d_notes.md).  Run from this folder:
%       >> run_E_field_2d
%
%   The 27 mm element is centred over the 30 mm gel: feed at x = 0, tip at
%   x = 27 mm, gel from x = -1.5 to 28.5 mm, z = 0 is the gel surface.

P_set  = 100;                               % [W] generator setting -> 20 W delivered
phases = (0:7)*pi/4;                        % [rad] omega*t for the eight panels
dx     = 0.25e-3;                           % [m] map resolution

% ------------------------------------------------------------------ MAP ---
x = (-1.5e-3 : dx : 28.5e-3);   z = (-30e-3 : dx : 0);
[X, Z] = meshgrid(x, z);
[Et, S] = E_field_2d(X, Z, phases, P_set);  % Et.x, Et.z : [nz x nx x 8]
Emag = hypot(Et.x, Et.z);                   % instantaneous magnitude

% ---------------------------------------------------------------- FIGURE --
f = figure('Name','E-field phases','Color','w','Position',[60 60 1500 760]);
tl = tiledlayout(f, 2, 4, 'TileSpacing','compact', 'Padding','compact');
title(tl, sprintf(['E-field of the 27 mm monopole in the gel plane, P_{set} = %g W ' ...
                   '(%g W delivered), f = 2.45 GHz'], P_set, S.P_del), 'FontWeight','bold');

lim   = [2 4.2];                            % log10 V/m colour range
qstep = round(2.5e-3/dx);                   % one arrow every 2.5 mm
iq = 1:qstep:numel(z);  jq = 1:qstep:numel(x);
for k = 1:numel(phases)
    ax = nexttile(tl);
    imagesc(ax, x*1e3, z*1e3, log10(Emag(:,:,k)));
    set(ax, 'YDir','normal');  axis(ax, 'image');  hold(ax, 'on');
    ax.Toolbar.Visible = 'off';
    clim(ax, lim);  colormap(ax, turbo);
    % unit-length arrows showing the direction of E at this phase
    ex = Et.x(iq,jq,k);  ez = Et.z(iq,jq,k);  m = hypot(ex, ez) + eps;
    quiver(ax, x(jq)*1e3, z(iq)*1e3, ex./m, ez./m, 0.45, 'w', 'LineWidth', 0.8);
    plot(ax, [0 27], [0.5 0.5], 'k-', 'LineWidth', 3);        % the element
    plot(ax, [x(1) 0]*1e3, [0.5 0.5], '-', 'Color', [.5 .5 .5], 'LineWidth', 3);  % braid
    ylim(ax, [-30 2]);
    title(ax, sprintf('\\omegat = %s', frac_pi(phases(k))));
    if k > 4, xlabel(ax, 'x (mm)'); end
    if mod(k,4) == 1, ylabel(ax, 'z (mm)'); end
end
cb = colorbar(ax, 'Ticks', log10([100 300 1000 3000 10000]), ...
              'TickLabels', {'100','300','1k','3k','10k'});
cb.Layout.Tile = 'east';  cb.Label.String = '|E(t)|  (V/m)';

exportgraphics(f, 'E_field_phases.png', 'Resolution', 150);

function s = frac_pi(p)
% '0', '\pi/4', '\pi/2', '3\pi/4', '\pi', ...
n = round(p/(pi/4));
switch n
    case 0, s = '0';
    case 1, s = '\pi/4';
    case 2, s = '\pi/2';
    case 4, s = '\pi';
    otherwise
        if mod(n,2) == 0, s = sprintf('%d\\pi/2', n/2);
        else,             s = sprintf('%d\\pi/4', n);
        end
end
end
