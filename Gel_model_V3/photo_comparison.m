%PHOTO_COMPARISON  Why the modelled heating pattern does not look like the
%                  50 W / 6 min photographs -- the tests behind
%                  gel_heating_notes.md s5-s6.  Prints tables, draws nothing.
%
%   Same gel, grid spacing, thermal properties and wall boundaries as
%   gel_heating_solve; 50 W set (10 W delivered); the rise above the no-power
%   steady state is reported at the times of the photographs (IMG_7908, 7910,
%   7918 = 29, 92, 332 s after IMG_7907).  The heat equation is linear, so the
%   rise does not depend on the bath temperature (37 degC here, 25 degC on the
%   photographed day) -- only the thermochromic threshold does.
%
%   The photographs (../50W CW 6 min, scale from the 30 mm block: 12 mm element,
%   105.8 px/mm in mirrored.jpg; feed x = 430 px, gel surface y = 355 px):
%     wire axis 1.1 mm above the gel surface
%     hottest (teal) spot     x ~ 3.1-3.6 mm along the element, 3.3-3.6 mm deep
%                             from 92 s on; 2.5-3 mm deep already at 29 s
%     cool top layer          1.1 / 0.8 / 0.5 / 0.4 / 0.3 mm at 29/62/92/182/332 s
%     second, diffuse lobe    x ~ 6-9 mm, 5-6 mm deep, visible from ~130 s
%     under the braid         faint lobe only
%   Metrics printed for each case:
%     x_pk, z_pk    position of the peak rise on the antenna plane
%     core          centroid of the region where the rise >= half the peak rise
%                   (the green/teal zone of a narrow-band pigment), and the
%                   depth of its top edge
%     s/p           rise at the surface node / peak rise, same column
%
%   ~6 min.  Run from this folder:  >> photo_comparison

rho = 1000;  c = 4180;  k = 0.60;  alpha = k/(rho*c);
h_wall = 7.9;  T_w = 37;  T_room = 25;  P_set = 50;  eta = 0.2;  a = 0.45e-3;
eps0 = 8.8541878128e-12;  mu0 = 1.25663706212e-6;  k0 = 2*pi*2.45e9*sqrt(mu0*eps0);
kc = k0*sqrt(43.98 - 1j*13.5);  kc = real(kc) - 1j*abs(imag(kc));   % gel, 344-52j
times = [29 92 332];

%      name                                              L      gap     E-field options                h_top   pinned  Ez-weight
cases = {
 'V3 as run: 27 mm, arm 27 mm, ke = k_gel, gap 0.5 mm',  27e-3, 0.5e-3, {},                              1460,   false,  1;
 '27 mm, no braid arm',                                  27e-3, 0.5e-3, {'L_sheath', 0},                 1460,   false,  1;
 '12 mm (photo), arm 12 mm, ke = k_gel, gap 0.65 mm',    12e-3, 0.65e-3, {},                             1460,   false,  1;
 '12 mm, arm 4 mm, ke = k_gel',                          12e-3, 0.65e-3, {'L_sheath', 4e-3},             1460,   false,  1;
 '12 mm, arm 4 mm, ke = k_gel/2 (insulated wire, water)', 12e-3, 0.65e-3, {'L_sheath', 4e-3, 'ke', kc/2}, 1460,  false,  1;
 '12 mm, arm 4 mm, ke = k_gel/2, surface pinned',        12e-3, 0.65e-3, {'L_sheath', 4e-3, 'ke', kc/2}, 1460,   true,   1;
 '12 mm, arm 4 mm, quarter-wave ke = pi/2L',             12e-3, 0.65e-3, {'L_sheath', 4e-3, 'ke', pi/24e-3}, 1460, false, 1;
 '12 mm, arm 4 mm, ke = k_gel, E_x ONLY (diagnostic)',   12e-3, 0.65e-3, {'L_sheath', 4e-3},             1460,   false,  0;
};

fprintf('%-58s %6s | %s\n', 'case (50 W set, 10 W delivered)', 'P_gel', 't [s]: peak rise K @ (x, z) mm | core z, top | s/p');
for ci = 1:size(cases, 1)
    [name, L_el, gap, opts, h_top, pinned, wz] = cases{ci, :};
    h_ax = a + gap;  dh = 0.5e-3;
    x = -6e-3:dh:(L_el + 6e-3);  y = -12e-3:dh:12e-3;  z = -24e-3:dh:0;
    nx = numel(x);  ny = numel(y);  nz = numel(z);
    [X3, Y3, Z3] = ndgrid(x, y, z);  iy0 = find(abs(y) < 1e-9, 1);  [Xp, Zp] = ndgrid(x, z);

    % ---- source from the E-field model (axisymmetric about the wire) ----
    rho_tab = logspace(log10(0.6e-3), log10(70e-3), 200);  [Xt, Rt] = ndgrid(x, rho_tab);
    [~, S] = E_field_2d(Xt, h_ax - Rt, 0, P_set, 'L', L_el, 'L_sheath', L_el, opts{:}, 'h', h_ax, 'eta', eta);
    qtab = 0.5*S.sig_eff*(abs(S.Ex).^2 + wz*abs(S.Ez).^2);
    q = interp2(log(rho_tab), x, qtab, log(max(hypot(Y3, Z3 - h_ax), rho_tab(1))), X3, 'linear');
    q(isnan(q)) = 0;  P_gel = sum(q(:))*dh^3;

    % ---- no-power steady state (as gel_heating_solve) ----
    g_top = 2*dh*h_top/k;  g_w = 2*dh*h_wall/k;
    [Lx, bx] = lap1d(nx, g_w, g_w, T_room, T_room);
    [Ly, by] = lap1d(ny, g_w, g_w, T_room, T_room);
    [Lz, bz] = lap1d(nz, g_w, g_top, T_room, T_w);
    Ix = speye(nx);  Iy = speye(ny);  Iz = speye(nz);
    L3 = kron(kron(Iz, Iy), Lx) + kron(kron(Iz, Ly), Ix) + kron(kron(Lz, Iy), Ix);
    B3 = bx(:) + reshape(by, 1, []) + reshape(bz, 1, 1, []);
    wx = ones(nx, 1);  wx([1 end]) = 0.5;  wy = ones(ny, 1);  wy([1 end]) = 0.5;
    wzz = ones(nz, 1);  wzz([1 end]) = 0.5;  W3 = kron(kron(wzz, wy), wx);
    A  = -spdiags(W3, 0, numel(W3), numel(W3))*L3;  Lc = ichol(A);
    [Tss, flag] = pcg(A, W3.*B3(:), 1e-9, 4000, Lc, Lc');  assert(flag == 0);
    Tss = reshape(Tss, nx, ny, nz);
    if pinned, Tss(:, :, end) = T_w; end

    % ---- FTCS with the source, report at the photo times ----
    T = Tss;  Tg = zeros(nx+2, ny+2, nz+2);  t = 0;
    dt0 = 0.9*dh^2/(alpha*(6 + g_top + 2*g_w));
    line = sprintf('%-58s %5.2f W |', name, P_gel);
    for tt = times
        n = ceil((tt - t)/dt0);  dt = (tt - t)/n;  r = alpha*dt/dh^2;  src = q/(rho*c)*dt;
        for s = 1:n
            Tg(2:end-1, 2:end-1, 2:end-1) = T;
            Tg(1,   2:end-1, 2:end-1) = T(2,     :, :) - g_w*(T(1,   :, :) - T_room);
            Tg(end, 2:end-1, 2:end-1) = T(end-1, :, :) - g_w*(T(end, :, :) - T_room);
            Tg(2:end-1, 1,   2:end-1) = T(:, 2,     :) - g_w*(T(:, 1,   :) - T_room);
            Tg(2:end-1, end, 2:end-1) = T(:, end-1, :) - g_w*(T(:, end, :) - T_room);
            Tg(2:end-1, 2:end-1, 1  ) = T(:, :, 2    ) - g_w  *(T(:, :, 1  ) - T_room);
            Tg(2:end-1, 2:end-1, end) = T(:, :, end-1) - g_top*(T(:, :, end) - T_w);
            T = T + r*( Tg(1:end-2,2:end-1,2:end-1) + Tg(3:end,2:end-1,2:end-1) ...
                      + Tg(2:end-1,1:end-2,2:end-1) + Tg(2:end-1,3:end,2:end-1) ...
                      + Tg(2:end-1,2:end-1,1:end-2) + Tg(2:end-1,2:end-1,3:end) - 6*T ) + src;
            if pinned, T(:, :, end) = T_w; end
        end
        t = tt;
        rise = squeeze(T(:, iy0, :)) - squeeze(Tss(:, iy0, :));
        [rpk, i] = max(rise(:));  [ix, iz] = ind2sub(size(rise), i);
        core = rise >= 0.5*rpk;
        line = sprintf('%s %3.0f: %4.1f @ (%4.1f,%5.1f) | %4.1f, %3.1f | %.2f |', line, tt, rpk, ...
                       x(ix)*1e3, z(iz)*1e3, -mean(Zp(core))*1e3, -max(Zp(core))*1e3, rise(ix, end)/rpk);
    end
    fprintf('%s\n', line);
end
fprintf(['\nphoto: hottest spot x 3.1-3.6, depth ~2.5-3 (29 s), 3.3-3.6 mm (92-332 s); top of coloured zone ' ...
         '1.1 -> 0.3 mm; s/p: surface below the pigment threshold\n']);

function [L, b] = lap1d(n, g_lo, g_hi, T_lo, T_hi)
e = ones(n, 1);  L = spdiags([e -2*e e], -1:1, n, n);
L(1, 2) = 2;  L(1, 1) = -(2 + g_lo);  L(n, n-1) = 2;  L(n, n) = -(2 + g_hi);
b = zeros(n, 1);  b(1) = g_lo*T_lo;  b(n) = g_hi*T_hi;
end
