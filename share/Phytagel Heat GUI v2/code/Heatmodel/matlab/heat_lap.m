function L = heat_lap(th, G)
%HEAT_LAP  Second-order Laplacian of theta on the grid G, with the ghost-node
%          boundary conditions described in HEAT_GRID.  Pinned (Dirichlet)
%          nodes return 0.
%
%   theta is the temperature RISE above the coolant temperature, so every
%   boundary condition is homogeneous:  Dirichlet theta = 0, adiabatic
%   dtheta/dn = 0, Robin k dtheta/ds = h theta at the surface.  That is what
%   makes the whole problem linear in the delivered power (see the write-up).

nx = G.nx; ny = G.ny; nz = G.nz; dh2 = G.dh^2;

% ---------------------------------------------------------------- x -------
xp = zeros(nx,ny,nz);  xm = zeros(nx,ny,nz);
xp(1:nx-1,:,:) = th(2:nx,:,:);
xm(2:nx,:,:)   = th(1:nx-1,:,:);
if strcmp(G.side,'adiabatic')
    xp(nx,:,:) = th(nx-1,:,:);
    xm(1,:,:)  = th(2,:,:);
end                                  % 'dirichlet' -> ghost = 0, already set
L = (xp + xm - 2*th)/dh2;
clear xp xm

% ---------------------------------------------------------------- y -------
yp = zeros(nx,ny,nz);  ym = zeros(nx,ny,nz);
yp(:,1:ny-1,:) = th(:,2:ny,:);
ym(:,2:ny,:)   = th(:,1:ny-1,:);
if strcmp(G.side,'adiabatic')
    yp(:,ny,:) = th(:,ny-1,:);
    ym(:,1,:)  = th(:,2,:);
end
L = L + (yp + ym - 2*th)/dh2;
clear yp ym

% ---------------------------------------------------------------- z -------
% index 1 = gel surface (z = 0), index nz = gel floor (z = -d)
zp = zeros(nx,ny,nz);  zm = zeros(nx,ny,nz);
zp(:,:,1:nz-1) = th(:,:,2:nz);
zm(:,:,2:nz)   = th(:,:,1:nz-1);

% surface ghost (above z = 0).  Robin: theta_g = theta_2 - 2 dh h theta_1/k.
% Adiabatic (h = 0) reduces to the mirror.  Pinned nodes are zeroed below.
hfin = G.hs;  hfin(isinf(hfin)) = 0;
zm(:,:,1) = th(:,:,2) - 2*G.dh*(hfin/G.k).*th(:,:,1);

% floor ghost (below z = -d)
if strcmp(G.bot,'adiabatic')
    zp(:,:,nz) = th(:,:,nz-1);
end                                  % 'dirichlet' -> ghost = 0, already set

L = L + (zp + zm - 2*th)/dh2;

L(G.pin) = 0;
end
