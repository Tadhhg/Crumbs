function G = heat_grid(Tc, x, y, z)
%HEAT_GRID  Build the finite-difference grid, boundary-condition masks and
%           control-volume weights for the 3-D transient conduction solve.
%
%   G = HEAT_GRID(Tc, x, y, z)  with Tc from THERMAL_CONFIG and x, y, z the
%   node coordinate vectors [m].  z runs from 0 (surface, index 1) DOWN to
%   -d (floor, index end), so index 1 is always the gel/blood interface.
%
%   DISCRETISATION.  Node-centred second-order finite differences on a uniform
%   grid.  Nodes lie ON the boundaries, so a boundary node owns half a control
%   volume in the normal direction.  Boundary conditions are imposed with
%   ghost nodes:
%
%     Dirichlet  theta = 0     -> node removed from the unknowns (G.pin)
%     adiabatic  dtheta/dn = 0 -> ghost = mirror of the first interior node
%     Robin      -k dT/dz = h (T - T_blood) at z = 0.  With s = -z (depth),
%                k dtheta/ds|_0 = h theta_0, so the ghost node at s = -dh is
%                     theta_g = theta_2 - 2 dh h theta_1 / k
%                and the z-Laplacian at node 1 becomes
%                     (2 theta_2 - 2 theta_1 - 2 dh h theta_1 / k)/dh^2 .
%                This is IDENTICAL to the finite-VOLUME energy balance on the
%                half cell (the factor 1/2 divides out of both sides), so the
%                scheme is conservative as well as second-order.
%
%   G.W is the control-volume weight of each node: 1/2 for each direction in
%   which the node sits on a boundary (so 1/2 on a face, 1/4 on an edge, 1/8
%   at a corner).  It is needed for (a) the energy audit and (b) making the
%   steady-state operator symmetric, since W.*Laplacian is symmetric while the
%   bare Laplacian is not.
%
%   THE CATHETER.  Tc.cath_bc controls the surface condition over the strip
%   |y| <= D/2, cath_x0 <= x <= cath_L.  It is a cooled body lying on the
%   surface, never a heat source:
%     'transparent' leaves the surface condition unchanged (exact when the
%                   surface is Dirichlet -- pinning a plane that a cooled
%                   object rests on changes nothing);
%     'cooled'      applies h_cath there (cylinder in cross flow in series
%                   with the catheter wall);
%     'insulating'  makes the strip adiabatic -- the pessimistic bracket.

nx = numel(x); ny = numel(y); nz = numel(z);
dh = Tc.dh;

G = struct();
G.x = x(:); G.y = y(:); G.z = z(:);
G.nx = nx; G.ny = ny; G.nz = nz; G.dh = dh;
G.k = Tc.k_gel;  G.rho = Tc.rho_gel;  G.c = Tc.c_gel;
G.alpha = Tc.alpha_gel;
G.side = lower(Tc.side_bc);
G.bot  = lower(Tc.bot_bc);

% ---------------------------------------------------- surface film array --
hs = zeros(nx,ny);
switch lower(Tc.surf_bc)
    case 'dirichlet', hs(:) = Inf;
    case 'robin',     hs(:) = Tc.h_surf;
    case 'adiabatic', hs(:) = 0;
    otherwise, error('heat_grid:surf','unknown surf_bc ''%s''',Tc.surf_bc);
end

% catheter footprint on the surface
cx = (x(:) >= Tc.cath_x0 - 1e-12) & (x(:) <= Tc.cath_L + 1e-12);
cy = abs(y(:)) <= Tc.cath_D/2 + 1e-12;
G.cath_mask = cx & cy.';
switch lower(Tc.cath_bc)
    case 'transparent'   % leave hs alone
    case 'cooled',       hs(G.cath_mask) = Tc.h_cath;
    case 'insulating',   hs(G.cath_mask) = 0;
    otherwise, error('heat_grid:cath','unknown cath_bc ''%s''',Tc.cath_bc);
end
G.hs = hs;

% ------------------------------------------------------- Dirichlet mask ---
pin = false(nx,ny,nz);
pin(:,:,1) = isinf(hs);                         % pinned surface nodes
if strcmp(G.bot,'dirichlet'),  pin(:,:,nz) = true;  end
if strcmp(G.side,'dirichlet')
    pin(1,:,:) = true;  pin(nx,:,:) = true;
    pin(:,1,:) = true;  pin(:,ny,:) = true;
end
G.pin  = pin;
G.free = ~pin;

% ------------------------------------------------ control-volume weights --
wx = ones(nx,1); wx([1 end]) = 0.5;
wy = ones(1,ny); wy([1 end]) = 0.5;
wz = ones(1,1,nz); wz(1) = 0.5; wz(nz) = 0.5;
G.W = (wx*wy).*wz;

% ------------------------------------------- diagonal of the SPD operator --
% For every node the ghost-node stencil contributes -2/dh^2 to the diagonal
% in each of the three directions, whatever the boundary type, plus an extra
% -2*h/(k*dh) at a Robin surface node.  See the header.
hrob = zeros(nx,ny,nz);
hfin = hs;  hfin(isinf(hfin)) = 0;
hrob(:,:,1) = hfin;
G.diagA = G.W.*(6 + 2*dh*hrob/G.k)/dh^2;

% -------------------------------------------- explicit stability time step -
G.dt_stable = dh^2/(G.alpha*(6 + 2*dh*max(hfin(:))/G.k));
end
