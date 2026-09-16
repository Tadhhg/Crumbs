function verify_heat_model(which)
%VERIFY_HEAT_MODEL  Validation suite for the thermal model.  RUN THIS FIRST.
%
%   verify_heat_model            all tests   (~10 min)
%   verify_heat_model(1:5)       a subset
%
%   TEST 1  3-D closed-form eigenfunction.  A source q0*sin*sin*sin in a box
%           with all six faces Dirichlet is an eigenfunction of the Laplacian,
%           so the transient has the exact solution
%               theta = q0/(k L) (1 - exp(-alpha L t)) sin sin sin,
%               L = pi^2 (1/Lx^2 + 1/Ly^2 + 1/Lz^2).
%           Checked against BOTH the continuum eigenvalue (shows second-order
%           convergence) and the discrete one (isolates the code from the
%           discretisation and should be exact to time-stepping error).
%
%   TEST 2  Robin boundary condition, 1-D limit.  Uniform source, adiabatic
%           sides and floor, convective surface.  The exact solution is the
%           quadratic  theta(s) = q0 d/h + (q0/k)(d s - s^2/2), which the
%           three-point stencil and the ghost-node Robin form reproduce
%           EXACTLY, so any error here is a coding error, not truncation.
%
%   TEST 3  Adiabatic box.  No boundary can lose heat, so theta must rise
%           uniformly as q0 t /(rho c).  Tests the mirror ghosts and the
%           energy bookkeeping to machine precision.
%
%   TEST 4  Symmetry and positive-definiteness of the steady operator, and
%           agreement between the steady solver and the long-time limit of the
%           explicit transient solver on the real problem.
%
%   TEST 5  Electromagnetic cross-validation.  The power deposition block
%           computed by SAR_VOLUME is compared, on the y = 0 plane, against
%           the independently stored solution of the field model
%           (results/field_model.mat and E_gel.m), and the delivered-power
%           normalisation is compared against the published field-model
%           budget.
%
%   TEST 6  Grid convergence of the full problem (dh = 1.0, 0.5, 0.25 mm).
%
%   TEST 7  Time-step convergence (CFL 0.9, 0.45, 0.225).
%
%   TEST 8  Linearity in delivered power -- the property the deliverable
%           function relies on.
%
%   TEST 9  Thermal dose (CEM43) against hand calculation.

if nargin < 1, which = 1:9; end
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(here); addpath(fullfile(fileparts(root),'matlab'));

fprintf('\n================ VALIDATION OF THE THERMAL MODEL ================\n');

% =========================================================== TEST 1 =======
if any(which==1)
fprintf('\n--- TEST 1: 3-D eigenfunction source, all faces Dirichlet ---\n');
for dh = [1e-3 0.5e-3 0.25e-3]
    Lb = 12e-3;
    Tc = thermal_config('dh',dh,'surf_bc','dirichlet','bot_bc','dirichlet', ...
                        'side_bc','dirichlet','h_surf',Inf,'h_cath',Inf, ...
                        'cath_bc','transparent');
    % y is centred on 0 so that the stored y = 0 plane cuts the PEAK of the
    % eigenfunction rather than one of its nodal faces
    xg = (0:dh:Lb).';  yg = (-Lb/2:dh:Lb/2).';  zg = -(0:dh:Lb).';
    G  = heat_grid(Tc, xg, yg, zg);
    q0 = 1e6;
    [Xg,Yg,Zg] = ndgrid(xg, yg + Lb/2, -zg);
    shape = sin(pi*Xg/Lb).*sin(pi*Yg/Lb).*sin(pi*Zg/Lb);
    q = q0*shape;

    Lam_c = 3*(pi/Lb)^2;                                    % continuum
    lam1  = 2*(1-cos(pi*dh/Lb))/dh^2;                       % discrete, 1-D
    Lam_d = 3*lam1;

    tt  = 4;                                                % [s]
    opt = struct('t_out',[0 tt],'verbose',false,'Nt3',2,'quiet',true);
    S   = heat_solve(q, G, Tc, opt);
    th  = double(S.Th_y0(:,:,end));                         % y = 0 plane
    iy0 = S.iy0;
    ana_c = q0/(G.k*Lam_c)*(1-exp(-G.alpha*Lam_c*tt))*squeeze(shape(:,iy0,:));
    ana_d = q0/(G.k*Lam_d)*(1-exp(-G.alpha*Lam_d*tt))*squeeze(shape(:,iy0,:));
    e_c = max(abs(th(:)-ana_c(:)))/max(ana_c(:));
    e_d = max(abs(th(:)-ana_d(:)))/max(ana_d(:));

    ths = heat_steady(q, G);
    ss_d = q0/(G.k*Lam_d)*shape;
    e_s = max(abs(ths(:)-ss_d(:)))/max(ss_d(:));
    fprintf(['  dh = %.3f mm | transient vs continuum %.3e | ' ...
             'vs discrete %.3e | steady vs discrete %.3e\n'], ...
            dh*1e3, e_c, e_d, e_s);
end
% The "discrete" column removes the SPATIAL error, so what is left is pure
% forward-Euler time-integration error.  It happens to fall as dh^2 as well,
% only because dt is tied to dh^2 by the stability limit.  Cutting the time
% step at fixed dh separates the two and must reduce it in proportion to dt:
dh = 1e-3;  Lb = 12e-3;  q0 = 1e6;  tt = 4;
for cfl = [0.9 0.1]
    Tc = thermal_config('dh',dh,'cfl',cfl,'surf_bc','dirichlet', ...
                        'bot_bc','dirichlet','side_bc','dirichlet', ...
                        'h_surf',Inf,'h_cath',Inf,'cath_bc','transparent');
    xg = (0:dh:Lb).';  yg = (-Lb/2:dh:Lb/2).';  zg = -(0:dh:Lb).';
    G  = heat_grid(Tc, xg, yg, zg);
    [Xg,Yg,Zg] = ndgrid(xg, yg + Lb/2, -zg);
    shape = sin(pi*Xg/Lb).*sin(pi*Yg/Lb).*sin(pi*Zg/Lb);
    Lam_d = 3*2*(1-cos(pi*dh/Lb))/dh^2;
    S = heat_solve(q0*shape, G, Tc, struct('t_out',[0 tt],'verbose',false, ...
                                           'Nt3',2,'quiet',true));
    th = double(S.Th_y0(:,:,end));
    ana_d = q0/(G.k*Lam_d)*(1-exp(-G.alpha*Lam_d*tt))*squeeze(shape(:,S.iy0,:));
    fprintf('  dh = 1.000 mm, CFL %.2f (dt = %.4f s): err vs discrete %.3e\n', ...
            cfl, S.dt_max, max(abs(th(:)-ana_d(:)))/max(ana_d(:)));
end
fprintf('  (continuum error must fall ~4x per halving of dh -- second order;\n');
fprintf('   discrete error is forward-Euler only, and falls in step with dt)\n');
end

% =========================================================== TEST 2 =======
if any(which==2)
fprintf('\n--- TEST 2: Robin surface, uniform source, 1-D exact quadratic ---\n');
dh = 0.5e-3;  dd = 20e-3;  hh = 1200;  q0 = 2e5;
Tc = thermal_config('dh',dh,'surf_bc','robin','h_surf',hh,'h_cath',hh, ...
                    'bot_bc','adiabatic','side_bc','adiabatic', ...
                    'cath_bc','transparent');
xg = (0:dh:5e-3).'; yg = xg; zg = -(0:dh:dd).';
G  = heat_grid(Tc, xg, yg, zg);
q  = q0*ones(numel(xg),numel(yg),numel(zg));
ths = heat_steady(q, G, 1e-13);
s   = -zg;
ana = q0*dd/hh + (q0/G.k)*(dd*s - s.^2/2);
got = squeeze(ths(round(end/2),round(end/2),:));
fprintf('  theta(surface)  exact %8.3f K   computed %8.3f K\n', ana(1), got(1));
fprintf('  theta(floor)    exact %8.3f K   computed %8.3f K\n', ana(end), got(end));
fprintf('  max rel err over the depth profile : %.3e\n', ...
        max(abs(got-ana))/max(ana));
% and the steady heat balance: everything deposited must leave through the top
Pin  = sum(G.W(:).*q(:))*dh^3;
Pout = hh*sum(sum(G.W(:,:,1)*2.*ths(:,:,1)))*dh^2;   % W(:,:,1) carries the 1/2
fprintf('  P_in = %.6f W ; P_out through the surface = %.6f W ; err %.2e\n', ...
        Pin, Pout, abs(Pin-Pout)/Pin);
end

% =========================================================== TEST 3 =======
if any(which==3)
fprintf('\n--- TEST 3: fully adiabatic box, uniform source ---\n');
dh = 1e-3;  q0 = 5e5;  tt = 30;
Tc = thermal_config('dh',dh,'surf_bc','adiabatic','h_surf',0,'h_cath',0, ...
                    'bot_bc','adiabatic','side_bc','adiabatic', ...
                    'cath_bc','transparent');
xg = (0:dh:10e-3).'; yg = xg; zg = -(0:dh:10e-3).';
G  = heat_grid(Tc, xg, yg, zg);
q  = q0*ones(numel(xg),numel(yg),numel(zg));
S  = heat_solve(q, G, Tc, struct('t_out',[0 tt],'verbose',false,'Nt3',2,'quiet',true));
ana = q0*tt/(G.rho*G.c);
fprintf('  theta exact %10.6f K  computed %10.6f K  rel err %.3e\n', ...
        ana, S.Th_max(end), abs(S.Th_max(end)-ana)/ana);
fprintf('  boundary loss (must be 0): %.3e W ; energy audit %.3e\n', ...
        S.P_loss(end), S.err_energy);
end

% =========================================================== TEST 4 =======
if any(which==4)
fprintf('\n--- TEST 4: operator symmetry and steady vs long transient ---\n');
dh = 1e-3;
Tc = thermal_config('dh',dh);
[xg,yg,zg] = roi_nodes(Tc);
G = heat_grid(Tc, xg, yg, zg);
nf = nnz(G.free);
rng(7);
u = randn(nf,1); v = randn(nf,1);
Au = applyA(u,G);  Av = applyA(v,G);
fprintf('  <Au,v> - <u,Av> = %.3e (relative %.3e)\n', ...
        v'*Au - u'*Av, abs(v'*Au - u'*Av)/abs(v'*Au));
fprintf('  <u,Au> = %.4e  (must be > 0)\n', u'*Au);

o  = build_source(Tc, xg, yg, zg);
q  = o.q;
ths = heat_steady(q, G, 1e-11);
tlong = 20000;
S = heat_solve(q, G, Tc, struct('t_out',[0 tlong],'verbose',false, ...
                                'Nt3',2,'quiet',true));
fprintf('  steady peak %.3f K vs transient at t = %g s peak %.3f K (%.2f %%)\n', ...
        max(ths(:)), tlong, S.Th_max(end), ...
        100*abs(max(ths(:))-S.Th_max(end))/max(ths(:)));
% Power deposited AT a node pinned to 37 degC is removed by the sink and
% never enters the solution, so the balance is over the free nodes.
Pin_all  = sum(G.W(:).*q(:))*G.dh^3;
Pin_free = sum(G.W(G.free).*q(G.free))*G.dh^3;
Lss = heat_lap(ths,G);
Pl  = -G.k*G.dh^3*sum(G.W(:).*Lss(:));
fprintf('  deposited in the block %.4f W, of which %.4f W lands on pinned nodes\n', ...
        Pin_all, Pin_all-Pin_free);
fprintf('  steady-state balance: P_in(free) %.4f W  P_loss %.4f W  err %.2e\n', ...
        Pin_free, Pl, abs(Pin_free-Pl)/Pin_free);
end

% =========================================================== TEST 5 =======
if any(which==5)
fprintf('\n--- TEST 5: electromagnetic cross-validation ---\n');
Tc = thermal_config();
[xg,yg,zg] = roi_nodes(Tc);
o = build_source(Tc, xg, yg, zg);
F = load(fullfile(fileparts(root),'results','field_model.mat'));
fprintf('  delivered power        this run %.3f W   field model %.3f W\n', ...
        o.Pgel+o.Pup+o.Pdn, F.Pgel+F.Pup+F.Pdn);
fprintf('  absorbed in the gel    this run %.3f W (%.2f %%)  field model %.3f W (%.2f %%)\n', ...
        o.Pgel, 100*o.Pgel/(o.Pgel+o.Pup+o.Pdn), F.Pgel, ...
        100*F.Pgel/(F.Pgel+F.Pup+F.Pdn));
fprintf('  current scale I_m      this run %.4f A     field model %.4f A  (%.3f %%)\n', ...
        o.scale, F.scale, 100*abs(o.scale-F.scale)/F.scale);
fprintf('  sigma_eff              this run %.5f       field model %.5f S/m\n', ...
        o.sigma1, F.sigma);
% point-by-point on the y = 0 plane, 1..29 mm depth.  Both grids are indexed
% by DEPTH s = -z so that the interpolation coordinates increase.
xs = [0 3 6 9 12 18].'*1e-3;  ss = (1:29).'*1e-3;
[XX,SS] = ndgrid(xs, ss);
Pv_ref = interp2(F.x, -F.z, F.Qv, XX, SS, 'spline');
iy0 = find(abs(yg)<1e-12,1);
qp  = squeeze(o.q(:,iy0,:));
Pv_new = interpn(xg, -zg, qp, XX, SS, 'spline');
rel = abs(Pv_new-Pv_ref)./Pv_ref;
fprintf('  q on y = 0, %d points 1-29 mm deep: max rel diff %.3e, rms %.3e\n', ...
        numel(rel), max(rel(:)), sqrt(mean(rel(:).^2)));
[~,Se] = E_gel(6e-3, -5e-3);
qpt = interpn(xg, -zg, qp, 6e-3, 5e-3, 'spline');
fprintf('  E_gel(6mm,-5mm).Pv = %.4e W/m^3 ; this model %.4e W/m^3 (%.2f %%)\n', ...
        Se.Pv, qpt, 100*abs(Se.Pv-qpt)/Se.Pv);
fprintf('  power inside the 30x30x30 mm block: %.3f W of %.3f W in the gel (%.2f %%)\n', ...
        o.Pblock, o.Pgel, 100*o.block_fraction);
end

% =========================================================== TEST 6 =======
if any(which==6)
fprintf('\n--- TEST 6: grid convergence of the full problem (t = 30 s) ---\n');
fprintf('   Two columns behave quite differently and that is the point: the\n');
fprintf('   UNRESTRICTED peak sits on the feed singularity and diverges under\n');
fprintf('   refinement, while a probe 6 mm away converges at second order.\n');
fprintf('     dh[mm]   nodes        peak(all)  peak(feed excl)  theta(6,0,-2mm)   P_in[W]\n');
for dh = [1e-3 0.5e-3 0.25e-3]
    Tc = thermal_config('dh',dh);
    [xg,yg,zg] = roi_nodes(Tc);
    o = build_source(Tc, xg, yg, zg);
    G = heat_grid(Tc, xg, yg, zg);
    [~,fi] = feed_regularise(o.q, xg, yg, zg, Tc);
    pr = [6e-3 0 -2e-3];
    S = heat_solve(o.q, G, Tc, struct('t_out',[0 30],'verbose',false, ...
                    'probes',pr,'Nt3',2,'quiet',true,'peakmask',fi.mask_peak));
    fprintf('   %6.3f  %4dx%3dx%3d  %11.2f %14.2f %16.3f %9.3f\n', dh*1e3, ...
            G.nx,G.ny,G.nz, S.Th_max_all(end), S.Th_max(end), ...
            S.Th_probe(1,end), S.P_in);
end
fprintf('   Repeating with the feed source capped (feed_fix = ''cap''):\n');
for dh = [1e-3 0.5e-3 0.25e-3]
    Tc = thermal_config('dh',dh,'feed_fix','cap');
    [xg,yg,zg] = roi_nodes(Tc);
    o = build_source(Tc, xg, yg, zg);
    G = heat_grid(Tc, xg, yg, zg);
    [qc,fi] = feed_regularise(o.q, xg, yg, zg, Tc);
    pr = [6e-3 0 -2e-3];
    S = heat_solve(qc, G, Tc, struct('t_out',[0 30],'verbose',false, ...
                    'probes',pr,'Nt3',2,'quiet',true,'peakmask',fi.mask_peak));
    fprintf('   %6.3f  %4dx%3dx%3d  %11.2f %14.2f %16.3f %9.3f\n', dh*1e3, ...
            G.nx,G.ny,G.nz, S.Th_max_all(end), S.Th_max(end), ...
            S.Th_probe(1,end), S.P_in);
end
end

% =========================================================== TEST 7 =======
if any(which==7)
fprintf('\n--- TEST 7: time-step convergence (dh = 0.5 mm, t = 30 s) ---\n');
Tc0 = thermal_config('dh',0.5e-3);
[xg,yg,zg] = roi_nodes(Tc0);
o = build_source(Tc0, xg, yg, zg);
tt = unique([0 logspace(-3, log10(30), 200)]);
for cfl = [0.9 0.45 0.225]
    Tc = thermal_config('dh',0.5e-3,'cfl',cfl);
    G  = heat_grid(Tc, xg, yg, zg);
    [~,fi] = feed_regularise(o.q, xg, yg, zg, Tc);
    S  = heat_solve(o.q, G, Tc, struct('t_out',tt,'verbose',false, ...
                    'Nt3',2,'quiet',true,'peakmask',fi.mask_peak));
    fprintf('   CFL %.3f  dt = %.4f s  %6d steps  theta_max = %.4f K  energy err %.2e\n', ...
            cfl, S.dt_max, S.nstep, S.Th_max(end), S.err_energy);
end
fprintf('   (the energy error is the closure of dU/dt = P_in - P_loss integrated\n');
fprintf('    over the run; it is limited by the trapezoid used on P_loss between\n');
fprintf('    output times, not by the solver, so it falls with more output times)\n');
end

% =========================================================== TEST 8 =======
if any(which==8)
fprintf('\n--- TEST 8: exact linearity in delivered power ---\n');
Tc = thermal_config('dh',1e-3);
[xg,yg,zg] = roi_nodes(Tc);
o = build_source(Tc, xg, yg, zg);
G = heat_grid(Tc, xg, yg, zg);
S1 = heat_solve(1.0*o.q, G, Tc, struct('t_out',[0 5 60],'verbose',false,'Nt3',2,'quiet',true));
S2 = heat_solve(3.7*o.q, G, Tc, struct('t_out',[0 5 60],'verbose',false,'Nt3',2,'quiet',true));
fprintf('   t =  5 s : 3.7*theta(P) = %.6f K ; theta(3.7P) = %.6f K ; rel %.2e\n', ...
        3.7*S1.Th_max(2), S2.Th_max(2), abs(3.7*S1.Th_max(2)-S2.Th_max(2))/S2.Th_max(2));
fprintf('   t = 60 s : 3.7*theta(P) = %.6f K ; theta(3.7P) = %.6f K ; rel %.2e\n', ...
        3.7*S1.Th_max(3), S2.Th_max(3), abs(3.7*S1.Th_max(3)-S2.Th_max(3))/S2.Th_max(3));
fprintf('   (this is exact by construction -- it is the property T_gel.m uses\n');
fprintf('    to turn one solve into a function of power; the check confirms\n');
fprintf('    no power-dependent term has crept into the implementation)\n');
end

% =========================================================== TEST 9 =======
if any(which==9)
fprintf('\n--- TEST 9: thermal dose against hand calculation ---\n');
t = linspace(0,600,1001);
chk = @(nm,got,want) fprintf('  %-34s %12.5g  expected %12.5g  rel %.1e\n', ...
        nm, got, want, abs(got-want)/abs(want));
chk('43 degC held for 10 min',  cem43(t, 43*ones(size(t))), 10);
chk('33 degC held for 10 min',  cem43(t, 33*ones(size(t))), 10*0.25^10);
t2 = linspace(0,60,201);
chk('53 degC held for 1 min',   cem43(t2, 53*ones(size(t2))), 0.5^(43-53));
chk('60 degC held for 1 min',   cem43(t2, 60*ones(size(t2))), 0.5^(43-60));
fprintf('  R = 0.5 above 43 degC and 0.25 below, so each degree above 43\n');
fprintf('  doubles the dose and each degree below quarters it.\n');
end

fprintf('\n================ END VALIDATION ================\n\n');
end

% =========================================================================
function [xg,yg,zg] = roi_nodes(Tc)
dh = Tc.dh;
xg = (Tc.xlim(1):dh:Tc.xlim(2)+1e-12).';
yg = (Tc.ylim(1):dh:Tc.ylim(2)+1e-12).';
zg = -(0:dh:Tc.d+1e-12).';
end

function o = build_source(Tc, xg, yg, zg)
%BUILD_SOURCE  q [W/m^3] on the thermal nodes, at Pdel_ref.
persistent cache
key = sprintf('%.9g_%.9g_%d', Tc.dh, Tc.Pdel_ref, numel(zg));
if ~isempty(cache) && isfield(cache, matlab.lang.makeValidName(key))
    o = cache.(matlab.lang.makeValidName(key));  return
end
C = model_config('Pdel', Tc.Pdel_ref);
% the EM solve always runs on the 0.25 mm FFT grid; only the depths and the
% x/y subsampling follow the thermal grid
step = round(Tc.dh/C.dx);
cfg = struct('f',C.f,'epsr',C.epsr,'d',C.d,'h',C.h,'L',C.L,'ke',C.ke, ...
             'Idist',C.Idist,'Nfft',C.Nfft,'dx',C.dx,'Pdel',C.Pdel, ...
             'zn',zg,'xlim',Tc.xlim,'ylim',Tc.ylim);
o = sar_volume(cfg);
o.q = o.q(1:step:end, 1:step:end, :);
o.x = o.x(1:step:end);  o.y = o.y(1:step:end);
if numel(o.x) ~= numel(xg) || numel(o.y) ~= numel(yg)
    error('build_source:grid','EM block %dx%d does not match thermal grid %dx%d', ...
          numel(o.x), numel(o.y), numel(xg), numel(yg));
end
cache.(matlab.lang.makeValidName(key)) = o;
end

function out = applyA(v, G)
th = zeros(G.nx,G.ny,G.nz);
th(G.free) = v;
L = heat_lap(th, G);
out = -G.W(G.free).*L(G.free);
end
