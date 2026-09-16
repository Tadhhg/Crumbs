function R = run_heat_model(varargin)
%RUN_HEAT_MODEL  Baseline transient heating of the Phytagel slab under the
%                12 mm 2.45 GHz monopole, with 37 degC blood flowing over the
%                surface.  Produces results/heat_model.mat (used by T_gel.m)
%                and the figures.
%
%   run_heat_model                       baseline (Dirichlet surface)
%   run_heat_model('surf_bc','robin')    physical convective surface
%   run_heat_model('tag','robin', ...)   choose the output file suffix
%
%   Any THERMAL_CONFIG field can be overridden by name/value.
%
%   WHAT IS COMPUTED
%     Theta(r,t)  the temperature RISE above 37 degC for the reference
%                 delivered power (150 W).  Everything else is obtained from
%                 it by the exact scaling  theta = (Pdel/150) * Theta.
%     Theta_ss    the steady state (t -> infinity).
%   Runtime is about 12 min at the default 0.25 mm grid and t_end = 600 s.

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(here); addpath(fullfile(fileparts(root),'matlab'));
if ~exist(fullfile(root,'results'),'dir'), mkdir(fullfile(root,'results')); end
if ~exist(fullfile(root,'figures'),'dir'), mkdir(fullfile(root,'figures')); end

% pull out the file tag before it reaches thermal_config
tag = '';
kv = varargin;
it = find(strcmpi(kv(1:2:end),'tag'));
if ~isempty(it), tag = kv{2*it}; kv([2*it-1 2*it]) = []; end

Tc = thermal_config(kv{:});
fprintf('\n=========== THERMAL MODEL ===========\n');
print_thermal_config(Tc);

dh = Tc.dh;
xg = (Tc.xlim(1):dh:Tc.xlim(2)+1e-12).';
yg = (Tc.ylim(1):dh:Tc.ylim(2)+1e-12).';
zg = -(0:dh:Tc.d+1e-12).';

% ---------------- 1. the electromagnetic source ---------------------------
C = model_config('Pdel', Tc.Pdel_ref);
cfg = struct('f',C.f,'epsr',C.epsr,'d',C.d,'h',C.h,'L',C.L,'ke',C.ke, ...
             'Idist',C.Idist,'Nfft',C.Nfft,'dx',C.dx,'Pdel',C.Pdel, ...
             'zn',zg,'xlim',Tc.xlim,'ylim',Tc.ylim);
t0 = tic;  o = sar_volume(cfg);
% the EM solve always runs on the 0.25 mm FFT grid; subsample onto the
% thermal grid if that is coarser
step = round(dh/C.dx);
if step > 1
    o.q = o.q(1:step:end, 1:step:end, :);
    o.x = o.x(1:step:end);  o.y = o.y(1:step:end);
end
assert(numel(o.x)==numel(xg) && numel(o.y)==numel(yg) && size(o.q,3)==numel(zg), ...
       'EM block does not match the thermal grid');
fprintf('EM solve %.1f s\n', toc(t0));
fprintf('\n--- microwave power deposition (at %.0f W delivered) ---\n', Tc.Pdel_ref);
fprintf('  sigma_eff                    %.4f S/m\n', o.sigma1);
fprintf('  absorbed in the gel          %8.3f W  (%5.2f %%)\n', ...
        o.Pgel, 100*o.Pgel/(o.Pgel+o.Pup+o.Pdn));
fprintf('  inside the %.0f x %.0f x %.0f mm block  %8.3f W  (%5.2f %% of the gel total)\n', ...
        Tc.Wroi*1e3, Tc.Wroi*1e3, Tc.d*1e3, o.Pblock, 100*o.block_fraction);
fprintf('  peak q                       %.3e W/m^3  (SAR %.3e W/kg)\n', ...
        max(o.q(:)), max(o.q(:))/Tc.rho_gel);

% ---------------- 2. the transient solve ----------------------------------
G = heat_grid(Tc, xg, yg, zg);
probes = [ 6e-3 0 -0.25e-3; 6e-3 0 -1e-3; 6e-3 0 -2e-3; 6e-3 0 -5e-3
           6e-3 0 -10e-3;   6e-3 0 -20e-3; 0 0 -2e-3;   12e-3 0 -2e-3];
% the feed singularity: mask it out of PEAK REPORTING, and optionally cap the
% source there (see feed_regularise.m for the derivation of the point charge)
[o.q, finfo] = feed_regularise(o.q, xg, yg, zg, Tc);
pmask = finfo.mask_peak;
fprintf('\n--- the feed singularity ---\n');
fprintf('  a point charge Q = I_m/omega sits at the feed because the\n');
fprintf('  prescribed current is truncated there while still at maximum.\n');
fprintf('  excluded from PEAK REPORTING: %.1f mm sphere about (0,0,0), %d of %d nodes\n', ...
        1e3*Tc.r_feed_excl, finfo.n_inside, numel(pmask));
fprintf('  power inside that sphere: %.3f W of %.3f W in the block (%.2f %%)\n', ...
        finfo.P_inside_before, o.Pblock, 100*finfo.P_inside_before/o.Pblock);
fprintf('  feed_fix = ''%s''', Tc.feed_fix);
if strcmpi(Tc.feed_fix,'cap')
    fprintf(' -> capped at %.3e W/m^3, sphere now holds %.3f W\n', ...
            finfo.qcap, finfo.P_inside_after);
else
    fprintf(' -> source left exactly as the field model gives it\n');
end
fprintf('\n--- transient solve ---\n');
fprintf('  grid %d x %d x %d = %.2f M nodes, dh = %.3f mm\n', ...
        G.nx, G.ny, G.nz, G.nx*G.ny*G.nz/1e6, dh*1e3);
fprintf('  stability limit dt = %.4f s, using dt = %.4f s\n', ...
        G.dt_stable, Tc.cfl*G.dt_stable);
S = heat_solve(o.q, G, Tc, struct('probes',probes,'verbose',true, ...
                                  'peakmask',pmask));
fprintf('  %d steps in %.0f s wall time; energy-balance error %.2e\n', ...
        S.nstep, S.walltime, S.err_energy);

% ---------------- 3. the steady state -------------------------------------
t0 = tic;
[Th_ss, ssinfo] = heat_steady(o.q, G, 1e-9);
fprintf('  steady state: pcg %d iters, relres %.1e, %.0f s\n', ...
        ssinfo.iter, ssinfo.relres, toc(t0));

% ---------------- 4. diagnostics ------------------------------------------
iy0 = S.iy0;
Th_ss_y0 = squeeze(Th_ss(:,iy0,:));
Th_ss_m  = Th_ss;  Th_ss_m(~pmask) = -Inf;
[thmax_ss, li] = max(Th_ss_m(:));
[ax,ay,az] = ind2sub(size(Th_ss), li);
[thmax_ss_all, li2] = max(Th_ss(:));
[bx,by,bz] = ind2sub(size(Th_ss), li2);

fprintf('\n--- the scaling law ---\n');
fprintf('  theta(r,t; P) = (P_delivered / %.0f W) * Theta(r,t), exactly.\n', Tc.Pdel_ref);
fprintf('  T(r,t) = %.0f + eta * P_generator / %.0f * Theta(r,t),  eta = %.3f\n', ...
        Tc.T_blood, Tc.Pdel_ref, Tc.eta);

fprintf('\n--- peak temperature rise, at the reference %.0f W generator power ---\n', ...
        Tc.Pfwd_ref);
fprintf('      t [s]     Theta_max [K]     T_max [degC]   located at (x,y,depth) mm\n');
tshow = [0.01 0.1 0.5 1 2 5 10 30 60 120 300 600];
for q = 1:numel(tshow)
    [~,i1] = min(abs(S.t - tshow(q)));
    ii = S.ij_max(:,i1);
    fprintf('   %8.2f   %13.2f   %14.1f   (%5.1f, %5.1f, %5.2f)\n', S.t(i1), ...
        S.Th_max(i1), Tc.T_blood+S.Th_max(i1), 1e3*xg(ii(1)), 1e3*yg(ii(2)), ...
        -1e3*zg(ii(3)));
end
fprintf('   %8s   %13.2f   %14.1f   (%5.1f, %5.1f, %5.2f)\n', 'steady', ...
        thmax_ss, Tc.T_blood+thmax_ss, 1e3*xg(ax), 1e3*yg(ay), -1e3*zg(az));

% Is the hot spot on the plane of symmetry?  The deposition is symmetric in y
% and so is the solution, so physically it must be.  It can nevertheless be
% reported off-plane, because the feed-exclusion sphere cuts the y = 0 plane
% and pushes the masked maximum sideways to just outside it.
noff = nnz(S.ij_max(2,2:end) ~= iy0);
fprintf('  hot spot on y = 0: %d of %d output times', numel(S.t)-1-noff, numel(S.t)-1);
if noff > 0
    yo = yg(S.ij_max(2, S.ij_max(2,:) ~= iy0 & (1:numel(S.t)) > 1));
    fprintf('; the %d exceptions are displaced by up to %.2f mm and lie on the\n', ...
            noff, max(abs(yo))*1e3);
    fprintf('  surface of the feed-exclusion sphere, not somewhere physical\n');
else
    fprintf('\n');
end
fprintf('  distinct hot-spot locations over the whole history: %d\n', ...
        size(unique(S.ij_max(:,2:end).','rows'),1));

fprintf('\n--- the same numbers WITHOUT the feed exclusion, for comparison ---\n');
for q = [find(abs(S.t-1)==min(abs(S.t-1)),1), find(abs(S.t-60)==min(abs(S.t-60)),1)]
    ii = S.ij_max_all(:,q);
    fprintf('  t = %6.1f s : Theta = %10.2f K at (%5.2f, %5.2f, %5.2f) mm\n', ...
        S.t(q), S.Th_max_all(q), 1e3*xg(ii(1)), 1e3*yg(ii(2)), -1e3*zg(ii(3)));
end
fprintf('  steady     : Theta = %10.2f K at (%5.2f, %5.2f, %5.2f) mm\n', ...
        thmax_ss_all, 1e3*xg(bx), 1e3*yg(by), -1e3*zg(bz));
fprintf('  ratio unrestricted / reported at steady state: %.1f x\n', ...
        thmax_ss_all/thmax_ss);
fprintf('  Everything quoted elsewhere uses the RESTRICTED peak.\n');

fprintf('\n--- depth of the peak (the cooled-surface signature) ---\n');
[~,ix6] = min(abs(xg - 6e-3));
prof_ss = Th_ss_y0(ix6,:);
[pk, ipk] = max(prof_ss);
fprintf('  on the vertical x = 6 mm the steady peak is %.2f mm BELOW the\n', -1e3*zg(ipk));
fprintf('  surface, not at it: the blood pins z = 0 at %.0f degC while the\n', Tc.T_blood);
fprintf('  deposition peaks at the surface, so the maximum is pushed inward.\n');
fprintf('  Theta at the surface %.3f K, at the peak %.2f K.\n', prof_ss(1), pk);

fprintf('\n--- where the heat goes ---\n');
fprintf('  deposited in the block, total                %8.3f W\n', S.P_in_total);
fprintf('  removed at once at the nodes pinned to 37 degC %7.3f W\n', ...
        S.P_in_total - S.P_in);
fprintf('  left to heat the gel                         %8.3f W\n', S.P_in);
for q = [1 find(abs(S.t-1)==min(abs(S.t-1)),1) find(abs(S.t-10)==min(abs(S.t-10)),1) ...
         find(abs(S.t-60)==min(abs(S.t-60)),1) numel(S.t)]
    fprintf('  t = %7.2f s: stored %8.3f J, leaving %8.3f W (%5.1f %% of input)\n', ...
        S.t(q), S.U(q), S.P_loss(q), 100*S.P_loss(q)/S.P_in);
end
Lss = heat_lap(Th_ss,G);
Pl_ss = -G.k*G.dh^3*sum(G.W(:).*Lss(:));
fprintf('  steady state : leaving %8.3f W (%5.2f %% of input) -- must be 100 %%\n', ...
        Pl_ss, 100*Pl_ss/S.P_in);
if isfield(Tc,'conv')
    fprintf('  coolant bulk temperature rise at %.3f W removed: %.3f K\n', ...
            S.P_in, S.P_in*Tc.conv.dT_per_W);
    fprintf('  (the brief says to ignore it; this confirms that is safe)\n');
end

% ---------------- 5. the practical inversion ------------------------------
fprintf('\n--- generator power allowed, by exposure time and temperature limit ---\n');
fprintf('  P_gen [W] such that the PEAK temperature just reaches the limit\n');
fprintf('     t [s]      43 degC     50 degC     60 degC    100 degC\n');
lims = [43 50 60 100];
for tq = [1 5 10 30 60 120 300 600]
    [~,i1] = min(abs(S.t - tq));
    row = Tc.Pdel_ref*(lims - Tc.T_blood)/(Tc.eta*S.Th_max(i1));
    fprintf('   %7.0f  %10.3f  %10.3f  %10.3f  %10.3f\n', S.t(i1), row);
end
row = Tc.Pdel_ref*(lims - Tc.T_blood)/(Tc.eta*thmax_ss);
fprintf('   %7s  %10.3f  %10.3f  %10.3f  %10.3f\n', 'steady', row);

% ---------------- 6. save --------------------------------------------------
R = struct();
R.Tc = Tc;  R.x = xg;  R.y = yg;  R.z = zg;  R.iy0 = iy0;
R.q = single(o.q);
R.q_y0 = squeeze(o.q(:,iy0,:));
R.sigma = o.sigma1;  R.Pgel = o.Pgel;  R.Pblock = o.Pblock;
R.block_fraction = o.block_fraction;
R.t = S.t;  R.Th_y0 = S.Th_y0;  R.Th_max = S.Th_max;  R.ij_max = S.ij_max;
R.Th3 = S.Th3;  R.t3 = S.t3;  R.x3 = S.x3;  R.y3 = S.y3;  R.z3 = S.z3;
R.Th_probe = S.Th_probe;  R.probes = probes;
R.P_in = S.P_in;  R.P_loss = S.P_loss;  R.U = S.U;
R.err_energy = S.err_energy;  R.nstep = S.nstep;
R.Th_ss = single(Th_ss);  R.Th_ss_y0 = Th_ss_y0;  R.Th_max_ss = thmax_ss;
R.Th_max_all = S.Th_max_all;  R.ij_max_all = S.ij_max_all;
R.Th_max_ss_all = thmax_ss_all;  R.ij_ss = [ax;ay;az];
R.pmask_y0 = squeeze(pmask(:,iy0,:));  R.r_feed_excl = Tc.r_feed_excl;
R.P_in_total = S.P_in_total;  R.feed = rmfield(finfo,'mask_peak');
R.Pdel_ref = Tc.Pdel_ref;  R.Pfwd_ref = Tc.Pfwd_ref;  R.eta = Tc.eta;
R.T_blood = Tc.T_blood;
R.cath_mask = G.cath_mask;
fn = fullfile(root,'results',['heat_model' tag '.mat']);
save(fn,'-struct','R','-v7.3');
d = dir(fn);
fprintf('\nsaved %s (%.0f MB)\n', fn, d.bytes/1e6);

if isempty(tag)
    make_heat_figures(R, root);
end
end

% =========================================================================
function print_thermal_config(Tc)
fprintf('  gel: rho = %.0f kg/m^3, c = %.0f J/(kg K), k = %.3f W/(m K)\n', ...
        Tc.rho_gel, Tc.c_gel, Tc.k_gel);
fprintf('       alpha = %.4e m^2/s -> diffusion length sqrt(4 a t) = %.2f mm at 10 s\n', ...
        Tc.alpha_gel, 1e3*sqrt(4*Tc.alpha_gel*10));
fprintf('  coolant: %s at %.0f degC, %.2f L/min\n', Tc.coolant, Tc.T_blood, Tc.Qflow*60000);
if isfield(Tc,'conv')
    H = Tc.conv;
    fprintf('       U = %.4f m/s, Re_L = %.0f, Re_Dh = %.0f, Pr = %.2f\n', ...
            H.U, H.ReL, H.ReD, H.Pr);
    fprintf('       entrance lengths: hydrodynamic %.2f m, thermal %.2f m ', ...
            H.L_hydro, H.L_therm);
    fprintf('(both >> %.0f mm, so the flat-plate estimate is the right one)\n', ...
            H.params.Lplate*1e3);
    fprintf('       h: plate-lam %.0f | plate-turb %.0f | duct-DB %.0f | duct-lam %.0f W/(m^2 K)\n', ...
            H.h.plate_lam, H.h.plate_turb, H.h.duct_dittus, H.h.duct_lam);
    fprintf('       catheter: Re %.0f, h_out %.0f, wall R %.2e -> h_eff %.0f W/(m^2 K)\n', ...
            H.cath.Re, H.cath.h_out, H.cath.R_wall, H.cath.h_eff);
    fprintf('       cooling length k_gel/h = %.3f mm\n', 1e3*Tc.k_gel/Tc.h_surf);
end
fprintf('  boundary conditions: surface ''%s'', catheter ''%s'', floor ''%s'', sides ''%s''\n', ...
        Tc.surf_bc, Tc.cath_bc, Tc.bot_bc, Tc.side_bc);
fprintf('  ROI: x %.0f..%.0f mm, y %.0f..%.0f mm, depth 0..%.0f mm\n', ...
        Tc.xlim*1e3, Tc.ylim*1e3, Tc.d*1e3);
fprintf('  power reference: P_gen = %.0f W, eta = %.3f -> P_del = %.0f W\n', ...
        Tc.Pfwd_ref, Tc.eta, Tc.Pdel_ref);
end
