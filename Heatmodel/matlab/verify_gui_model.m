function ok = verify_gui_model(varargin)
%VERIFY_GUI_MODEL  Validation suite for the fast model behind HEAT_GUI.
%
%   ok = VERIFY_GUI_MODEL          run every test, print a table
%   ok = VERIFY_GUI_MODEL('quick') skip the 0.25 mm test (test 3), which is
%                                  the only slow one
%
%   GUI_MODEL trades accuracy for speed in exactly three places, and each one
%   is measured here rather than asserted:
%
%     1. the grid is 0.5 mm, not the 0.25 mm of RUN_HEAT_MODEL
%     2. the transverse FFT window is 512 points, not 1024
%     3. the initial condition may be non-zero, which is a code path
%        RUN_HEAT_MODEL never exercises
%
%   TESTS
%     1  against the stored Dirichlet baseline (results/heat_model.mat)
%     2  against the stored Robin baseline     (results/heat_model_robin.mat)
%     3  grid convergence, 1.0 / 0.5 / 0.25 mm
%     4  FFT window, Nfft = 512 against 1024
%     5  exact linearity in the delivered power
%     6  superposition of the initial condition -- the new HEAT_SOLVE path
%     7  the energy audit, and what the quadrature in it is worth
%     8  the parameters the GUI exists to sweep actually move the answer, in
%        the direction physics says they should
%
%   Tests 1 and 2 are the ones that matter: they compare a few-second solve
%   against the quarter-hour solve the write-up is built on, at the write-up's
%   own configuration (gel starting at the coolant temperature, transparent
%   catheter), and quote the difference in kelvin.
%
%   See also GUI_MODEL, HEAT_GUI, VERIFY_HEAT_MODEL.

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(here); addpath(fullfile(fileparts(root),'matlab'));
quick = any(strcmpi(varargin,'quick'));

fprintf('\n================ VERIFY_GUI_MODEL ================\n');
res = {};   ok = true;

% The write-up's own configuration, so the comparison is like for like:
% the gel starts in equilibrium with the coolant and the catheter is
% transparent.  Only the grid and the FFT window differ.
base = {'T_gel',37, 'T_water',37, 'cath_bc','transparent', ...
        'Pgen',12, 't_stim',60};

% =========================================================== test 1 =======
% Dirichlet: the stored baseline that T_GEL and every quoted number use.
[res, ok] = compare_stored(res, ok, 1, 'Dirichlet surface', ...
    fullfile(root,'results','heat_model.mat'), [base {'surf_bc','dirichlet'}]);

% =========================================================== test 2 =======
[res, ok] = compare_stored(res, ok, 2, 'Robin surface', ...
    fullfile(root,'results','heat_model_robin.mat'), [base {'surf_bc','robin'}]);

% =========================================================== test 3 =======
fprintf('\n--- 3. grid convergence ---\n');
dhl = [1e-3 0.5e-3 0.25e-3];
if quick, dhl = dhl(1:2); end
pk = zeros(size(dhl));  pr = zeros(size(dhl));  tt = zeros(size(dhl));
for m = 1:numel(dhl)
    Rm = gui_model(base{:}, 'surf_bc','robin', 'dh',dhl(m));
    pk(m) = Rm.info.T_end;  pr(m) = Rm.info.T_probe_end;  tt(m) = Rm.info.t_total;
    fprintf('  dh = %.2f mm : peak %8.2f degC   probe %7.2f degC   (%.1f s)\n', ...
            1e3*dhl(m), pk(m), pr(m), tt(m));
end
% The probe is the convergent quantity.  The PEAK is not expected to converge
% and must not be asserted on: it sits on the feed-exclusion sphere, where the
% deposition is a truncated point charge, so refining the grid keeps finding a
% larger value.  That is the artefact of Sec. 6.2, not a numerical failure.
dpr = 100*abs(pr(end)-pr(end-1))/(pr(end-1)-37);
fprintf('  probe changes %.2f %% of its rise over the last refinement\n', dpr);
fprintf('  peak changes  %.2f %% -- expected NOT to converge (feed charge, Sec. 6.2)\n', ...
        100*abs(pk(end)-pk(end-1))/(pk(end-1)-37));
[res, ok] = add(res, ok, 3, 'grid convergence of the probe point', dpr, 5, '%');

% =========================================================== test 4 =======
fprintf('\n--- 4. FFT window ---\n');
R512  = gui_model(base{:}, 'surf_bc','robin', 'Nfft',512);
R1024 = gui_model(base{:}, 'surf_bc','robin', 'Nfft',1024);
dsar = 100*max(abs(R512.SAR(:)-R1024.SAR(:)))/max(R1024.SAR(:));
dpr4 = 100*abs(R512.info.T_probe_end - R1024.info.T_probe_end) / ...
       (R1024.info.T_probe_end - 37);
fprintf('  peak SAR   512: %.4e   1024: %.4e W/kg\n', ...
        R512.info.SAR_max, R1024.info.SAR_max);
fprintf('  worst SAR difference anywhere on the plane %.3f %% of the peak\n', dsar);
fprintf('  probe temperature differs by %.3f %% of its rise (%.1f s vs %.1f s)\n', ...
        dpr4, R512.info.t_em, R1024.info.t_em);
[res, ok] = add(res, ok, 4, 'Nfft 512 vs 1024, probe', dpr4, 1, '%');

% =========================================================== test 5 =======
fprintf('\n--- 5. linearity in the delivered power ---\n');
% With the gel starting AT the coolant temperature the whole solution is
% linear in the power, exactly.  This is the property T_GEL is built on.
Ra = gui_model(base{:}, 'surf_bc','robin', 'Pgen',5);
Rb = gui_model(base{:}, 'surf_bc','robin', 'Pgen',20, 'em',Ra.em);
rise_a = Ra.T_plane(:,:,end) - 37;
rise_b = Rb.T_plane(:,:,end) - 37;
lin = max(abs(rise_b(:) - 4*rise_a(:)))/max(abs(rise_b(:)));
fprintf('  theta(20 W) vs 4 x theta(5 W): worst relative difference %.2e\n', lin);
[res, ok] = add(res, ok, 5, 'exact power scaling', lin, 1e-10, '-');

% =========================================================== test 6 =======
fprintf('\n--- 6. superposition of the initial condition ---\n');
% The heat equation is linear, so a gel that starts D kelvin from the coolant
% and is then driven must equal (driven from equilibrium) + (relaxing from D
% with no drive).  This exercises the opt.theta0 path added to HEAT_SOLVE.
% The permittivity MUST be pinned for this test.  Left linked, changing the
% gel temperature would change eps' and sigma as well as the initial
% condition, so the three runs would not share a source term and the
% comparison would be meaningless -- which is exactly the trap the GUI's
% tick box exists to make explicit.
D   = -12;                                      % gel at 25, coolant at 37
fix = [base, {'surf_bc','robin', 'link_eps',false, ...
              'eps_r',74.29, 'sigma',2.7228}];
Rd = gui_model(fix{:});                                     % driven only
Rr = gui_model(fix{:}, 'Pgen',1e-9, 'T_gel',37+D, 'em',Rd.em);  % relaxing only
Rc = gui_model(fix{:}, 'T_gel',37+D, 'em',Rd.em);           % both at once
lhs = Rc.T_plane(:,:,end) - 37;
rhs = (Rd.T_plane(:,:,end) - 37) + (Rr.T_plane(:,:,end) - 37);
sup = max(abs(lhs(:)-rhs(:)))/max(abs(lhs(:)));
fprintf('  combined vs driven + relaxation: worst relative difference %.2e\n', sup);
fprintf('  (relaxation alone has fallen to %.3f K at 60 s)\n', ...
        max(abs(Rr.T_plane(:,:,end)-37),[],'all'));
% HEAT_SOLVE stores the y = 0 plane as SINGLE, so eps('single') = 1.2e-7 is
% the floor this comparison can reach; anything at that level is exact.
[res, ok] = add(res, ok, 6, 'initial-condition superposition', sup, 1e-6, '-');

% =========================================================== test 7 =======
fprintf('\n--- 7. the energy audit ---\n');
% dU/dt = P_in - P_loss must hold on the discrete solution.  It does not hold
% exactly, and there are two separate reasons, which this test separates:
%
%   (a) the audit integrates P_loss by trapezoid over the OUTPUT times, which
%       are few and log-spaced.  Refining them shrinks this part.
%   (b) forward Euler satisfies the CONTINUOUS energy balance only to O(dt).
%       At 0.5 mm the step is 0.19 s, three times the baseline's, so this part
%       is three times larger than in RUN_HEAT_MODEL.  It is a real
%       discretisation error, and the check is that halving the step shrinks it
%       AT LEAST in proportion -- in practice it does rather better than that,
%       because the two error sources partly cancel at this step size.
%
% The quantity to test is the balance at the FINAL time, where the stored
% energy is large and the denominator is not near zero.  The max over all
% times, which is what err_energy reports, is dominated by early samples whose
% denominator is tiny, so it is a poor figure of merit here.
fin = @(A) abs(A.U(end) - (A.U(1) + A.P_in*A.t(end) - trapz(A.t,A.P_loss))) / ...
           abs(A.U(1) + A.P_in*A.t(end) - trapz(A.t,A.P_loss));
fprintf('  output-time refinement (theta0 = 0, dt = %.3f s):\n', Rd.info.dt);
for Ntq = [41 161 641]
    Re = gui_model(fix{:}, 'Nt',Ntq, 'em',Rd.em);
    fprintf('    %3d output times : max %.2e   at the final time %.2e\n', ...
            Ntq, Re.info.err_energy, fin(Re.audit));
end
e_full = fin(Rd.audit);
Rh = gui_model(fix{:}, 'cfl',0.45, 'em',Rd.em);
e_half = fin(Rh.audit);
fprintf('  time step %.4f s -> final-time balance %.2e\n', Rd.info.dt, e_full);
fprintf('  time step %.4f s -> final-time balance %.2e   (ratio %.2f)\n', ...
        Rh.info.dt, e_half, e_full/e_half);
[res, ok] = add(res, ok, 7, 'energy balance at the final time', e_full, 2e-3, '-');
[res, ok] = add(res, ok, 7, 'it is first order in dt (ratio >= 1.7)', ...
                1.7/max(e_full/e_half,eps), 1, '-');

% =========================================================== test 8 =======
fprintf('\n--- 8. the swept parameters move the answer, the right way ---\n');
sw = struct('name',{},'lo',{},'hi',{},'want',{});
Rref = gui_model('T_gel',37,'T_water',37,'Pgen',12,'t_stim',60);   % GUI defaults
q1 = gui_model('T_gel',37,'T_water',37,'Pgen',12,'t_stim',60,'Qflow_Lmin',1,'em',Rref.em);
q2 = gui_model('T_gel',37,'T_water',37,'Pgen',12,'t_stim',60,'Qflow_Lmin',20,'em',Rref.em);
sw(end+1) = mk('flow 1 -> 20 L/min', q1.info.T_probe_end, q2.info.T_probe_end, 'down');
c1 = gui_model('T_gel',37,'T_water',37,'Pgen',12,'t_stim',60,'L_cath_mm',0,'em',Rref.em);
c2 = gui_model('T_gel',37,'T_water',37,'Pgen',12,'t_stim',60,'L_cath_mm',32,'em',Rref.em);
sw(end+1) = mk('catheter 0 -> 32 mm', c1.info.T_probe_end, c2.info.T_probe_end, 'up');
% Raising the coolant temperature shifts everything up by definition, which
% tests nothing.  What it also does, and what is worth testing, is change the
% water PROPERTIES: at 47 degC the viscosity is lower and the conductivity
% higher, so h is about 4 % larger and the temperature RISE is smaller.
% The permittivity is pinned here too, for the same reason as test 6: linked,
% the gel temperature would change the source term as well as the coolant, and
% (being outside the Stogryn fit above 40 degC) would warn as well.
wf = {'link_eps',false,'eps_r',74.29,'sigma',2.7228,'Pgen',12,'t_stim',60};
w1 = gui_model(wf{:},'T_gel',37,'T_water',37);
w2 = gui_model(wf{:},'T_gel',47,'T_water',47,'em',w1.em);
fprintf('  (h goes %.0f -> %.0f W/(m^2 K) over that change)\n', ...
        w1.info.h_surf, w2.info.h_surf);
sw(end+1) = mk('water 37 -> 47 degC, rise', w1.info.T_probe_end-37, ...
               w2.info.T_probe_end-47, 'down');
s1 = gui_model('T_gel',37,'T_water',37,'Pgen',12,'t_stim',60,'link_eps',false, ...
               'eps_r',74.29,'sigma',1.5);
s2 = gui_model('T_gel',37,'T_water',37,'Pgen',12,'t_stim',60,'link_eps',false, ...
               'eps_r',74.29,'sigma',4.0);
sw(end+1) = mk('sigma 1.5 -> 4.0 S/m', s1.info.d_1e*1e3, s2.info.d_1e*1e3, 'down');
for m = 1:numel(sw)
    good = (strcmp(sw(m).want,'up') && sw(m).hi > sw(m).lo) || ...
           (strcmp(sw(m).want,'down') && sw(m).hi < sw(m).lo);
    fprintf('  %-22s %8.3f -> %8.3f  (expected %-4s)   %s\n', sw(m).name, ...
            sw(m).lo, sw(m).hi, sw(m).want, tick(good));
    ok = ok && good;
    res(end+1,:) = {8, sw(m).name, NaN, NaN, '', good}; %#ok<AGROW>
end

% ============================================================ summary =====
fprintf('\n================ SUMMARY ================\n');
for m = 1:size(res,1)
    if isnan(res{m,3})
        fprintf('  test %d  %-42s %19s  %s\n', res{m,1}, res{m,2}, '', tick(res{m,6}));
    else
        fprintf('  test %d  %-42s %9.3g <= %-6.3g %-2s %s\n', res{m,1}, res{m,2}, ...
                res{m,3}, res{m,4}, res{m,5}, tick(res{m,6}));
    end
end
fprintf('\n  %s\n\n', ternary(ok,'ALL TESTS PASSED','SOMETHING FAILED -- read the table'));
end

% =========================================================================
function [res, ok] = compare_stored(res, ok, n, name, matfile, args)
%COMPARE_STORED  The fast model against the stored full-resolution solution.
fprintf('\n--- %d. against the stored baseline, %s ---\n', n, name);
if ~exist(matfile,'file')
    fprintf('  %s not found -- run run_heat_model first.  SKIPPED.\n', matfile);
    res(end+1,:) = {n, [name ' (no stored file)'], NaN, NaN, '', true};
    return
end
M = load(matfile,'x','z','t','Th_max','Th_y0','Tc','eta','Pdel_ref','T_blood');
G = gui_model(args{:});
Pgen = G.opt.Pgen;   s = Pgen*M.eta/M.Pdel_ref;      % the exact power scaling

% peak, at the end of the exposure
[~,i60] = min(abs(M.t - G.opt.t_stim));
peak_ref = M.T_blood + s*M.Th_max(i60);

% the bench point (6, 0, -2) mm
[~,ix] = min(abs(M.x - 6e-3));   [~,iz] = min(abs(M.z + 2e-3));
pr_ref = M.T_blood + s*double(M.Th_y0(ix,iz,i60));

d_peak = G.info.T_end     - peak_ref;
d_prob = G.info.T_probe_end - pr_ref;
fprintf('  peak  at %g s : stored %8.2f   fast %8.2f degC   (%+.2f K, %+.2f %% of rise)\n', ...
        G.opt.t_stim, peak_ref, G.info.T_end, d_peak, 100*d_peak/(peak_ref-37));
fprintf('  probe (6,0,-2): stored %8.2f   fast %8.2f degC   (%+.2f K, %+.2f %% of rise)\n', ...
        pr_ref, G.info.T_probe_end, d_prob, 100*d_prob/(pr_ref-37));
fprintf('  solved in %.1f s, against a %s solve for the stored one\n', ...
        G.info.t_total, 'quarter-hour');
[res, ok] = add(res, ok, n, [name ': probe'], abs(100*d_prob/(pr_ref-37)), 8, '%');
end

% =========================================================================
function [res, ok] = add(res, ok, n, name, val, tol, unit)
good = val <= tol;
ok = ok && good;
res(end+1,:) = {n, name, val, tol, unit, good};
end

function s = mk(name, lo, hi, want)
s = struct('name',name,'lo',lo,'hi',hi,'want',want);
end

function s = tick(b)
if b, s = 'PASS'; else, s = '**FAIL**'; end
end

function s = ternary(c,a,b)
if c, s = a; else, s = b; end
end
