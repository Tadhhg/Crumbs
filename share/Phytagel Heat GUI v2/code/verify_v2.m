function verify_v2()
%VERIFY_V2  Checks on the bath-above-the-gel change, before any result from it
%           is believed.
%
%   1  tangential E is continuous across z = 0   -- pins the TM sign S_TM
%   2  'air' above reproduces v1 EXACTLY         -- nothing else moved
%   3  the cut height Zw does not matter         -- the construction as a whole
%   4  the power budget closes                   -- fractions sum to 1
%   5  uniform medium, wire on z = 0, is up/down symmetric
%   6  the braid return arm removes the feed point charge, and L_sheath = 0
%      still reproduces the truncated current bit for bit
%   7  duty cycling delivers exactly the energy it claims to
%
%   Test 5 is the only one that puts the NEW bath field expressions and the
%   OLD gel ones on the same footing, so it is the one that would catch an
%   error the continuity test at a single plane could miss.
%
%   Test 6b is the one that matters for the feed, and it tests the MECHANISM
%   rather than a symptom.  A first version tried to catch the singularity by
%   refining the grid and watching the peak deposition run away.  It does not
%   run away: the gel starts 0.5 mm below the wire, so no amount of refinement
%   ever samples closer than that to the point charge and the peak converges --
%   to a value 56 times too large, but it converges.  The test would have
%   passed a model with the charge still in it.  What is unambiguous is the
%   spectral tail: a jump of dI at the feed makes kx*Itilde(kx) -> dI, and a
%   continuous I makes it -> 0.  No grid, no threshold to argue about.
%
%   Test 6c is the consequence that matters for the deliverable: with the
%   charge gone the PEAK TEMPERATURE stops depending on the grid.  The
%   truncated model's peak moves by tens of degrees between grids and is in any
%   case only the value on the edge of an arbitrary 2 mm sphere.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'matlab'), fullfile(here,'Heatmodel','matlab'));
v1 = fileparts(here);
fprintf('\n================ verify_v2 ================\n');
pass = true;

C  = model_config('Pdel',1,'upper','water','bath_NaCl_pct',0.9,'T_bath_C',25);
zn = -(0:0.5e-3:30e-3).';
base = {'f',C.f,'d',C.d,'h',C.h,'L',C.L,'ke',C.ke,'Idist',C.Idist, ...
        'Nfft',512,'dx',C.dx,'Pdel',1,'zn',zn, ...
        'xlim',[-10e-3 30e-3],'ylim',[-10e-3 10e-3],'Zw',C.Zw,'dzw',C.dzw};
wet = [C.epsr_upper, C.epsr_gel, C.epsr_air];

% ---- 1. interface continuity --------------------------------------------
cfg = struct(base{:}, 'epsr',wet, 'debug_interface',true);
s = sar_volume(cfg);  D = s.dbg;
for nm = {'Ex','Ey'}
    b = D.([nm{1} '_bath']);  g = D.([nm{1} '_gel']);
    e = norm(b(:)-g(:))/norm(g(:));
    ok = e < 1e-12;  pass = pass && ok;
    fprintf('  1  tangential %s continuous at z=0 : rel err %8.2e  %s\n', nm{1}, e, tick(ok));
end

% ---- 2. air above must reproduce v1 -------------------------------------
cfg_air = struct(base{:}, 'epsr',[C.epsr_air, C.epsr_gel, C.epsr_air]);
s_air = sar_volume(cfg_air);
old = fullfile(v1,'Heatmodel','matlab');
have_v1 = exist(fullfile(old,'sar_volume.m'),'file') == 2;
if ~have_v1
    % The shared package ships v2 on its own.  Do NOT fall through to the
    % comparison here: addpath on a missing folder only WARNS, so sar_volume
    % would resolve to v2's own copy and this test would compare v2 against
    % itself and pass vacuously.  Say it was skipped instead.
    disp('  2  v1 is not beside this copy -- comparison skipped, not run');
    disp('     (at build time: Pgel, Pup, Pdn and the whole deposition');
    disp('      map identical to machine precision)');
else
addpath(old,'-begin');                                  % v1 solver wins
s_v1 = sar_volume(rmfield(cfg_air,{'Zw','dzw'}));
rmpath(old);
for f = {'Pgel','Pup','Pdn'}
    a = s_air.(f{1});  b = s_v1.(f{1});
    e = abs(a-b)/max(abs(b), s_v1.Pgel);
    ok = e < 1e-12;  pass = pass && ok;
    fprintf('  2  air above %-5s v2 %.7e  v1 %.7e  rel %7.1e  %s\n', f{1},a,b,e,tick(ok));
end
dq = max(abs(s_air.q(:)-s_v1.q(:)))/max(abs(s_v1.q(:)));
ok = isequal(size(s_air.q),size(s_v1.q)) && dq < 1e-12;  pass = pass && ok;
fprintf('  2  air above, deposition map identical      max rel %7.1e  %s\n', dq, tick(ok));
end
ok = s_air.Pwat == 0;  pass = pass && ok;
fprintf('  2  air above, bath absorbs exactly zero                     %s\n', tick(ok));

% ---- 3. the cut height must not matter ----------------------------------
% What leaves the top of the integrated column is counted as Pup instead of
% Pwat, so where the column is cut must not change what reaches the gel.
fprintf('  3  Zw independence:\n');
g0 = [];
for Zw = [20e-3 40e-3 80e-3]
    c = struct(base{:}, 'epsr',wet);  c.Zw = Zw;
    z = sar_volume(c);  T = z.Pgel+z.Pwat+z.Pup+z.Pdn;
    fr = z.Pgel/T;
    if isempty(g0), g0 = fr; end
    e = abs(fr-g0)/g0;  ok = e < 2e-3;  pass = pass && ok;
    fprintf('     Zw=%2.0f mm  gel %.4f  bath %.4f  escaping %.2e   drift %6.1e  %s\n', ...
            1e3*Zw, fr, z.Pwat/T, z.Pup/T, e, tick(ok));
end

% ---- 4. budget closes ----------------------------------------------------
T = s.Pgel+s.Pwat+s.Pup+s.Pdn;
fr = [s.Pgel s.Pwat s.Pup s.Pdn]/T;
ok = abs(sum(fr)-1) < 1e-12;  pass = pass && ok;
fprintf('  4  fractions sum to 1                                       %s\n', tick(ok));
fprintf('     saline bath: gel %.1f%%  bath %.1f%%  up %.2f%%  down %.2f%%\n', 100*fr);

% ---- 5. uniform medium, wire ON the plane -------------------------------
% All three layers the same and h = 0: the half-spaces are mirror images, so
% the power absorbed above must equal that below.
%   An earlier version of this test left the wire at h = 0.5 mm and expected
%   a ratio near 1.  That expectation was wrong, not the code -- at h > 0 the
%   wire's own near field sits inside the UPPER region, so that side
%   legitimately dominates (it came out 16x).  Putting the source on the
%   interface is what makes the two sides comparable.
cfg_h = struct(base{:}, 'epsr',[C.epsr_gel C.epsr_gel C.epsr_gel]);
cfg_h.h = 0;  cfg_h.Zw = 30e-3;  cfg_h.dzw = 0.5e-3;
sh = sar_volume(cfg_h);
r = sh.Pwat/sh.Pgel;  e = abs(r-1);  ok = e < 5e-3;  pass = pass && ok;
fprintf('  5  uniform medium, wire on z=0: above/below = %.6f          %s\n', r, tick(ok));

% ---- 6. the braid return arm --------------------------------------------
% 6a  L_sheath = 0 must leave the old answer untouched, or every earlier
%     result in the write-ups silently moves.
c_off = struct(base{:}, 'epsr',wet);
c_zero = c_off;  c_zero.L_sheath = 0;
s_off = sar_volume(c_off);  s_zero = sar_volume(c_zero);
dq = max(abs(s_off.q(:)-s_zero.q(:)));
ok = dq == 0;  pass = pass && ok;
fprintf('  6a L_sheath=0 vs absent: max |dq| = %8.2e (want exactly 0)   %s\n', ...
        dq, tick(ok));

% 6b  The point charge, measured directly.  |kx * Itilde(kx)| at large kx is
%     the size of the jump in I at the feed: 1.0 (in units of the peak element
%     current) if the current is truncated there, 0 if it is continuous.
c_on = c_off;  c_on.L_sheath = C.L_sheath;
tl = zeros(1,2);  k = 0;
for cc = {c_off, c_on}
    k = k + 1;
    c = cc{1};  c.debug_current = true;
    s6 = sar_volume(c);
    top = abs(s6.kv) > 0.9*max(abs(s6.kv));         % the top decade of kx
    tl(k) = mean(abs(s6.kv(top).*s6.It(top)));      % -> |jump| or -> 0
end
ok = abs(tl(1) - 1) < 0.05 && tl(2) < 0.02;  pass = pass && ok;
fprintf(['  6b jump in I at the feed, from the spectral tail:\n' ...
         '     truncated %.4f (the current itself)   braid arm %.4f      %s\n'], ...
        tl(1), tl(2), tick(ok));

% 6c  ... and the consequence: the peak temperature stops depending on the
%     grid.  This is the number the model is actually asked for.
fprintf('  6c peak temperature vs grid (50 W, 20 s, everything at 25 C):\n');
gp = {'T_gel',25,'T_water',25,'Pgen',50,'t_stim',20, ...
      'surf_bc','robin','cath_bc','cooled'};
Tt = [];  Ts = [];
for dh = [1e-3 0.5e-3 0.25e-3]
    Tt(end+1) = gui_model(gp{:},'dh',dh,'L_sheath_mm',0).info.T_end;  %#ok<AGROW>
    Ts(end+1) = gui_model(gp{:},'dh',dh,'L_sheath_mm',1e3*C.L_sheath).info.T_end; %#ok<AGROW>
    fprintf('     dh %4.2f mm   truncated %7.1f C   braid arm %7.1f C\n', ...
            1e3*dh, Tt(end), Ts(end));
end
st = (max(Tt)-min(Tt))/mean(Tt);   ss = (max(Ts)-min(Ts))/mean(Ts);
ok = ss < 0.02 && st > 0.10;  pass = pass && ok;
fprintf('     spread: truncated %.1f %%, braid arm %.2f %%                 %s\n', ...
        100*st, 100*ss, tick(ok));

% 6d  Nothing may be created or destroyed by adding the arm.
s_on = sar_volume(c_on);
T6 = s_on.Pgel + s_on.Pwat + s_on.Pup + s_on.Pdn;
e = abs(T6/s_on.cfg.Pdel - 1);
ok = e < 1e-12;  pass = pass && ok;
fprintf('  6d budget still closes with the arm on: rel %8.2e            %s\n', ...
        e, tick(ok));

% ---- 7. duty cycling ----------------------------------------------------
% The gate is applied inside the time loop, so the thing that can go wrong is
% a time step straddling a pulse edge and delivering the wrong energy.  These
% check the energy actually delivered, not the temperature it produced.
p = {'T_gel',37,'T_water',37,'Pgen',40,'t_stim',60,'dh',1e-3, ...
     'surf_bc','robin','cath_bc','cooled'};
Rcw = gui_model(p{:});
R1  = gui_model(p{:}, 'pulse_period',1, 'duty',1);
e = max(abs(R1.T_max - Rcw.T_max));
ok = e == 0;  pass = pass && ok;
fprintf('  7a duty = 1 reproduces continuous: max |dT| = %8.2e K        %s\n', ...
        e, tick(ok));

worst = -Inf;  emax = 0;
fprintf('  7b delivered energy, and the peak at a fixed 12 W mean:\n');
for Tp = [0.2 1 5 30]
    Rp = gui_model(p{:}, 'pulse_period',Tp, 'duty',0.3);
    want = 0.3*Rp.info.P_in*60;
    e = abs(Rp.audit.E_in(end) - want)/want;   emax = max(emax, e);
    fprintf('     period %5.1f s : E_in %9.4f J want %9.4f J (rel %7.1e), peak %6.2f C\n', ...
            Tp, Rp.audit.E_in(end), want, e, Rp.info.T_worst);
    ok = Rp.info.T_worst >= worst - 1e-9;   % monotone in the period
    pass = pass && ok;
    worst = Rp.info.T_worst;
end
ok = emax < 1e-10;  pass = pass && ok;
fprintf('     energy exact to %8.2e, and the peak rises with the period   %s\n', ...
        emax, tick(ok));

Rmean = gui_model(p{:}, 'Pgen',12);
ok = Rmean.info.T_worst <= worst + 1e-9;  pass = pass && ok;
fprintf(['  7c continuous at the same 12 W mean is the floor: %.2f C ' ...
         'vs %.2f C   %s\n'], Rmean.info.T_worst, worst, tick(ok));

fprintf('\n  %s\n\n', pick(pass,'ALL CHECKS PASSED','SOMETHING FAILED'));
end
function s = tick(ok),   if ok, s='ok'; else, s='FAIL'; end, end
function s = pick(c,a,b), if c, s=a; else, s=b; end, end
