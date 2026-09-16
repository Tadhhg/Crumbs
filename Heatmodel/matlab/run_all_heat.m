function run_all_heat()
%RUN_ALL_HEAT  Everything, in the right order, with the console captured to
%              results/full_heat_log.txt.
%
%   1. verify_heat_model      validation suite
%   2. run_heat_model         baseline: surface pinned at 37 degC   (~15 min)
%   3. run_heat_model robin   the physical convective surface       (~17 min)
%   4. run_heat_sensitivity   22-case parameter study               (~10 min)
%
%   Total about 45 minutes.

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(here); addpath(fullfile(fileparts(root),'matlab'));
if ~exist(fullfile(root,'results'),'dir'), mkdir(fullfile(root,'results')); end

logf = fullfile(root,'results','full_heat_log.txt');
if exist(logf,'file'), delete(logf); end
diary(logf); diary on;
c0 = tic;

verify_heat_model;

run_heat_model;

fprintf('\n\n########## the same model with a CONVECTIVE surface ##########\n');
run_heat_model('surf_bc','robin','cath_bc','cooled','tag','_robin');

run_heat_sensitivity;

% ---- exercise the deliverable functions on the saved solution -----------
clear T_gel P_for_T time_to_T
fprintf('\n\n########## the deliverable, T_gel ##########\n');
fprintf('  T_gel(P, t)                peak temperature in the block\n');
fprintf('  T_gel(P, t, x, z)          on the y = 0 plane\n');
fprintf('  T_gel(P, t, x, y, z)       anywhere in the block\n\n');
fprintf('   P_gen[W]   t[s]      T_peak[C]   T(6,0,-2mm)[C]   T(6,0,-10mm)[C]\n');
for P = [5 10 20 50]
    for tt = [1 10 60]
        fprintf('   %7.1f %6.1f  %12.2f %14.2f %16.2f\n', P, tt, ...
            T_gel(P,tt), T_gel(P,tt,6e-3,-2e-3), T_gel(P,tt,6e-3,-10e-3));
    end
end
fprintf('\n  steady state:\n');
for P = [5 10 20]
    fprintf('   %7.1f    Inf  %12.2f %14.2f %16.2f\n', P, ...
        T_gel(P,Inf), T_gel(P,Inf,6e-3,-2e-3), T_gel(P,Inf,6e-3,-10e-3));
end

fprintf('\n  off the y = 0 plane, at t = 60 s, P = 10 W:\n');
for yv = [0 1 2 5]*1e-3
    fprintf('   y = %4.1f mm : T(6mm, y, -2mm) = %.2f C\n', 1e3*yv, ...
            T_gel(10, 60, 6e-3, yv, -2e-3));
end

fprintf('\n  thermal dose and flags, 60 s at 12 W:\n');
tt = linspace(0,60,600);
[Tt,S] = T_gel(12, tt);
fprintf('   T_peak(60 s) = %.2f C, CEM43 = %.3e equivalent minutes\n', ...
        Tt(end), S.CEM43);
fprintf('   hot spot at (%.2f, %.2f, %.2f) mm\n', 1e3*S.location(end).x, ...
        1e3*S.location(end).y, 1e3*S.location(end).z);
for q = 1:numel(S.flags), fprintf('   flag: %s\n', S.flags{q}); end

fprintf('\n  inverse functions:\n');
for tq = [1 10 60 Inf]
    fprintf('   P_for_T(100 degC, %5.0f s) = %8.3f W generator\n', tq, ...
            P_for_T(100, tq));
end
for Pq = [10 20 50 100]
    fprintf('   time_to_T(%3.0f W, 100 degC) = %8.3f s\n', Pq, time_to_T(Pq,100));
end

fprintf('\n  pulsed drive, 30 %% duty at 1 Hz, 60 s, same average power:\n');
Pavg = P_for_T(100, 60);
Tcw  = T_gel(Pavg, 60);
Tp   = T_gel(@(tt) (Pavg/0.3)*(mod(tt,1) < 0.3), linspace(0,60,4000));
fprintf('   CW at %.2f W                 -> T_peak = %.2f C\n', Pavg, Tcw);
fprintf('   %.2f W at 30 %% duty (same mean) -> T_peak = %.2f C\n', Pavg/0.3, max(Tp));

fprintf('\n\nALL DONE in %.1f min\n', toc(c0)/60);
diary off;
end
