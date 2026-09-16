function extra_diagnostics()
%EXTRA_DIAGNOSTICS  Numbers quoted in the write-up that the main run does not
%   print: lesion dimensions, isotherm depths, per-watt sensitivities and the
%   Dirichlet/convective comparison.  Appends to results/full_heat_log.txt.

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(here); addpath(fullfile(fileparts(root),'matlab'));
diary(fullfile(root,'results','full_heat_log.txt')); diary on;

fprintf('\n\n########## extra diagnostics ##########\n');

M = load(fullfile(root,'results','heat_model.mat'));
R = load(fullfile(root,'results','heat_model_robin.mat'));
xm = M.x(:).'*1e3;  sm = -M.z(:).'*1e3;  t = M.t;
Tb = M.T_blood;  eta = M.eta;  Pref = M.Pdel_ref;

fprintf('\n--- temperature rise per watt of GENERATOR power ---\n');
for tq = [1 10 60 120 Inf]
    if isinf(tq), th = M.Th_max_ss; lbl = 'steady';
    else, [~,i1] = min(abs(t-tq)); th = M.Th_max(i1); lbl = sprintf('%g s',tq); end
    fprintf('  %-8s : %7.3f K/W  (peak)\n', lbl, th*eta/Pref);
end

fprintf('\n--- time to reach a temperature, at the field model reference 200 W ---\n');
for Tl = [43 50 60 100]
    fprintf('  %3d degC : %8.3f s\n', Tl, time_to_T(200, Tl));
end

fprintf('\n--- Dirichlet (as specified) vs the convective surface it stands for ---\n');
fprintf('   t [s]   Theta_max Dirichlet   Theta_max Robin(h=1460)   ratio   P100 D / P100 R\n');
for tq = [1 10 60 120 600]
    [~,i1] = min(abs(t-tq));
    a = M.Th_max(i1);  b = R.Th_max(i1);
    fprintf('  %6.0f   %17.1f   %22.1f   %5.2f   %6.2f / %.2f W\n', t(i1), a, b, b/a, ...
        Pref*63/(eta*a), Pref*63/(eta*b));
end
fprintf('  steady   %17.1f   %22.1f   %5.2f   %6.2f / %.2f W\n', ...
    M.Th_max_ss, R.Th_max_ss, R.Th_max_ss/M.Th_max_ss, ...
    Pref*63/(eta*M.Th_max_ss), Pref*63/(eta*R.Th_max_ss));
[~,ix6d] = min(abs(M.x - 6e-3));
fprintf('  surface rise at x = 6 mm, steady state: Dirichlet 0 K by construction,\n');
fprintf('  convective %.1f K per 150 W delivered = %.2f K per generator watt.\n', ...
        double(R.Th_ss_y0(ix6d,1)), double(R.Th_ss_y0(ix6d,1))*eta/Pref);
fprintf('  That film resistance is exactly what the Dirichlet condition throws away.\n');

fprintf('\n--- lesion size, 60 s exposure ---\n');
fprintf('  The y = 0 plane is masked to exclude the 2 mm feed sphere, so these\n');
fprintf('  dimensions describe the zone the RADIATING ELEMENT produces.\n');
[~,i60] = min(abs(t-60));
P100 = Pref*(100-Tb)/(eta*M.Th_max(i60));
kk = t <= t(i60);
pm = M.pmask_y0;
for Pg = [0.5 0.75 1.0]*P100
    s = Pg*eta/Pref;
    Tser = Tb + s*double(M.Th_y0(:,:,kk));
    D = cem43(t(kk), Tser);
    D(~pm) = 0;
    Tend = Tb + s*double(M.Th_y0(:,:,i60));
    Tend(~pm) = 0;
    fprintf('  P_gen = %5.2f W (T_peak = %6.1f degC):\n', Pg, max(Tend(:)));
    for lev = [240 60]
        m = D >= lev;
        if any(m(:))
            [ii,jj] = find(m);
            fprintf('     CEM43 >= %3d min : depth %.2f mm, x from %.1f to %.1f mm (width %.1f mm)\n', ...
                lev, max(sm(jj)), min(xm(ii)), max(xm(ii)), max(xm(ii))-min(xm(ii)));
        else
            fprintf('     CEM43 >= %3d min : none\n', lev);
        end
    end
    for lev = [43 50 60]
        m = Tend >= lev;
        if any(m(:))
            [ii,jj] = find(m);
            fprintf('     T >= %3d degC     : depth %.2f mm, x from %.1f to %.1f mm\n', ...
                lev, max(sm(jj)), min(xm(ii)), max(xm(ii)));
        end
    end
end

fprintf('\n--- how fast the peak equilibrates ---\n');
for tq = [1 10 30 60 120 300 600]
    [~,i1] = min(abs(t-tq));
    fprintf('  t = %5.0f s : %.1f %% of the steady-state rise\n', t(i1), ...
            100*M.Th_max(i1)/M.Th_max_ss);
end

fprintf('\n--- pulsed drive at a fixed mean power (60 s, mean = P100@60s) ---\n');
Pavg = P100;
fprintf('  CW at %.2f W                        -> T_peak = %6.2f degC\n', ...
        Pavg, T_gel(Pavg, 60));
for duty = [0.5 0.3 0.1]
    for fpulse = [1 0.1 0.02]
        Tp = T_gel(@(tt) (Pavg/duty)*(mod(tt,1/fpulse) < duty/fpulse), ...
                   linspace(0,60,6000));
        fprintf('  %5.1f %% duty at %5.3f Hz (%6.2f W peak) -> T_peak = %6.2f degC\n', ...
                100*duty, fpulse, Pavg/duty, max(Tp));
    end
end

fprintf('\n--- depth of the 43 and 50 degC isotherms under the element, 60 s ---\n');
[~,ix6] = min(abs(M.x - 6e-3));
for Pg = [5 10 12.42 20 50]
    s = Pg*eta/Pref;
    prof = Tb + s*double(squeeze(M.Th_y0(ix6,:,i60)));
    d43 = last_depth(sm, prof, 43);  d50 = last_depth(sm, prof, 50);
    fprintf('  P_gen = %5.2f W : T(surface..) peak %6.1f degC, 43 degC to %5.2f mm, 50 degC to %5.2f mm\n', ...
            Pg, max(prof), d43, d50);
end

diary off;
end

function d = last_depth(sm, prof, lev)
i = find(prof >= lev, 1, 'last');
if isempty(i), d = 0; return; end
if i == numel(prof), d = sm(end); return; end
w = (prof(i)-lev)/(prof(i)-prof(i+1));
d = sm(i) + w*(sm(i+1)-sm(i));
end
