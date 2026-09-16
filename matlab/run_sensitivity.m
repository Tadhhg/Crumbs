function T = run_sensitivity()
%RUN_SENSITIVITY  Parameter study for the 2.45 GHz Phytagel field model.
%
%   Quantifies how the modelled field depends on each assumption that is not
%   fixed by the problem statement:
%     A  the current-wave hypothesis on the monopole  (largest uncertainty)
%     B  the standoff height h of the wire above the gel
%     C  the gel complex permittivity (salt, temperature, formulation)
%     D  the gel depth
%     E  what backs the gel (air vs a metal dish)
%
%   Metrics reported for each case:
%     Rin        implied input resistance at the feed [ohm]
%     %gel       share of delivered power absorbed in the gel
%     E2, E10, E25   |E|_rms at x = 6 mm and 2 / 10 / 25 mm depth [V/m]
%     z1e        depth at which |E| falls to 1/e of its near-surface value

here = fileparts(mfilename('fullpath'));  root = fileparts(here);
if ~exist(fullfile(root,'results'),'dir'), mkdir(fullfile(root,'results')); end

rows = {};
% ---------------- A. current-wave hypothesis -----------------------------
rows{end+1} = {'A baseline: sinusoid, L = lambda_e/4', {'ke_case','quarterwave'}};
rows{end+1} = {'A sinusoid, bare wire on interface',   {'ke_case','bare_iface'}};
rows{end+1} = {'A sinusoid, wire sees air only',       {'ke_case','bare_air'}};
rows{end+1} = {'A uniform current',                    {'Idist','uniform'}};
rows{end+1} = {'A triangular current',                 {'Idist','triangular'}};
% ---------------- B. standoff --------------------------------------------
for hh = [0.1 0.25 1 2 5]
    rows{end+1} = {sprintf('B standoff h = %.2f mm', hh), {'h',hh*1e-3}}; %#ok<AGROW>
end
% ---------------- C. medium ----------------------------------------------
rows{end+1} = {'C 1.1 pct NaCl, 25 C',            {'NaCl_pct',1.1}};
rows{end+1} = {'C 0.9 pct NaCl, 20 C',            {'T_C',20}};
rows{end+1} = {'C 0.9 pct NaCl, 37 C',            {'T_C',37}};
rows{end+1} = {'C muscle-equivalent target',     {'epsr_override', muscle_epsr()}};
rows{end+1} = {'C pure water (no salt), 25 C',   {'NaCl_pct',1e-6}};
% ---------------- D. gel depth -------------------------------------------
rows{end+1} = {'D gel depth 20 mm',              {'d',20e-3}};
rows{end+1} = {'D gel depth 40 mm',              {'d',40e-3}};
% ---------------- E. backing ---------------------------------------------
rows{end+1} = {'E metal (PEC) dish under gel',   {'backing','pec'}};

fprintf('\n%-38s %7s %6s %9s %9s %9s %8s\n', 'case','Rin[ohm]','%gel', ...
        'E@2mm','E@10mm','E@25mm','z1e[mm]');
fprintf('%s\n', repmat('-',1,92));

T = struct('name',{},'Rin',{},'pgel',{},'E2',{},'E10',{},'E25',{},'z1e',{});
for q = 1:numel(rows)
    [m, nm] = run_case(rows{q}{1}, rows{q}{2});
    fprintf('%-38s %7.1f %6.1f %9.1f %9.1f %9.1f %8.2f\n', ...
        nm, m.Rin, m.pgel, m.E2, m.E10, m.E25, m.z1e);
    T(end+1) = m; %#ok<AGROW>
end
save(fullfile(root,'results','sensitivity.mat'),'T');
fprintf('\nsaved results/sensitivity.mat\n');
end

% =========================================================================
function [m, nm] = run_case(nm, ov)
if strcmp(ov{1},'epsr_override')
    C = model_config();
    C.epsr_gel = ov{2};
    C.epsr = [C.epsr_air, C.epsr_gel, C.epsr_bot];
else
    C = model_config(ov{:});
end
C.z = -(C.dz/2 : C.dz : C.d - C.dz/2);
cfg = struct('f',C.f,'epsr',C.epsr,'d',C.d,'h',C.h,'L',C.L,'ke',C.ke, ...
             'Idist',C.Idist,'Nfft',C.Nfft,'dx',C.dx,'z',C.z,'Pdel',C.Pdel);
o = layered_field(cfg);

E2m = abs(o.Ex).^2 + abs(o.Ey).^2 + abs(o.Ez).^2;
Erms = sqrt(E2m/2);
xg = o.x(:).'; zg = o.z(:);

switch lower(C.Idist)
    case 'sinusoid', I0 = sin(C.ke*C.L);
    otherwise,       I0 = 1;
end
if abs(I0) < 1e-6, I0 = 1e-6; end

[~,ix] = min(abs(xg - 6e-3));
gv = @(zq) Erms(argmin(abs(zg + zq)), ix);
m.name = nm;
m.Rin  = 2*o.Ptot/abs(I0)^2;
m.pgel = 100*o.Pgel/o.Ptot;
m.E2   = gv(2e-3);
m.E10  = gv(10e-3);
m.E25  = gv(25e-3);

prof = Erms(:,ix);
i = find(prof <= prof(1)/exp(1), 1);
if isempty(i), m.z1e = NaN;
elseif i == 1, m.z1e = 0;
else
    w = (prof(i-1)-prof(1)/exp(1))/(prof(i-1)-prof(i));
    m.z1e = -1e3*(zg(i-1) + w*(zg(i)-zg(i-1)));
end
end

function i = argmin(v)
[~,i] = min(v);
end

function e = muscle_epsr()
%MUSCLE_EPSR  Skeletal muscle at 2.45 GHz from the four-pole Cole-Cole model
%   of Gabriel et al. (1996); see TISSUE_COLE_COLE for the full source chain.
%   Returns eps' = 52.73, sigma_eff = 1.738 S/m.
%   Used here only as the muscle-equivalence TARGET that the 30 % glycerol /
%   70 % saline Phytagel formulation of Qian et al. was designed to hit.
e = tissue_cole_cole(2.45e9, 'muscle');
end
