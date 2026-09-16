function SS = run_heat_sensitivity(varargin)
%RUN_HEAT_SENSITIVITY  Parameter study for the thermal model.
%
%   Every case is run on the SAME 0.5 mm grid to 120 s, plus its exact steady
%   state, so the cases are comparable with each other.  The grid error of
%   0.5 mm relative to the 0.25 mm baseline is quantified separately in
%   VERIFY_HEAT_MODEL, test 6, and is small compared with the spreads below.
%
%   Reported for each case:
%     Theta at 1 s, 10 s, 120 s and at steady state, in K per the reference
%     150 W delivered, at the hot spot (feed sphere excluded); the hot-spot
%     location; and P100, the generator power at which the peak just reaches
%     100 degC after 60 s -- the number a bench operator actually needs.
%
%   Writes results/heat_sensitivity.mat and results/heat_sensitivity_log.txt,
%   and figures/fig9_sensitivity.png.

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(here); addpath(fullfile(fileparts(root),'matlab'));

dh = 0.5e-3;  t_end = 120;
base = {'dh',dh,'t_end',t_end,'Nt',81};
if nargin > 0, base = [base varargin]; end

cases = {
 'A  baseline (surface pinned at 37 C)',            {}
 'B  Robin surface, h = 1460 (plate laminar)',      {'surf_bc','robin','h_surf',1460}
 'B  Robin surface, h =  673 (duct Dittus-Boelter)',{'surf_bc','robin','h_surf',673}
 'B  Robin surface, h =  300',                      {'surf_bc','robin','h_surf',300}
 'B  Robin surface, h =  103 (laminar duct, worst)',{'surf_bc','robin','h_surf',103}
 'B  Robin surface, h = 5000 (near-Dirichlet)',     {'surf_bc','robin','h_surf',5000}
 'C  catheter cooled, h_cath = 857 (with h = 1460)',{'surf_bc','robin','h_surf',1460,'cath_bc','cooled','h_cath',857}
 'C  catheter INSULATING (with h = 1460)',          {'surf_bc','robin','h_surf',1460,'cath_bc','insulating'}
 'C  catheter INSULATING (pinned surface)',         {'cath_bc','insulating'}
 'D  floor held at 37 C (water-cooled dish)',       {'bot_bc','dirichlet'}
 'E  ROI side walls adiabatic',                     {'side_bc','adiabatic'}
 'F  k_gel = 0.50 W/(m K)',                         {'k_gel',0.50}
 'F  k_gel = 0.55 W/(m K)',                         {'k_gel',0.55}
 'F  k_gel = 0.628 W/(m K) (pure water at 37 C)',   {'k_gel',0.628}
 'G  c_gel = 4000 J/(kg K)',                        {'c_gel',4000}
 'G  muscle thermal properties',                    {'k_gel',0.49,'c_gel',3421,'rho_gel',1090}
 'H  gel permittivity taken at 37 C, not 25 C',     {'T_medium',37}
 'I  feed exclusion radius 1 mm',                   {'r_feed_excl',1e-3}
 'I  feed exclusion radius 3 mm',                   {'r_feed_excl',3e-3}
 'I  feed source regularised (capped on 2 mm shell)',{'feed_fix','cap'}
 'J  ROI 20 mm wide',                               {'Wroi',20e-3}
 'J  ROI 40 mm wide',                               {'Wroi',40e-3}
};

nC = size(cases,1);
SS = struct('name',cases(:,1),'args',cases(:,2));
tout = unique([0 logspace(-2, log10(t_end), 60) 1 10 60 t_end]);
i1  = find(tout==1);  i10 = find(tout==10);
i60 = find(tout==60); iE  = numel(tout);

fid = fopen(fullfile(root,'results','heat_sensitivity_log.txt'),'w');
hdr = sprintf(['%-48s %8s %8s %8s %9s   %-18s %8s\n'], 'case', ...
       'Th(1s)','Th(10s)','Th(120s)','Th_ss','hot spot (x,y,d) mm','P100@60s');
line = [repmat('-',1,120) newline];
emit(fid, [newline hdr line]);

for c = 1:nC
    Tc = thermal_config(base{:}, cases{c,2}{:});
    [xg,yg,zg] = roi_nodes(Tc);
    q  = em_source(Tc, xg, yg, zg);
    G  = heat_grid(Tc, xg, yg, zg);
    [q, finfo] = feed_regularise(q, xg, yg, zg, Tc);
    pmask = finfo.mask_peak;
    S = heat_solve(q, G, Tc, struct('t_out',tout,'verbose',false, ...
                                    'Nt3',2,'quiet',true,'peakmask',pmask));
    ths = heat_steady(q, G, 1e-9);
    thm = ths;  thm(~pmask) = -Inf;
    [thss, li] = max(thm(:));
    [a,b,cc] = ind2sub(size(ths), li);
    P100 = Tc.Pdel_ref*(100-Tc.T_blood)/(Tc.eta*S.Th_max(i60));

    SS(c).Th   = S.Th_max(2:end);
    SS(c).t    = S.t(2:end);
    SS(c).Th_ss= thss;
    SS(c).loc  = [xg(a) yg(b) zg(cc)];
    SS(c).P100 = P100;
    SS(c).P_in = S.P_in;
    SS(c).Th_ss_y0 = squeeze(ths(:,S.iy0,:));
    SS(c).x = xg;  SS(c).z = zg;

    emit(fid, sprintf('%-48s %8.1f %8.1f %8.1f %9.1f   (%5.1f,%4.1f,%5.2f) %8.3f\n', ...
        cases{c,1}, S.Th_max(i1), S.Th_max(i10), S.Th_max(iE), thss, ...
        1e3*xg(a), 1e3*yg(b), -1e3*zg(cc), P100));
end

% ---- the nonlinear check: temperature-dependent conductivity -------------
emit(fid, [newline '--- K  temperature-dependent gel properties (breaks linearity) ---' newline]);
emit(fid, sprintf(['  The Stogryn model for 0.9 %% saline at 2.45 GHz gives\n' ...
    '  sigma_eff = 2.782 S/m at 37 C rising to 3.18 at 60 C (+0.62 %%/K) and\n' ...
    '  eps'''' falling at -0.38 %%/K.  Re-running with q scaled pointwise by\n' ...
    '  sigma(T)/sigma(37) (|E|^2 held fixed) bounds the error in the exact\n' ...
    '  power scaling that T_gel relies on.\n']));
Tc = thermal_config(base{:});
[xg,yg,zg] = roi_nodes(Tc);
q0 = em_source(Tc, xg, yg, zg);
G  = heat_grid(Tc, xg, yg, zg);
[q0, finfo] = feed_regularise(q0, xg, yg, zg, Tc);
pmask = finfo.mask_peak;
for Pg = [5 10 20]
    s  = Pg*Tc.eta/Tc.Pdel_ref;
    lin = heat_solve(s*q0, G, Tc, struct('t_out',[0 60],'verbose',false, ...
                                'Nt3',2,'quiet',true,'peakmask',pmask)); %#ok<*AGROW>
    nl  = nonlinear_solve(s*q0, G, Tc, 60, pmask);
    emit(fid, sprintf(['  P_gen = %5.1f W : T_max(60 s) linear %7.2f C, ' ...
        'temperature-dependent %7.2f C  (%+.2f %%)\n'], Pg, ...
        Tc.T_blood+lin.Th_max(2), Tc.T_blood+nl, ...
        100*((nl)/(lin.Th_max(2))-1)));
end

emit(fid, [newline]);
fclose(fid);
save(fullfile(root,'results','heat_sensitivity.mat'),'SS');
make_sensitivity_figure(SS, root);
fprintf('\nwrote results/heat_sensitivity_log.txt and figures/fig9_sensitivity.png\n');
end

% =========================================================================
function emit(fid, s)
fprintf('%s', s);  fprintf(fid, '%s', s);
end

function [xg,yg,zg] = roi_nodes(Tc)
dh = Tc.dh;
xg = (Tc.xlim(1):dh:Tc.xlim(2)+1e-12).';
yg = (Tc.ylim(1):dh:Tc.ylim(2)+1e-12).';
zg = -(0:dh:Tc.d+1e-12).';
end

function q = em_source(Tc, xg, yg, zg)
persistent cache
key = matlab.lang.makeValidName(sprintf('k%g_%g_%g_%g', Tc.dh, Tc.Wroi, ...
                                        Tc.T_medium, Tc.Pdel_ref));
if ~isempty(cache) && isfield(cache,key), q = cache.(key); return; end
C = model_config('Pdel',Tc.Pdel_ref,'T_C',Tc.T_medium);
cfg = struct('f',C.f,'epsr',C.epsr,'d',C.d,'h',C.h,'L',C.L,'ke',C.ke, ...
             'Idist',C.Idist,'Nfft',C.Nfft,'dx',C.dx,'Pdel',C.Pdel, ...
             'zn',zg,'xlim',Tc.xlim,'ylim',Tc.ylim);
o = sar_volume(cfg);
st = round(Tc.dh/C.dx);
q = o.q(1:st:end, 1:st:end, :);
assert(isequal(size(q),[numel(xg) numel(yg) numel(zg)]),'grid mismatch');
cache.(key) = q;
end

function thmax = nonlinear_solve(q, G, Tc, tend, pmask)
%NONLINEAR_SOLVE  Same explicit scheme, but the deposition follows sigma(T).
%   |E|^2 is held fixed, so this captures the change in absorbed power but not
%   the change in field distribution; it is a bound, not a coupled solution.
eps0 = 8.8541878128e-12;  f = 2.45e9;  w = 2*pi*f;
N = 0.9*10/58.443;
Tt = (20:5:120).';
sg = arrayfun(@(TT) w*eps0*(-imag(saline_stogryn(f,TT,N))), Tt);
sg = sg/interp1(Tt, sg, Tc.T_blood);            % normalise to 37 degC
th = zeros(G.nx,G.ny,G.nz);
dt0 = Tc.cfl*G.dt_stable;  t = 0;  rc = G.rho*G.c;
while t < tend - 1e-12
    dt = min(dt0, tend-t);
    fac = interp1(Tt, sg, min(max(Tc.T_blood+th,20),120), 'linear');
    th = th + dt*(G.alpha*heat_lap(th,G) + fac.*q/rc);
    th(G.pin) = 0;
    t = t + dt;
end
thm = th;  thm(~pmask) = -Inf;
thmax = max(thm(:));
end
