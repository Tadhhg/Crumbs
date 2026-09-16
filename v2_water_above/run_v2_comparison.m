function run_v2_comparison()
%RUN_V2_COMPARISON  The bench run that was photographed -- 50 W CW, 6 minutes,
%                   everything at 25 C -- through every version of the model.
%
%   v1     air above the gel, current truncated at the feed
%   v2     the antenna immersed in the bath, current still truncated
%   v2.1   ... and the feed current returning along the coax braid
%
%   THE SHAPE TEST AT THE END DOES NOT SEPARATE THE FEED MODELS, and saying so
%   is more useful than dressing it up.  It asks what temperature the model
%   would have to be calling the edge of its blob for that blob to be the size
%   the photograph shows, once from the length and once from the depth.  Two
%   different answers would mean the proportions are wrong.  Every version,
%   including v1, answers within 3 K -- so the photograph constrains the
%   MAGNITUDE and says nothing about the feed.  That is consistent with what v2
%   already found ("the shape barely moves"), and it means the case for the
%   braid arm has to rest on the things the photograph cannot see: that there
%   is no point charge (verify_v2 test 6b), that the peak stops depending on
%   the grid (6c), and that nothing has to be blanked out of the map.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'matlab'), fullfile(here,'Heatmodel','matlab'));

% Log-spaced output times, not the GUI's fixed increment.  This script asks
% "when did the peak cross 100 C", and on a 360 s run the GUI's 10 s steps can
% only ever answer "some time in the first 10 seconds".  GUI_MODEL takes the
% grid explicitly for exactly this.
b = {'T_gel',25,'T_water',25,'Pgen',50,'t_stim',360, ...
     'surf_bc','robin','cath_bc','cooled', ...
     't_out',unique([0 logspace(-3, log10(360), 240)])};
fprintf('\n=========== 50 W CW, 360 s, gel and bath at 25 C ===========\n\n');

R = struct();
R.v1   = gui_model(b{:}, 'upper','air',   'L_sheath_mm',0);
R.di   = gui_model(b{:}, 'upper','water', 'bath_NaCl_pct',0,   'L_sheath_mm',0);
R.sal  = gui_model(b{:}, 'upper','water', 'bath_NaCl_pct',0.9, 'L_sheath_mm',0);
R.di2  = gui_model(b{:}, 'upper','water', 'bath_NaCl_pct',0);
R.sal2 = gui_model(b{:}, 'upper','water', 'bath_NaCl_pct',0.9);
nm = {'v1    air above (the old one)', 'v2    deionised/tap bath', ...
      'v2    0.9% saline bath',        'v2.1  + braid arm, DI bath', ...
      'v2.1  + braid arm, saline'};
fn = {'v1','di','sal','di2','sal2'};
N  = numel(fn);

fprintf('%-30s %8s %8s %8s %9s %9s\n', ...
        '', 'gel %', 'bath %', 'Pgel W', 'peak C', 'probe C');
for k = 1:N
    i = R.(fn{k}).info;
    fprintf('%-30s %8.1f %8.1f %8.2f %9.1f %9.1f\n', nm{k}, ...
            100*i.frac_gel, 100*i.frac_bath, i.Pgel, i.T_end, i.T_probe_end);
end

fprintf('\n  bath permittivity used:\n');
for k = 2:3
    e = R.(fn{k}).info.epsr_upper;
    fprintf('    %-26s eps'' %6.2f   sigma %.3f S/m\n', nm{k}, ...
            real(e), -imag(e)*2*pi*2.45e9*8.8541878128e-12);
end

fprintf('\n  time for the peak to cross a threshold:\n');
fprintf('    %-28s %8s %8s %8s\n','','43 C','60 C','100 C');
for k = 1:N
    r = R.(fn{k});  s = sprintf('    %-28s', nm{k});
    for lev = [43 60 100]
        j = find(r.T_max >= lev, 1);
        if isempty(j), s = [s sprintf(' %8s','never')]; %#ok<AGROW>
        else,          s = [s sprintf(' %8.2f', r.t(j))]; end %#ok<AGROW>
    end
    fprintf('%s\n', s);
end

% ---- is the peak a prediction at all? -----------------------------------
fprintf('\n  is the reported peak a number, or an artefact?\n');
for k = 1:N
    i = R.(fn{k}).info;
    if i.r_feed_excl > 0
        why = sprintf(['blanked over %.0f mm; the peak is the value on the ' ...
                       'edge of that sphere'], 1e3*i.r_feed_excl);
    else
        why = 'nothing blanked; the peak is where the model puts it';
    end
    fprintf('    %-28s %s\n', nm{k}, why);
end

fprintf('\n  ---- against the photographs ----\n');
fprintf('  observed: colour first visible between t=0 and t=29 s; gel intact,\n');
fprintf('            no boiling, no bubbles; coloured core ~9.5 x 6.3 mm at 332 s.\n\n');
for k = 1:N
    i = R.(fn{k}).info;
    fprintf('  %-28s predicts peak %6.1f C  -> %s\n', nm{k}, i.T_end, ...
            verdict(i.T_end));
end

% ---- the shape test ------------------------------------------------------
% Peak temperature depends on eta, which nobody measured, so it cannot falsify
% anything on its own.  The SIZE of the coloured region depends on the shape of
% the deposition instead, so ask about that -- both ways round.
obsL = 9.5;  obsD = 6.3;  tobs = 332;
fprintf(['\n  the coloured core is %.1f mm long and %.1f mm deep at t = %g s.\n' ...
         '  (a) what would the model have to be calling the edge of its blob?\n\n'], ...
        obsL, obsD, tobs);
fprintf('    %-28s %10s %10s %10s\n', '', 'from len', 'from depth', 'disagree');
for k = 1:N
    [tL, tD] = shape_levels(R.(fn{k}), tobs, obsL, obsD);
    if isnan(tL) || isnan(tD), dis = NaN; else, dis = abs(tL-tD); end
    fprintf('    %-28s %8.1f C %8.1f C %8.1f K\n', nm{k}, tL, tD, dis);
end
fprintf(['\n    "disagree" would be large if the proportions were wrong.  It is ' ...
         'small\n    for EVERY version, v1 included, so the photograph cannot ' ...
         'tell the feed\n    models apart on shape and never could.  What it ' ...
         'does say is that\n    the proportions were right all along and only ' ...
         'the magnitude was not.\n']);

fprintf(['\n  (b) the other way round -- how big is the model''s blob at a ' ...
         'temperature\n      an indicator might actually trip at?  (length x ' ...
         'depth, mm)\n\n']);
fprintf('    %-28s %14s %14s %14s\n', '', '40 C', '60 C', '100 C');
for k = 1:N
    r = R.(fn{k});
    [~, it] = min(abs(r.t - tobs));
    A = r.T_plane(:,:,it);  A(~r.pmask_y0) = -Inf;
    s = sprintf('    %-28s', nm{k});
    for lev = [40 60 100]
        [Lk, Dk] = blob(A, r.x*1e3, -r.z*1e3, lev);
        s = [s sprintf(' %6.1f x%5.1f ', Lk, Dk)]; %#ok<AGROW>
    end
    fprintf('%s\n', s);
end
fprintf(['\n    Every version puts a 40 C region several times the size of the ' ...
         'coloured\n    core.  That is the same over-prediction the peak shows, ' ...
         'with the same\n    unmeasured eta behind it.  MEASURE THE FORWARD AND ' ...
         'REFLECTED POWER.\n']);

% how much power the bench behaviour implies actually reached the gel
d = R.sal2.info.T_probe_end - 25;
fprintf('\n  v2.1 saline: %.2f K at the probe per generator watt.\n', d/50);
fprintf('  A measured probe rise of 10-25 K would imply %.1f-%.1f W into the gel.\n', ...
        10/(d/50), 25/(d/50));
fprintf('\nDONE\n');
end

% =========================================================================
function [tL, tD] = shape_levels(r, tobs, obsL, obsD)
%SHAPE_LEVELS  The contour levels whose blob is obsL long and obsD deep.
[~, it] = min(abs(r.t - tobs));
A  = r.T_plane(:,:,it);   A(~r.pmask_y0) = -Inf;
xm = r.x*1e3;  zm = -r.z*1e3;
lv = linspace(min(r.opt.T_water, r.opt.T_gel) + 0.5, max(A(:)) - 0.1, 400);
Ls = nan(size(lv));  Ds = nan(size(lv));
for j = 1:numel(lv)
    [Ls(j), Ds(j)] = blob(A, xm, zm, lv(j));
end
tL = cross_at(lv, Ls, obsL);
tD = cross_at(lv, Ds, obsD);
end

% =========================================================================
function [len, dep] = blob(A, xm, zm, level)
%BLOB  Bounding box of {T >= level} on the y = 0 plane, in mm.
m = A >= level;
if ~any(m(:)), len = 0; dep = 0; return; end
ix = any(m,2);  iz = any(m,1);
len = xm(find(ix,1,'last')) - xm(find(ix,1));
dep = zm(find(iz,1,'last')) - zm(find(iz,1));
end

% =========================================================================
function t = cross_at(lv, sz, target)
%CROSS_AT  The level at which the blob dimension SZ falls through TARGET.
%   sz decreases with the level, so take the first crossing from above.
j = find(sz <= target, 1);
if isempty(j) || j == 1, t = NaN; return; end
t = interp1(sz([j-1 j]), lv([j-1 j]), target);
end

% =========================================================================
function s = verdict(T)
if     T >= 100, s = 'boils -- ruled out by the photos';
elseif T >= 60,  s = 'still far hotter than the gel looks';
else,            s = 'physically plausible';
end
end
