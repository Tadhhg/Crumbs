function heat_gui()
%HEAT_GUI  Move a slider, watch the SAR and temperature maps change.
%
%   heat_gui
%
%   Pick a starting point from the top drop-down, then change anything.  The
%   maps redraw by themselves -- about a second for a change of power or time,
%   about five if the gel itself changed.  The line above the maps says what
%   the answer is in words.
%
%   WHAT YOU CAN CHANGE
%     gel        permittivity, conductivity, temperature
%     coolant    temperature, flow rate, and how the surface is cooled
%     catheter   how much of it lies on the gel, and what it does thermally
%     drive      generator power, exposure time, and pulsing (period and duty)
%     feed       how much of the coax braid carries the return current
%
%   eps' and sigma follow the gel temperature through the saline model until
%   you type a value into either box, which releases the tie.  BOTH BOXES ARE
%   ALWAYS EDITABLE -- the tick box only decides whether the pair is recomputed
%   when the gel temperature moves, never whether you may type in it.
%
%   IT OPENS ON 50 W CONTINUOUS WITH ROOM-TEMPERATURE COOLANT
%   Gel and water both start at 25 degC (the bench temperature the field model
%   is characterised at), the drive is continuous, and the exposure is 60 s.
%   That is far past boiling and the sentence above the maps says so; the point
%   of the preset is the shape of the maps and the time to a limit, not the
%   60 s number.  "bench baseline" in the drop-down restores 12 W with 37 degC
%   blood, and "the photographed run" restores 50 W for 6 minutes at 25 degC,
%   which is the case the bath photographs show.
%
%   PULSING.  Set "drive" to pulsed and the generator power becomes the power
%   DURING a pulse; the mean is duty x that, and the line above the maps gives
%   both.  Compare pulsed against continuous at equal MEAN power, not at equal
%   generator power -- at a fixed mean, a slower pulse always reaches a higher
%   peak, and that is the whole reason to look at it.
%
%   THE TIME BAR moves in fixed increments of seconds.  The frames it steps
%   through are 0, dt, 2 dt, ... with dt a round number chosen from the
%   exposure time, so half way along the bar is half way through the run.
%
%   TWO DEFAULTS THAT DIFFER FROM THE WRITE-UP BASELINE
%   This opens with realistic flow cooling and a catheter that is cooled by
%   that flow, not the write-up's perfect heat sink and ignored catheter.  It
%   has to: under the write-up baseline the flow rate and the catheter length
%   do nothing at all (see GUI_MODEL).  Both baselines are in the drop-downs,
%   and "write-up baseline" is one of the starting points.
%
%   THE FEED, AND WHERE THE GREY WEDGE WENT
%   Earlier versions blanked a 2 mm sphere at the feed and drew it as a grey
%   semicircle, because the prescribed antenna current was truncated there and
%   the deposition was singular.  It was an artefact of the source model, and
%   the bath photographs show it is not what the gel does: the heated zone is a
%   smooth blob a couple of millimetres along and a couple of millimetres down,
%   reaching back past the junction, with no hole in it.
%
%   The current now RETURNS along the outside of the coax braid instead of
%   stopping, which is what it does on a bench with no ground plane and no
%   balun -- the thing is an asymmetric dipole, not a monopole.  "Feed / braid
%   arm" is how much of the braid carries it, 12 mm by default, which makes the
%   arm the mirror of the element.  There is then no point charge, nothing is
%   blanked, and the peak the window reports is a real prediction.  Set the arm
%   to 0 for the old truncated current, grey wedge and all.
%
%   The probe point (6 mm along, 2 mm deep) is still quoted beside the peak,
%   and it is still the one to compare against a thermocouple -- it is a place
%   a thermocouple can actually go.
%
%   See also GUI_MODEL, VERIFY_GUI_MODEL, HEAT_SCENARIO, T_GEL.

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(here);      % the rest of the heating model, beside this file
hm_addpath();       % and the field model, wherever it has been put

R  = [];  em = [];
busy     = false;    % a solve is running
pending  = false;    % something changed while it was running
starting = true;     % suppress solves while the window is still being built
deg      = char(176);

% the drop-down wording is plain language; these are the codes GUI_MODEL wants
SURF = {'robin','dirichlet','adiabatic'};
UPPER = {'water','air'};
CATH = {'cooled','transparent','insulating'};
DHV  = [0.5e-3 1e-3 0.25e-3];

% starting points: name, then the fields it sets
P = @(varargin) struct(varargin{:});
% drv 1 = continuous, 2 = pulsed;  per = pulse period [s];  dut = duty [%];
% lf = braid arm length [mm], 0 = the old truncated current.
PRESET = { ...
 'room-temp water  --  50 W continuous, 60 s', ...
     P('pow',50,'tim',60,'tgel',25,'tw',25,'q',5,'lc',32,'surf',1,'cath',1, ...
       'drv',1,'per',1,'dut',30,'lf',12)
 'the photographed run  --  50 W for 6 min at 25 C', ...
     P('pow',50,'tim',360,'tgel',25,'tw',25,'q',5,'lc',32,'surf',1,'cath',1, ...
       'drv',1,'per',1,'dut',30,'lf',12)
 'bench baseline  --  12 W for 60 s', ...
     P('pow',12,'tim',60,'tgel',37,'tw',37,'q',5,'lc',32,'surf',1,'cath',1, ...
       'drv',1,'per',1,'dut',30,'lf',12)
 'gentle test  --  5 W for 60 s, stays testable', ...
     P('pow',5, 'tim',60,'tgel',37,'tw',37,'q',5,'lc',32,'surf',1,'cath',1, ...
       'drv',1,'per',1,'dut',30,'lf',12)
 'high power  --  40 W for 10 s', ...
     P('pow',40,'tim',10,'tgel',37,'tw',37,'q',5,'lc',32,'surf',1,'cath',1, ...
       'drv',1,'per',1,'dut',30,'lf',12)
 'pulsed  --  40 W at 30 % duty, 1 s period, 60 s', ...
     P('pow',40,'tim',60,'tgel',37,'tw',37,'q',5,'lc',32,'surf',1,'cath',1, ...
       'drv',2,'per',1,'dut',30,'lf',12)
 'slow pulsing  --  50 W at 50 % duty, 30 s period', ...
     P('pow',50,'tim',60,'tgel',37,'tw',37,'q',5,'lc',32,'surf',1,'cath',1, ...
       'drv',2,'per',30,'dut',50,'lf',12)
 'cold gel from the bench  --  gel at 25 C', ...
     P('pow',12,'tim',60,'tgel',25,'tw',37,'q',5,'lc',32,'surf',1,'cath',1, ...
       'drv',1,'per',1,'dut',30,'lf',12)
 'write-up baseline  --  perfect heat sink', ...
     P('pow',12,'tim',60,'tgel',37,'tw',37,'q',5,'lc',32,'surf',2,'cath',2, ...
       'drv',1,'per',1,'dut',30,'lf',12)
 'v1/v2 feed  --  truncated current, grey wedge back', ...
     P('pow',50,'tim',60,'tgel',25,'tw',25,'q',5,'lc',32,'surf',1,'cath',1, ...
       'drv',1,'per',1,'dut',30,'lf',0)};

% ======================================================== the window ======
ss = get(groot,'ScreenSize');
W  = min(1420, ss(3)-60);   H = min(860, ss(4)-80);
fig = uifigure('Name','2.45 GHz monopole on Phytagel  --  SAR and temperature', ...
               'Position',[ss(1)+30 ss(2)+40 W H], 'Color',[0.94 0.94 0.96]);
main = uigridlayout(fig,[1 2]);
main.ColumnWidth = {352,'1x'};
main.RowHeight   = {'1x'};
main.Padding     = [8 8 8 8];
main.ColumnSpacing = 10;

ctrlpan = uipanel(main,'Title','Parameters','Scrollable','on','FontWeight','bold');
right   = uigridlayout(main,[4 1]);
right.RowHeight  = {50,'1x',44,96};
right.Padding    = [0 0 0 0];
right.RowSpacing = 4;

% the answer, in words, above the maps
head = uilabel(right,'Text','','FontSize',15,'WordWrap','on', ...
               'HorizontalAlignment','center');

maps = uigridlayout(right,[1 2]);
maps.Padding = [0 0 20 0];  maps.ColumnSpacing = 14;
axS = uiaxes(maps);   axT = uiaxes(maps);

tbar = uigridlayout(right,[1 3]);
tbar.ColumnWidth = {150,'1x',200};
tbar.Padding = [6 0 6 0];
uilabel(tbar,'Text','time through the run','HorizontalAlignment','right');
% The bar is in SECONDS, not in frame number.  It used to be an index into a
% log-spaced list of output times, which meant the first quarter of its travel
% covered the first second of a sixty second run and the last quarter covered
% half a minute: dragging it felt broken because it was, as a scrub bar,
% broken.  GUI_MODEL now returns frames on a fixed increment and the bar is
% linear in time over them, with ticks on round seconds and one minor tick per
% frame, so equal movements are equal times.
tsl = uislider(tbar,'Limits',[0 1],'Value',1,'MajorTicks',[],'MinorTicks',[], ...
               'Enable','off');
tsl.ValueChangedFcn = @(s,e) draw_T();
tlab = uilabel(tbar,'Text','--','FontWeight','bold');

rdo = uitextarea(right,'Editable','off','FontName','Consolas','FontSize',11, ...
                 'Value',{''});

% ==================================================== the control panel ===
gl = uigridlayout(ctrlpan,[40 2]);
gl.ColumnWidth = {150,'1x'};
% Slider rows are 30 and carry no tick labels: the number is in the box beside
% them, so printing the range under every slider was clutter that also clipped.
gl.RowHeight = {24, ...                              % 1      starting point
                18, 24,30, 24,30, 24, 24,30, 24,30, ...  % 2-11   drive
                18, 24,30, 22, 24, 24, ...           % 12-17  gel
                18, 24,30, 24,30, 24, ...            % 18-23  coolant
                18, 24, 24, ...                      % 24-26  bath  (v2)
                18, 24,30, ...                       % 27-29  feed junction
                18, 24,30, 24, ...                   % 30-33  catheter
                24, ...                              % 34     auto / advanced
                0, 0, 0, ...                         % 35-37  advanced, hidden
                30, 30, 26};                         % 38-40  status, buttons
gl.RowSpacing = 2;  gl.Padding = [8 8 8 8];
r = 0;

% ----- starting point -----------------------------------------------------
[dPre, r] = droprow(gl, r, 'start from', PRESET(:,1).', 1, @() apply_preset());
dPre.Tooltip = 'a sensible set of parameters to start from; change anything afterwards';

% ----- drive --------------------------------------------------------------
r = header(gl, r, 'Drive');
[fPow, sPow, r] = pair(gl, r, 'generator power  (W)', 50, [0.1 200], '%.4g', ...
    @() changed(), 'forward power at the generator; the gel receives 75 % of it');
[fTim, sTim, r] = pair(gl, r, 'exposure time  (s)', 60, [0.1 600], '%.4g', ...
    @() changed(), ['how long the run lasts.  With pulsing on this is the ' ...
    'whole run, pulses and gaps together, not the on-time.']);
[dDrv, r] = droprow(gl, r, 'drive', {'continuous','pulsed'}, 1, @() on_drive());
dDrv.Tooltip = ['pulsed: the power above is the power DURING a pulse.  The ' ...
    'mean is duty x that, and both are quoted above the maps.'];
[fPer, sPer, r] = pair(gl, r, 'pulse period  (s)', 1, [0.1 60], '%.4g', ...
    @() changed(), ['one on-time plus one off-time.  0.1 s is 10 Hz, 30 s is ' ...
    'a very slow cycle -- and the slow ones are the ones that hurt.']);
[fDut, sDut, r] = pair(gl, r, 'duty cycle  (%)', 30, [1 100], '%.4g', ...
    @() changed(), ['the fraction of each period the power is on.  100 % is ' ...
    'the same as continuous.']);

% ----- gel ----------------------------------------------------------------
r = header(gl, r, 'Gel');
[fTgel, sTgel, r] = pair(gl, r, 'gel temperature  (C)', 25, [0 40], '%.4g', ...
    @() on_gelT(), ['the temperature the gel STARTS at.  With the box below ' ...
    'ticked it also sets eps'' and sigma.']);
r = r+1;
cLink = uicheckbox(gl,'Text',' tie eps'', sigma to the gel temperature', ...
                   'Value',true,'FontSize',10.5);
cLink.Layout.Row = r;  cLink.Layout.Column = [1 2];
cLink.ValueChangedFcn = @(s,e) on_link();
cLink.Tooltip = ['ticked: both numbers follow the gel temperature through the ' ...
                 'saline model.  You do NOT have to untick this to edit them -- ' ...
                 'just type in either box and the tie releases itself.'];
[fEps, r] = fieldrow(gl, r, 'permittivity  eps''', 74.29, [1.01 120], '%.3f', ...
    @() on_eps(), ['real part of the relative permittivity of the gel.  Type a ' ...
    'measured value in here whenever you like -- it releases the tie to the ' ...
    'gel temperature for you.']);
[fSig, r] = fieldrow(gl, r, 'conductivity  (S/m)', 2.723, [0 20], '%.4f', ...
    @() on_eps(), ['effective conductivity: this is what absorbs the microwave.  ' ...
    'Type a measured value in here whenever you like -- it releases the tie.']);

% ----- coolant ------------------------------------------------------------
r = header(gl, r, 'Coolant  (the blood analogue)');
[fTw, sTw, r] = pair(gl, r, 'water temperature  (C)', 25, [17 60], '%.4g', ...
    @() changed(), 'the temperature the flowing water holds the surface towards');
[fQ, sQ, r] = pair(gl, r, 'flow rate  (L/min)', 5, [0.1 50], '%.4g', ...
    @() changed(), ['sets how hard the surface is cooled.  Only does anything ' ...
    'with realistic flow cooling selected below.']);
[dSurf, r] = droprow(gl, r, 'surface cooling', ...
      {'realistic flow cooling', 'perfect heat sink', 'none'}, 1, @() changed());
dSurf.Tooltip = ['realistic: h computed from the flow rate.  Perfect heat sink: ' ...
    'the surface is pinned at the water temperature -- what the brief specified, ' ...
    'and the LOWEST temperature the gel can reach.'];

% ----- bath (v2) ----------------------------------------------------------
% This is the whole point of v2.  v1 had AIR above the gel, so essentially all
% the delivered power went into the gel; the bench photographs show the
% antenna immersed, and a bath whose permittivity is close to the gel's own
% takes most of it instead.  Leaving these two controls out would ship a
% window that looks exactly like v1 and cannot show what changed.
r = header(gl, r, 'Bath  (what the antenna sits in)');
[dUp, r] = droprow(gl, r, 'above the gel', ...
      {'water  (v2)', 'air  (v1 behaviour)'}, 1, @() on_upper());
dUp.Tooltip = ['water: the antenna is immersed, and the bath absorbs its own ' ...
    'share.  air: reproduces v1 exactly, where 99.6 % of the power reached ' ...
    'the gel.'];
[fSalt, r] = fieldrow(gl, r, 'bath salt  (% w/v)', 0, [0 5], '%.2f', ...
    @() changed(), ['0 = deionised or tap, 0.9 = physiological saline.  Pure ' ...
    'water is already lossy at 2.45 GHz, so this matters less than you would ' ...
    'expect.']);

% ----- feed ---------------------------------------------------------------
% The reason the maps no longer have a grey semicircle punched out of them.
% The prescribed antenna current used to stop dead at the feed, which put a
% point charge there and made the deposition -- and so the temperature --
% singular at that one node; the only defence was to blank a sphere around it
% and quote the peak outside.  The current does not stop, it RETURNS along the
% bare braid, and with that arm in place there is nothing left to blank.
r = header(gl, r, 'Feed  (the coax braid)');
[fLf, sLf, r] = pair(gl, r, 'braid arm  (mm)', 12, [0 30], '%.4g', ...
    @() changed(), ['how much of the bare braid behind the feed carries the ' ...
    'RETURN current.  12 mm makes it the mirror of the element, which is what ' ...
    'a dipole with no ground plane and no balun does.  Shorter piles charge ' ...
    'up at the junction; 0 = the old truncated current, singular at the feed, ' ...
    'so a 2 mm sphere gets blanked and the peak stops meaning anything.']);

% ----- catheter -----------------------------------------------------------
r = header(gl, r, 'Catheter');
[fLc, sLc, r] = pair(gl, r, 'length on the gel  (mm)', 32, [0 33], '%.4g', ...
    @() changed(), ['how much catheter lies on the gel, measured back from the ' ...
    'antenna tip.  Only does anything if it is cooled or insulating.']);
[dCath, r] = droprow(gl, r, 'catheter is', ...
      {'cooled by the flow', 'ignored', 'blocking the flow'}, 1, @() changed());
dCath.Tooltip = ['cooled: it gets its own film coefficient, which is lower than ' ...
    'the bare surface, so a longer catheter runs slightly hotter.'];

% ----- auto-update and advanced ------------------------------------------
r = r+1;
orow = uigridlayout(gl,[1 2]);
orow.Layout.Row = r;  orow.Layout.Column = [1 2];
orow.Padding = [0 0 0 0];  orow.ColumnSpacing = 4;
cAuto = uicheckbox(orow,'Text',' update as I change things','Value',true, ...
                   'FontSize',10.5);
cAdv  = uicheckbox(orow,'Text',' advanced','Value',false,'FontSize',10.5, ...
                   'ValueChangedFcn',@(s,e) on_advanced());
cAuto.Tooltip = ['off: nothing recomputes until you press the button.  Turned off ' ...
                 'automatically if you pick the slow 0.25 mm grid.'];

[dGrid, r] = droprow(gl, r, 'grid', {'0.5 mm  (a few seconds)', ...
       '1.0 mm  (fast, coarse)', '0.25 mm  (slow, the baseline grid)'}, 1, ...
       @() on_grid());
[dSar, r]  = droprow(gl, r, 'SAR scale', {'log10','linear'}, 1, @() draw_SAR());
[dTsc, r]  = droprow(gl, r, 'T scale', {'auto','fixed 20-100 C'}, 1, @() draw_T());
advrows = [35 36 37];

% ----- status and buttons -------------------------------------------------
r = r+1;
stat = uilabel(gl,'Text','','WordWrap','on','FontSize',11);
stat.Layout.Row = r;  stat.Layout.Column = [1 2];
r = r+1;
bGo = uibutton(gl,'Text','update now','FontWeight','bold','FontSize',12, ...
               'BackgroundColor',[0.20 0.45 0.75],'FontColor','w', ...
               'ButtonPushedFcn',@(s,e) compute());
bGo.Layout.Row = r;  bGo.Layout.Column = [1 2];
r = r+1;
brow = uigridlayout(gl,[1 2]);
brow.Layout.Row = r;  brow.Layout.Column = [1 2];
brow.Padding = [0 0 0 0];  brow.ColumnSpacing = 6;
uibutton(brow,'Text','start over','ButtonPushedFcn',@(s,e) apply_preset());
uibutton(brow,'Text','save PNG','ButtonPushedFcn',@(s,e) save_png());

on_upper();                 % starting = true, so this only greys the salt box
on_drive();                 % ... and this only greys the pulse boxes
on_link();                  % starting = true, so this only fills the two boxes
starting = false;
stat.Text = 'first solve...';
drawnow;
compute();

% =========================================================================
%                               callbacks
% =========================================================================
    function changed()
        % one entry point for every parameter: recompute, or say it is stale
        if starting, return; end
        if busy
            pending = true;         % compute() will notice on the way out
        elseif cAuto.Value
            compute();
        else
            dirty();
        end
    end

    function on_gelT()
        if cLink.Value, show_stogryn(); end
        changed();
    end

    function on_link()
        % The two boxes stay editable whatever this checkbox says.  Locking
        % them was the old behaviour and it read as "these numbers are dead":
        % the only way in was a tick box whose label does not say "untick me
        % to edit", so people reasonably concluded the fields were broken.
        % The tick now decides ONE thing -- whether the pair is recomputed
        % from the gel temperature -- and typing releases it (see ON_EPS).
        %
        % Leaving them locked also desynced: APPLY_PRESET re-ticks the box
        % without coming through here, so the fields stayed editable while
        % the tie was on and anything typed was silently overwritten.
        if cLink.Value, show_stogryn(); end
        changed();
    end

    function on_eps()
        % Typing a measured eps' or sigma is an explicit override.  Release
        % the tie, rather than silently discarding what was just typed the
        % next time the gel temperature moves.
        cLink.Value = false;
        changed();
    end

    function on_upper()
        % with air above there is no bath to salt
        on = strcmp(UPPER{idx(dUp)},'water');
        fSalt.Editable = on;
        if on, fSalt.BackgroundColor = [1 1 1];
        else,  fSalt.BackgroundColor = [0.94 0.94 0.94]; end
        changed();
    end

    function on_drive()
        grey_pulse();
        changed();
    end

    function grey_pulse()
        % continuous means the period and the duty say nothing, so grey them
        % rather than leaving two live-looking boxes that change no answer.
        % Kept separate from ON_DRIVE so APPLY_PRESET can re-grey without
        % triggering a second solve on top of its own.
        en = 'off';
        if idx(dDrv) == 2, en = 'on'; end
        fPer.Enable = en;  sPer.Enable = en;
        fDut.Enable = en;  sDut.Enable = en;
    end

    function on_grid()
        % the 0.25 mm grid takes a minute and a half; auto-update would make
        % the window feel broken, so turn it off and say why
        if idx(dGrid) == 3 && cAuto.Value
            cAuto.Value = false;
            stat.FontColor = [0.75 0.35 0];
            stat.Text = 'auto-update off: the 0.25 mm grid takes ~90 s. Press update now.';
            return
        end
        changed();
    end

    function on_advanced()
        h = 24*double(cAdv.Value);
        for k = advrows, gl.RowHeight{k} = h; end
    end

    function apply_preset()
        p = PRESET{idx(dPre),2};
        setpair(fPow,sPow,p.pow);    setpair(fTim,sTim,p.tim);
        setpair(fTgel,sTgel,p.tgel); setpair(fTw,sTw,p.tw);
        setpair(fQ,sQ,p.q);          setpair(fLc,sLc,p.lc);
        setpair(fPer,sPer,p.per);    setpair(fDut,sDut,p.dut);
        setpair(fLf,sLf,p.lf);
        cLink.Value = true;
        dSurf.Value = dSurf.Items{p.surf};
        dCath.Value = dCath.Items{p.cath};
        dDrv.Value  = dDrv.Items{p.drv};
        grey_pulse();
        show_stogryn();
        changed();
    end

    function show_stogryn()
        [e, gi] = saline_stogryn(2.45e9, min(max(fTgel.Value,0),40), 0.9*10/58.443);
        fEps.Value = real(e);
        fSig.Value = gi.sigma_eff_Spm;
    end

    function dirty()
        stat.FontColor = [0.75 0.35 0];
        stat.Text = 'the maps are out of date -- press update now';
    end

    function compute()
        if busy, return; end
        busy = true;
        bGo.Enable = 'off';
        stat.FontColor = [0 0 0];  stat.Text = 'solving...';  drawnow;
        try
            if idx(dDrv) == 2
                per = fPer.Value;   dut = fDut.Value/100;
            else
                per = [];           dut = 1;
            end
            R = gui_model('eps_r',fEps.Value, 'sigma',fSig.Value, ...
                          'link_eps',cLink.Value, ...
                          'T_gel',fTgel.Value, 'T_water',fTw.Value, ...
                          'Qflow_Lmin',fQ.Value, 'L_cath_mm',fLc.Value, ...
                          'Pgen',fPow.Value, 't_stim',fTim.Value, ...
                          'pulse_period',per, 'duty',dut, ...
                          'L_sheath_mm',fLf.Value, ...
                          'surf_bc',SURF{idx(dSurf)}, 'cath_bc',CATH{idx(dCath)}, ...
                          'upper',UPPER{idx(dUp)}, 'bath_NaCl_pct',fSalt.Value, ...
                          'dh',DHV(idx(dGrid)), 'em',em);
            em = R.em;
            if cLink.Value
                fEps.Value = R.info.eps_r;   fSig.Value = R.info.sigma;
            end
            set_timebar();
            draw_SAR();  draw_T();
            stat.FontColor = [0.1 0.4 0.1];
            stat.Text = sprintf('up to date  (%.1f s)', R.info.t_total);
        catch err
            stat.FontColor = [0.8 0 0];
            stat.Text = ['FAILED: ' err.message];
            head.Text = 'the solve failed';
            rdo.Value = {err.message};
        end
        bGo.Enable = 'on';
        busy = false;
        if pending          % a control moved while we were solving
            pending = false;
            if cAuto.Value, compute(); else, dirty(); end
        end
    end

    function set_timebar()
        %SET_TIMEBAR  Re-range the time bar in SECONDS for the run just solved.
        %   The ticks are cleared first: a leftover tick outside the new limits
        %   is an error, and the exposure time can shrink between solves.
        tend = R.t(end);
        tsl.MajorTicks = [];  tsl.MinorTicks = [];
        tsl.Limits = [0 tend];
        tsl.Value  = tend;
        tk = 0:nice_tick(tend):tend;
        if tk(end) < tend, tk(end+1) = tend; end
        tsl.MajorTicks      = tk;
        tsl.MajorTickLabels = arrayfun(@(v) sprintf('%g',v), tk, 'uni',0);
        tsl.MinorTicks      = R.t;          % one per frame -- the increments
        tsl.Enable = 'on';
    end

    function draw_SAR()
        if isempty(R), return; end
        [xm, sm, msk] = geom();
        A = R.SAR;
        if strcmp(dSar.Value,'log10')
            top  = log10(max(A(msk)));
            V    = log10(max(A, 10^(top-4)));
            lims = [top-4 top];   ttl = 'log_{10} SAR   (W/kg)';
        else
            V = A;  lims = [0 max(A(msk))];   ttl = 'SAR   (W/kg)';
        end
        cla(axS);
        imagesc(axS, xm, sm, V.', 'AlphaData', double(msk.'));
        axS.YDir = 'reverse';  axS.Color = [0.45 0.45 0.50];
        clim(axS, lims);  colormap(axS, turbo);
        hold(axS,'on');  decorate(axS,'w');  hold(axS,'off');
        cb = colorbar(axS);  cb.Label.String = ttl;
        frame(axS, xm, sm);
        xlabel(axS,'x  (mm)');  ylabel(axS,'depth below the surface  (mm)');
        title(axS, sprintf('where the power goes  --  %.4g W in, %.4g W delivered', ...
                           R.opt.Pgen, R.info.Pdel));
        report();
    end

    function draw_T()
        if isempty(R), return; end
        [~, it] = min(abs(R.t - tsl.Value));
        tsl.Value = R.t(it);            % snap to the frame actually drawn
        [xm, sm, msk] = geom();
        A = R.T_plane(:,:,it);
        B = A;  B(~msk) = NaN;
        if strcmp(dTsc.Value,'auto')
            lo = min(R.opt.T_water, R.opt.T_gel);
            hi = max(max(B(:)), lo+1);
        else
            lo = 20;  hi = 100;
        end
        cla(axT);
        imagesc(axT, xm, sm, A.', 'AlphaData', double(msk.'));
        axT.YDir = 'reverse';  axT.Color = [0.45 0.45 0.50];
        clim(axT, [lo hi]);  colormap(axT, hot);
        hold(axT,'on');
        lv = [43 50 60 80 100];  lv = lv(lv > lo & lv < hi);
        if ~isempty(lv)
            [cc,hh] = contour(axT, xm, sm, B.', lv, 'w-','LineWidth',0.9);
            if ~isempty(cc), clabel(cc,hh,'Color','w','FontSize',7); end
        end
        decorate(axT,'c');  hold(axT,'off');
        cb = colorbar(axT);  cb.Label.String = 'T   (\circC)';
        frame(axT, xm, sm);
        xlabel(axT,'x  (mm)');  ylabel(axT,'depth below the surface  (mm)');
        title(axT, sprintf('how hot it gets  --  t = %.3g s', R.t(it)));
        tlab.Text = sprintf('t = %.3g s  of  %.4g s', R.t(it), R.opt.t_stim);
        report();
    end

    function [xm, sm, msk] = geom()
        xm = R.x*1e3;  sm = -R.z*1e3;  msk = R.pmask_y0;
    end

    function frame(ax, xm, sm)
        %FRAME  Square pixels, and the axes cropped to the data.
        %   The limits are set EXPLICITLY.  axis('image') leaves a margin of
        %   background either side once a colorbar has taken its width, and the
        %   axes background is the grey that means "blanked" -- so an empty
        %   margin read as a blanked strip at the edge of the block, which is
        %   exactly the confusion the feed wedge used to cause.  Cropping to
        %   the data means every grey pixel in the maps is a real exclusion.
        axis(ax,'image');
        xlim(ax, [xm(1) xm(end)]);
        ylim(ax, [sm(1) sm(end)]);
    end

    function decorate(ax, barcol)
        % The dashed semicircle is the blanked feed sphere.  It only exists
        % when the braid arm is off; with the arm on nothing is blanked, so
        % drawing a circle there would claim a boundary that is not in the
        % answer.
        rex = 1e3*R.Tc.r_feed_excl;
        if rex > 0
            th = linspace(0,pi,120);
            plot(ax, rex*cos(th), rex*sin(th), 'w--','LineWidth',1);
        end
        % the braid arm, drawn back from the feed at half the weight of the
        % element so the driven arm and the return arm are not confused
        if R.info.L_sheath > 0
            plot(ax, [-1e3*R.info.L_sheath 0], [0 0], barcol, 'LineWidth',2);
        end
        if ~strcmp(CATH{idx(dCath)},'transparent')
            x0 = max(1e3*R.Tc.cath_x0, R.x(1)*1e3);
            plot(ax, [x0 1e3*R.Tc.cath_L], [0 0], 'w:','LineWidth',3);
        end
        plot(ax, [0 1e3*R.Tc.L_ant], [0 0], barcol, 'LineWidth',4);
    end

    function report()
        if isempty(R), return; end
        i = R.info;  o = R.opt;

        % ---- the answer, in one sentence ------------------------------
        % judge on the WORST of the run, not on where it happened to finish
        Tv = max(i.T_end, i.T_worst);
        if Tv >= 100
            verdict = 'that is past boiling -- this power cannot be used';
        elseif Tv >= 60
            verdict = 'hot enough to ablate, and past where the gel may soften';
        elseif Tv >= 43
            verdict = 'above the 43 C damage threshold';
        else
            verdict = 'stays below the 43 C damage threshold';
        end
        if i.pulsed
            drv = sprintf('%.4g W pulsed at %.4g %% duty, %.4g s period (%.4g W mean)', ...
                          o.Pgen, 100*i.duty, i.period, i.Pgen_mean);
        else
            drv = sprintf('%.4g W continuous', o.Pgen);
        end
        peaktxt = sprintf('peak %.0f %sC', i.T_end, deg);
        if i.T_worst > i.T_end + 0.5
            % a pulsed run can end in an off-period, so say what the worst was
            peaktxt = sprintf('peak %.0f %sC at the end, %.0f %sC at its worst', ...
                              i.T_end, deg, i.T_worst, deg);
        end
        head.Text = sprintf(['%s for %.4g s  ->  %s,  %.0f %sC at the ' ...
            'probe point (6 mm along, 2 mm deep).  %s.'], ...
            drv, o.t_stim, peaktxt, i.T_probe_end, deg, verdict);
        if Tv >= 100,     head.FontColor = [0.75 0.10 0.10];
        elseif Tv >= 43,  head.FontColor = [0.70 0.40 0.00];
        else,             head.FontColor = [0.10 0.45 0.15];
        end

        % ---- the numbers behind it ------------------------------------
        L = {};
        L{end+1} = sprintf(['gel     eps'' %.2f, sigma %.3f S/m (%s); half the power ' ...
            'lands above %.1f mm deep, nine tenths above %.1f mm'], ...
            i.eps_r, i.sigma, i.eps_src, 1e3*i.d50, 1e3*i.d90);
        L{end+1} = sprintf(['cooling h = %.0f W/(m^2 K) at the surface, %.0f under the ' ...
            'catheter; the water itself would warm %.2f K'], ...
            i.h_surf, i.h_cath, i.dT_coolant);
        L{end+1} = sprintf(['power   %.1f %% absorbed in the gel, %.1f %% of that in the ' ...
            '30 mm block (%.2f W); peak SAR %.3g W/kg'], ...
            100*i.frac_gel, 100*i.frac_block, i.P_in, i.SAR_max);
        L{end+1} = sprintf(['bath    %.1f %% of the delivered power is absorbed in the ' ...
            'bath (%.2f W), eps'' %.1f; with air above it was 0'], ...
            100*i.frac_bath, i.Pwat, real(i.epsr_upper));
        if i.L_sheath > 0
            L{end+1} = sprintf(['feed    the return current runs %.1f mm back ' ...
                'along the braid, so nothing is blanked and the peak below ' ...
                'is a prediction'], 1e3*i.L_sheath);
        else
            L{end+1} = ['feed    truncated current (v1/v2): singular at the ' ...
                'feed, so a 2 mm sphere is blanked and the peak below sits ' ...
                'on its edge'];
        end
        L{end+1} = sprintf(['peak    %.1f %sC at x = %.1f mm, %.2f mm deep ' ...
            '(worst over the whole run %.1f %sC)'], ...
            i.T_end, deg, 1e3*i.xyz_end(1), -1e3*i.xyz_end(3), i.T_worst, deg);
        if i.pulsed
            L{end+1} = sprintf(['drive   %.4g W on for %.3g s, off for %.3g s, ' ...
                '%.0f cycles; %.4g W mean into the generator, %.2f W mean into ' ...
                'the gel'], o.Pgen, i.duty*i.period, (1-i.duty)*i.period, ...
                floor(o.t_stim/i.period), i.Pgen_mean, i.P_in_mean);
        end
        for k = 1:numel(i.flags)
            L{end+1} = ['  !     ' i.flags{k}]; %#ok<AGROW>
        end
        rdo.Value = L(:);
    end

    function save_png()
        if isempty(R), return; end
        stamp = char(datetime('now','Format','yyyyMMdd_HHmmss'));
        try
            fn = fullfile(pick_figdir(), ['gui_' stamp '.png']);
            % exportapp, not exportgraphics: it captures the whole window,
            % controls included, so the image records the parameters that
            % produced the maps.
            exportapp(fig, fn);
            stat.FontColor = [0.1 0.4 0.1];
            stat.Text = ['saved  ' fn];
        catch err
            stat.FontColor = [0.8 0 0];
            stat.Text = ['could not save: ' err.message];
        end
    end

    function d = pick_figdir()
        %PICK_FIGDIR  The first writable folder to put a saved PNG in.
        %   A shared copy of this model may sit somewhere that cannot be
        %   written to -- a network share, a read-only synced folder, inside
        %   an installed toolbox -- so the folder the code was launched from
        %   is a preference, not an assumption.  Falling back keeps "save PNG"
        %   working instead of failing in someone else's hands.
        persistent cached
        if ~isempty(cached), d = cached;  return; end
        for c = {fullfile(root,'figures'), ...
                 fullfile(userpath,'PhytagelHeatGUI'), ...
                 fullfile(tempdir,'PhytagelHeatGUI')}
            d = c{1};
            if isempty(d) || isempty(fileparts(d)), continue; end
            try
                if ~exist(d,'dir'), mkdir(d); end
                t   = fullfile(d,'.hm_writetest');
                fid = fopen(t,'w');
                if fid < 0, continue; end
                fclose(fid);  delete(t);
                cached = d;  return
            catch
                % not this one -- try the next
            end
        end
        error('heat_gui:savedir', ...
              'no writable folder found for the PNG (tried the model folder, %s and %s)', ...
              userpath, tempdir);
    end
end

% =========================================================================
%                          small layout helpers
% =========================================================================
function d = nice_tick(T)
%NICE_TICK  A round tick spacing giving about six labels across T seconds.
raw = T/6;
p   = 10^floor(log10(raw));
m   = raw/p;
if     m <= 1, d = 1*p;
elseif m <= 2, d = 2*p;
elseif m <= 5, d = 5*p;
else,          d = 10*p;
end
end

function r = header(gl, r, txt)
r = r+1;
h = uilabel(gl,'Text',txt,'FontWeight','bold','FontColor',[0.15 0.30 0.55]);
h.Layout.Row = r;  h.Layout.Column = [1 2];
end

function [f, s, r] = pair(gl, r, txt, v, lim, fmt, cb, tip)
%PAIR  label + box on one row, a slider for it on the next, kept in step.
[f, r] = fieldrow(gl, r, txt, v, lim, fmt, [], tip);
r = r+1;
s = uislider(gl,'Limits',lim,'Value',v,'MajorTicks',[],'MinorTicks',[]);
s.Layout.Row = r;  s.Layout.Column = [1 2];
f.ValueChangedFcn = @(o,e) sync(f, s, f.Value, cb);
s.ValueChangedFcn = @(o,e) sync(f, s, s.Value, cb);
end

function sync(f, s, v, cb)
v = min(max(v, s.Limits(1)), s.Limits(2));
f.Value = v;  s.Value = v;
if ~isempty(cb), cb(); end
end

function setpair(f, s, v)
v = min(max(v, s.Limits(1)), s.Limits(2));
f.Value = v;  s.Value = v;
end

function [f, r] = fieldrow(gl, r, txt, v, lim, fmt, cb, tip)
r = r+1;
l = uilabel(gl,'Text',txt,'HorizontalAlignment','right','FontSize',11);
l.Layout.Row = r;  l.Layout.Column = 1;
f = uieditfield(gl,'numeric','Value',v,'Limits',lim,'ValueDisplayFormat',fmt);
f.Layout.Row = r;  f.Layout.Column = 2;
if nargin >= 7 && ~isempty(cb), f.ValueChangedFcn = @(o,e) cb(); end
if nargin >= 8 && ~isempty(tip), l.Tooltip = tip;  f.Tooltip = tip; end
end

function [d, r] = droprow(gl, r, txt, items, k, cb)
r = r+1;
l = uilabel(gl,'Text',txt,'HorizontalAlignment','right','FontSize',11);
l.Layout.Row = r;  l.Layout.Column = 1;
d = uidropdown(gl,'Items',items,'Value',items{k},'FontSize',10.5);
d.Layout.Row = r;  d.Layout.Column = 2;
if ~isempty(cb), d.ValueChangedFcn = @(o,e) cb(); end
end

function k = idx(d)
k = find(strcmp(d.Items, d.Value), 1);
end
