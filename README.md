# 2.45 GHz monopole on Phytagel — field and heating models

Two coupled models of a 12 mm 2.45 GHz monopole lying on a 30 mm Phytagel slab, with 37 °C blood
flowing over the surface.

| start here | what it answers |
|---|---|
| **[`E_field_model_writeup.md`](E_field_model_writeup.md)** | where the electric field goes |
| **[`Heatmodel/Heat_model_writeup.md`](Heatmodel/Heat_model_writeup.md)** | how hot the gel gets, and how much power you may use |
| **[`Heatmodel/GUI_instructions.md`](Heatmodel/GUI_instructions.md)** | how to drive `heat_gui` — the interactive parameter sweep |
| **[`share/`](share/)** | a self-contained copy of the GUI to send to other people |
| **[`v2_water_above/README.md`](v2_water_above/README.md)** | the corrections since: the antenna is in the bath, and the feed current returns along the braid |
| **[`Gel_model_V3/E_field_2d_notes.md`](Gel_model_V3/E_field_2d_notes.md)**, **[`Gel_model_V3/gel_heating_notes.md`](Gel_model_V3/gel_heating_notes.md)** | V3 (Sep 2026): closed-form field in one lossy medium, 3-D conduction with convective faces, GUI + surface-sparing optimiser — and why the pattern differs from the 50 W photographs (12 mm element, insulated-wire wavenumber, surface pinning) |

> **The GUI in `share/Phytagel Heat GUI v2/` is ahead of this document and of both write-ups.**
> It puts the antenna in the **bath** rather than in air, which is what the photographs show and
> which moves most of the delivered power out of the gel; it carries the feed current **back along
> the coax braid** instead of truncating it, which removes the feed singularity, the blanked grey
> wedge with it, and about 20 % of the peak; and it can **pulse** the drive. Every heating number
> below is from the air-above model and reads several times too high.
> `v2_water_above/README.md` has the physics, the verification and what is still unexplained.

Each document carries its own derivation, every constant with its source, a literature screen,
validation against analytic solutions, results and limitations.

## The two answers in one line each

**Field.** The magnitude ripples at $2f$ between the semi-axes of the local polarisation ellipse:

$$|\mathbf E(x,z,t)| = \sqrt{\tfrac12\Big[\big|\tilde{\mathbf E}(x,z)\big|^{2} + |Q(x,z)|\cos(2\omega t + \arg Q)\Big]},
\qquad Q = \tilde E_x^2+\tilde E_y^2+\tilde E_z^2$$

$\tilde{\mathbf E}$ has no elementary closed form and is computed numerically.

**Heating.** The heat equation is linear and the deposition is proportional to delivered power, so
one solve answers every power exactly:

$$T(\mathbf r,t) = 37\ \text{°C} + \frac{\eta\,P_\text{gen}}{150\ \text{W}}\,\Theta(\mathbf r,t),
\qquad \eta = \frac{P_\text{del}}{P_\text{gen}} = 0.75$$

Note $|\mathbf E|\propto\sqrt{P}$ but $\Delta T\propto P$, because heating goes with $E^2$.

## Quick start

```matlab
addpath('matlab'); addpath('Heatmodel/matlab');
run_model                            % ~30 s   -> results/field_model.mat + 9 figures
run_heat_model                       % ~15 min -> Heatmodel/results/heat_model.mat + 8 figures

E_gel(6e-3, -5e-3)                   % time-RMS |E| at x = 6 mm, 5 mm deep        -> V/m
T_gel(12, 60)                        % peak gel temperature, 12 W generator, 60 s -> 97.7 °C
T_gel(12, 60, 6e-3, -2e-3)           % at a point on the y = 0 plane              -> °C
P_for_T(100, 60)                     % power that just boils it in a minute       -> 12.46 W
time_to_T(50, 100)                   % seconds until 100 °C at 50 W               ->  3.39 s

% one bench run -- power [W], pulse period [s] ([] = continuous), stimulation time [s]
heat_scenario(12, [], 60)                        % draws the figures, writes a CSV
heat_scenario(40,  1, 60, 'duty', 0.3)           % 40 W, 30 % duty at 1 Hz
heat_scenario(12, [], 60, 'model', 'robin')      % convective surface, for the bracket
heat_scenario(12, [], 60, 'measured', myLog)     % overlay a thermocouple record

% sweep the PHYSICS, not just the drive -- a window with the SAR and temperature maps
heat_gui                                         % ~5 s per solve
```

`heat_gui` is the one to use for asking *what if the gel is not what we assumed?* Everything else
in this repository reuses one stored solution and so takes the physics as fixed; `heat_gui`
re-solves the whole chain, electromagnetics included, and shows the SAR map and the temperature map
side by side. Editable: **gel permittivity, gel conductivity, gel temperature, water temperature,
water flow rate, catheter length**, plus power, exposure time and the grid. The v2 window in
`share/` adds **the medium above the gel and its salinity**, **the length of braid carrying the
return current**, and **pulsing** — period and duty cycle.

Pick a starting point from the top drop-down and move a slider — the maps redraw themselves, in
about a second for a change of power or time and about five if the gel changed, and the line above
them gives the answer in words. It opens on a convective surface with a cooled catheter rather than
this document's
Dirichlet baseline, because under the baseline the flow rate and the catheter length have no effect
at all — both baselines stay in the drop-downs. What the speed costs is measured, not assumed:
`verify_gui_model` puts it at 1.1 % of the temperature rise against the Dirichlet baseline and
3.8 % against the Robin one (heat writeup §11.5).

`heat_scenario` is the one to use for validating against a bench measurement: it writes
`figures/scenario_*_traces.png` (drive, peak and probe traces with the cool-down, depth profiles,
summary), `figures/scenario_*_maps.png` (temperature maps through the run plus thermal dose) and
`results/scenario_*.csv` (the traces, ready to plot next to a measurement). It needs no new solve —
every scenario is a superposition of the one stored step response — so a run takes under a second.

**Units are SI.** $x$ runs along the antenna (feed at $x=0$, tip at $x=12$ mm); $z$ is height, so
**$z$ is negative inside the gel** (surface $z=0$, floor $z=-30$ mm). `T_gel` takes the **generator
forward power**; the gel sees $\eta P$.

## Sending the GUI to someone else

`share/Phytagel Heat GUI.zip` (103 kB) is the whole GUI, standalone. They need MATLAB R2022b or
newer and nothing else — no toolboxes, no downloads, and nothing to type. They unzip it, read
`START HERE.txt`, and pick one of two routes:

- double-click **`Phytagel Heat GUI.mltbx`** — installs into MATLAB, after which the GUI is a button
  in the **Apps** tab (labelled `heat_gui`) for good; or
- double-click **`Open the GUI (Windows).bat`** / **`(Mac).command`** — starts MATLAB and opens the
  window, installing nothing.

The zip carries only the GUI's 13-file dependency closure, so `T_gel`, `heat_scenario` and the rest
of the project are *not* in it: those reuse the stored quarter-hour solve, which the GUI does not
need. Rebuild the zip after changing the model — nothing does that automatically.

## Layout

| path | what |
|---|---|
| `E_field_model_writeup.md` | **the field document** — theory, sources, validation, results, limitations |
| `matlab/E_gel.m` | the field deliverable |
| `matlab/layered_field.m` | the field solver (exact 3-layer spectral / Sommerfeld, FFT-based) |
| `matlab/model_config.m` | every electromagnetic constant, each with its source |
| `matlab/verify_layered_field.m` | field validation suite |
| `figures/`, `results/` | field figures and outputs |
| `Heatmodel/Heat_model_writeup.md` | **the heating document** |
| `Heatmodel/matlab/T_gel.m` | **the heating deliverable** — power and time in, temperature out |
| `Heatmodel/matlab/heat_scenario.m` | **one bench run** — power, pulse period, stimulation time in; figures and a CSV out |
| `Heatmodel/matlab/heat_gui.m` | **the parameter sweep** — a window showing the SAR and temperature maps for gel permittivity, conductivity and temperature, water temperature and flow rate, and catheter length |
| `Heatmodel/matlab/gui_model.m` | the fast forward solve behind it (headless, scriptable) |
| `Heatmodel/matlab/verify_gui_model.m` | what the fast solve costs against the stored baseline — 12 checks |
| `Heatmodel/matlab/P_for_T.m`, `time_to_T.m` | the inverses: power for a limit, time to a limit |
| `Heatmodel/matlab/sar_volume.m` | 3-D power deposition, from the field solver |
| `Heatmodel/matlab/heat_solve.m`, `heat_steady.m` | the transient and steady conduction solvers |
| `Heatmodel/matlab/thermal_config.m` | every thermal constant and boundary condition, each sourced |
| `Heatmodel/matlab/verify_heat_model.m` | heating validation suite — **run this first** |
| `Heatmodel/figures/`, `Heatmodel/results/` | heating figures and outputs |
| `Heatmodel/matlab/hm_addpath.m` | finds both code folders wherever they have been copied to |
| `share/` | the shareable GUI package and the zip built from it |
| `OLD_dont_use/` | earlier notes — left untouched |

## Headline results

### Field (150 W delivered, 0.9 % saline Phytagel, 25 °C)

- $\varepsilon_r = 74.29 - j19.98$, $\sigma_\text{eff} = 2.72$ S/m, $\lambda_\text{gel} = 14.07$ mm
- **99.7 % of the delivered power is absorbed in the gel**
- Field falls to $1/e$ by **2.7 mm** and $1/10$ by **12.1 mm** below the antenna — much shallower
  than the 16.9 mm plane-wave skin depth, because the near field is spreading-dominated
- 96 % of in-plane power deposition is above 10 mm depth

### Heating (surface held at 37 °C, as specified)

- At **200 W generator** the gel reaches 43 °C in **0.04 s** and 100 °C in **0.72 s**.
  **The phantom cannot absorb anything like full power.**
- **12.4 W** gives 100 °C after a 60 s exposure; **10.2 W** can be run indefinitely; **1.0 W**
  holds the peak below 43 °C for ever
- The peak rises **5.07 K per generator watt at 60 s**, and reaches 82 % of its final value by
  one minute
- The hottest gel is **2–4 mm below the surface, never at it** — the blood pins the surface and the
  catheter, cooled by the same flow, adds no hot spot of its own
- At 12.4 W for 60 s the ablation-equivalent zone (CEM43 ≥ 240 min) is **9 mm deep, 15 mm long**
- **Pulsing does not buy headroom.** At the same mean power, 30 % duty at 1 Hz raises the peak by
  2 °C and at 0.02 Hz by 79 °C

## The four things to measure before trusting the numbers

1. **The standoff gap** between the conductor and the gel surface (field writeup §9B). It is the
   largest sensitivity in the field model — and it matters *twice as much* thermally, because
   $\Delta T\propto|\mathbf E|^2$.
2. **The flow geometry over the phantom** (heat writeup §3, §10B). "Surface held at 37 °C" is the
   perfect-heat-sink limit; a realistic film coefficient costs up to a factor 2.6 in allowed power.
   Put the gel at the leading edge of a fast, thin channel and the assumption becomes true.
3. **The gel permittivity**, with an open-ended coaxial probe if you can get access (field
   writeup §4, §11). It is currently inferred from a validated saline model, not measured.
4. **The gel–sol transition of the Phytagel itself** (heat writeup §12.5). Low-acyl gellan gels are
   reported to soften between 30 and 50 °C at low ionic strength. The phantom may stop being a
   phantom well below the temperatures these models report.

## Known artefact you should be aware of

The field model prescribes the antenna current as $\sin[k_e(L-x)]$, which is at its **maximum**
where it is truncated at the feed. That step puts a fictitious point charge at $(0,0,h)$ and makes
the deposition — and hence the temperature — singular there. It never mattered for the field maps,
which were quoted away from the feed; it matters for a temperature model, whose natural output is a
maximum. Diagnosis, evidence and treatment are in `Heatmodel/matlab/feed_regularise.m` and heat
writeup §6.2. Everything more than about 3 mm from the feed is unaffected to better than 1 %.
A real feed junction *is* a hot spot — **instrument it on the bench.**
