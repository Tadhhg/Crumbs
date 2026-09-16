> ## ⚠ Read this first — this is the **v2.1** window
>
> This document was written for **v1**. The controls, the layout and everything about *driving* the
> window are still right; the **temperatures quoted in the body are v1 numbers and read several
> times too high**. §5's preset table has been brought up to date; the rest has not.
>
> Three things changed since it was written.
>
> **v2 — the antenna is in the bath.** v1 put **air** above the gel and so sent 99.6 % of the
> delivered power into it. The bench photographs show the antenna fully immersed, and a bath of
> nearly the same permittivity takes most of the power instead. New: a **Bath** section (`above the
> gel`, `bath salt (% w/v)`) and a **bath** line in the report. Setting it to **air** reproduces v1.
>
> **v2.1 — the feed current returns along the braid.** The grey wedge at the feed is gone: §4 below
> describes a blanked region that no longer exists, and has been rewritten. New: a **Feed** section
> (`braid arm (mm)`). Setting it to **0** restores the v1/v2 truncated current, wedge and all.
>
> **v2.1 — the drive can be pulsed, and the time bar is linear in seconds.** New: `drive`,
> `pulse period (s)`, `duty cycle (%)`, and a time bar labelled in seconds whose steps are equal.
>
> On the opening preset (50 W, 60 s, everything at 25 °C) v2.1 gives **peak 117 °C, probe 81 °C**,
> against v2's 154 / 128 and v1's 400 / 239.
>
> `What changed in v2.md` has the physics, the verification checks, and what is still **not**
> explained.

# `heat_gui` — how to use it

A window with two maps: **where the microwave power goes** (SAR) and **how hot the gel gets**. Change
the gel, the coolant, the catheter or the drive, and both maps redraw themselves.

Everything else in this project reuses one stored solution and so takes the physics as fixed.
This is the only tool that re-solves from the electromagnetics down, which is why it is the one to
use when the question is *what if the gel is not what we assumed?*

---

## 1. Opening it

In MATLAB:

```matlab
cd 'C:\Users\Tadhg Stuckey\OneDrive\Desktop\TABLE\BMET3921 - Biomedical Design and Technology\Edward\Heatmodel\matlab'
heat_gui
```

The window takes about **13 seconds** to appear — it runs the first solve before drawing, so it opens
with maps already on screen rather than blank.

It opens on **50 W continuous into room-temperature coolant** — gel and water both at 25 °C, 60 s
exposure. **50 W is far more than this phantom can take**: the peak passes 43 °C at 0.33 s and
100 °C at 4.47 s, and the 60 s figure the window reports — peak 400 °C — is a warning, not a
temperature, because the model has no phase change in it. Drag the *exposure time* slider down to
the few seconds that are physical, or read the time-to-limit off `time_to_T(50, 100)`. Pick **bench
baseline** from the drop-down for the 12 W / 37 °C blood case the rest of the project uses.

`heat_gui` puts the rest of the model on the MATLAB path itself, so that is all you need. Two other
routes to the same place:

```matlab
% from the project root
cd 'C:\...\Edward';  addpath('Heatmodel\matlab');  heat_gui
```

or browse to `Edward\Heatmodel\matlab` in MATLAB's file browser, double-click `heat_gui.m`, and
press **Run** (F5).

Tested on MATLAB R2025b. No toolboxes beyond base MATLAB.

**Giving it to someone else.** `share/Phytagel Heat GUI.zip` in the project root is a standalone
copy — they double-click either the `.mltbx` (installs it into the **Apps** tab) or the `.bat` /
`.command` launcher (installs nothing). Neither route needs a MATLAB command typed, and neither
needs this folder. See the README.

---

## 2. The thirty-second version

1. Pick something from **start from** at the top.
2. Move a slider.
3. Read the sentence above the maps.

That is the whole workflow. There is no "run" step — the maps update on their own about a second
after you let go of a slider.

---

## 3. What you are looking at

### The sentence above the maps

> **50 W for 60 s → peak 400 °C, 239 °C at the probe point (6 mm along, 2 mm deep). That is past
> boiling — this power cannot be used.**

Green below 43 °C, amber above it, red past boiling. Two numbers matter and they are different
things:

- **peak** — the hottest cell anywhere in the 30 mm block. A bound, not something you can measure.
- **probe point** — the temperature at (6 mm along the antenna, 2 mm deep). **This is the number to
  compare against a thermocouple**, and §4 explains why the peak is not.

### The left map — where the power goes

SAR in W/kg on a log scale, four decades from the maximum. This is the microwave deposition, and it
does **not** change with time — only with power and with the gel's permittivity and conductivity.

### The right map — how hot it gets

Temperature in °C at the time shown by the bar underneath. White contour lines are drawn at 43,
50, 60, 80 and 100 °C. Drag the bar to watch the run; it re-reads the stored solution, so scrubbing
is instant and does not re-solve anything.

**The time bar is in seconds and its steps are equal.** The frames sit on a round fixed increment
chosen from the exposure time — 2 s for a 60 s run, 10 s for six minutes, 20 s for ten — the ticks
are labelled with them, and the minor ticks are the frames themselves. Half way along the bar is
half way through the run.

### Markings on both maps

| what | meaning |
|---|---|
| thick cyan / white bar, x = 0 to 12 mm | the monopole |
| thin bar, x = −12 to 0 mm | the braid arm carrying the return current (§4) |
| white dotted line | the catheter footprint (only drawn when the catheter is doing something) |
| **grey wedge, dashed arc** | only appears with the braid arm set to 0 — the blanked feed
  singularity of v1/v2, **not a prediction**. See §4 |

Nothing else in the maps is grey. If you see grey, something has been blanked.

---

## 4. The feed, and where the grey wedge went

**What used to happen.** The field model prescribed the antenna current as a sine that is still at
its maximum where it was cut off at the feed. That step puts a fictitious point charge at the
origin, so the deposition there blew up, so a 2 mm sphere about the feed was blanked grey and
excluded from the reported peak — and the peak then sat on the boundary of that sphere, growing
every time the grid was refined.

**Why that was wrong.** The photographs of the 50 W run show no hole and no spike: the gel at the
feed is heated smoothly, the hottest part of the blob is a couple of millimetres along the element
and a couple down, and the heated region reaches *back past the junction, under the braid*. The
current arriving at the feed does not stop and does not pile up as charge — **it returns**. With no
ground plane and no balun on that bench, the return path is the outside of the bare, immersed coax
braid, so the structure is an asymmetric **dipole** whose second arm is the coax. Truncating the
current at x = 0 was the statement *that arm does not exist*, and the point charge was its price.

**What happens now.** The current continues behind the feed as
`I(x) = I(0)·cos(πx / 2L_sheath)`, going smoothly to zero at the end of the arm. `I` is continuous,
so there is no point charge, nothing is blanked, and the peak the window reports is a prediction. It
also stops depending on the grid — 1.2 % across three grid sizes, against 32 % before.

| control | range | notes |
|---|---|---|
| braid arm (mm) | 0 – 30 | how much bare braid behind the feed carries the return current. **12 mm** (default) makes it the mirror of the element. Shorter piles charge up at the junction; **0** restores the v1/v2 truncated current, grey wedge and all |

**This is a regularisation with a physical argument, not a solved sheath wave.** The arm length is
asserted, not derived, which is exactly why it is a slider — sweep it and see how much of your answer
depends on it. Between 6 mm and 24 mm the peak on the photographed run moves from 167 °C to 117 °C,
so it is not a small dependence.

**The probe point is still the number to compare against a thermocouple.** The peak is now a real
prediction but it sits about 3 mm down beside the feed junction, which is not somewhere a
thermocouple can go. Put the probe under the middle of the element, around x = 6 mm, and read the
probe figure.

Full diagnosis of the original artefact: `Heat_model_writeup.md` §6.2.

---

## 5. The controls

### Start from

Ten starting points. Pick one, then change whatever you like.

Numbers below are **v2.1** — water bath, 12 mm braid arm — at the 0.5 mm grid.

| preset | peak | probe | what it is for |
|---|---|---|---|
| **room-temp water** — 50 W continuous, 60 s | 117 °C | 81 °C | the one the window opens on: gel and water at 25 °C |
| **the photographed run** — 50 W for 6 min at 25 °C | 151 °C | 105 °C | the bench run in `50W CW 6 min/`. Still boils; the gel in the photographs did not |
| **bench baseline** — 12 W for 60 s | 60 °C | 51 °C | 37 °C blood, the case the rest of the project uses |
| **gentle test** — 5 W for 60 s | 46 °C | 43 °C | the one to use for a validation run |
| **high power** — 40 W for 10 s | 73 °C | 54 °C | short and hard |
| **pulsed** — 40 W at 30 % duty, 1 s period | 59 / 60 °C | 51 °C | fast pulsing: almost identical to 12 W continuous |
| **slow pulsing** — 50 W at 50 % duty, 30 s period | 72 / **102** °C | 57 °C | the one worth looking at — see §5 *Drive* |
| **cold gel from the bench** — gel at 25 °C | 55 °C | 45 °C | gel colder than the coolant |
| **write-up baseline** — perfect heat sink | 56 °C | 48 °C | the brief's Dirichlet surface |
| **v1/v2 feed** — truncated current | 154 °C | 128 °C | the old feed model, grey wedge and all, for comparison |

Where two peaks are given they are *at the end of the run* and *at its worst* — see **Drive** below.

**start over** (bottom left) re-applies whichever preset is currently selected, including the feed
and pulsing settings.

### Drive

| control | range | notes |
|---|---|---|
| generator power (W) | 0.1 – 200 | forward power at the generator; the gel receives 75 % of it. With pulsing on, this is the power **during a pulse** |
| exposure time (s) | 0.1 – 600 | the whole run, pulses and gaps together — not the on-time |
| drive | — | continuous, or pulsed |
| pulse period (s) | 0.1 – 60 | one on-time plus one off-time. Greyed out unless the drive is pulsed |
| duty cycle (%) | 1 – 100 | the fraction of each period the power is on. 100 % is the same as continuous |

**Compare pulsed against continuous at equal MEAN power, never at equal generator power.** The
headline gives both: *50 W pulsed at 30 % duty, 1 s period (15 W mean)*. At a fixed mean power a
slower pulse always reaches a higher peak — 12 W mean for 60 s gives 59.5 °C continuous and
59.5 / 59.7 / 61.2 / 76.9 °C at periods of 0.2 / 1 / 5 / 30 s.

**Watch for the second temperature in the headline.** A slow pulse train can finish deep in an
off-period, so the temperature at the end of the run is not the worst the gel saw. The *slow
pulsing* preset reads **72 °C at the end, 102 °C at its worst**, and it is the second number that
decides whether the phantom survives. The window tracks the peak at every pulse-off edge, not just
at the frames the time bar steps through, so it does not miss one between frames.

Pulsing costs time: the solver has to stop at every pulse edge, so a 0.1 s period over a 600 s run
takes several times longer than continuous. The maps are unchanged in shape — pulsing scales the
deposition in time, not in space.

### Gel

| control | range | notes |
|---|---|---|
| gel temperature (C) | 0 – 40 | the temperature the gel **starts** at |
| permittivity eps′ | 1.01 – 120 | real part of the relative permittivity. **Type in it whenever you like** |
| conductivity (S/m) | 0 – 20 | **this is what absorbs the microwave.** Type in it whenever you like |

**Both boxes are always editable.** To enter a measured pair, just type in them — you do not need
to untick anything first, and the tie releases itself when you do.

**The tick box "tie eps′, sigma to the gel temperature"** is on by default, and controls one thing:
whether the pair is *recomputed* from the saline model each time you move the gel temperature. Type
a value into either box and it unticks itself, because typing a measured number is an explicit
override. Tick it again to go back to the saline model at whatever gel temperature is currently set.

The report line under the maps always says which you are looking at — `(Stogryn, 0.9% NaCl at
25.0 degC)` or `(entered by hand)`.

> ⚠ **This is the one genuine trap in the window.** With the box ticked, moving the gel temperature
> changes the *source term* as well as the starting temperature — so it is not a clean experiment.
> Untick it if you want to vary one without the other. (Two of the model's own validation tests were
> wrong on first writing for exactly this reason.)

### Coolant (the blood analogue)

| control | range | notes |
|---|---|---|
| water temperature (C) | 17 – 60 | the temperature the surface is driven towards, and what sets the water properties. The floor is where `water_properties` is tabulated (290 K), not a physical limit — room temperature sits comfortably inside it |
| flow rate (L/min) | 0.1 – 50 | sets how hard the surface is cooled |
| surface cooling | — | see below |

**surface cooling**

- **realistic flow cooling** — the film coefficient *h* is computed from the flow rate. The default.
- **perfect heat sink** — the surface is pinned at the water temperature. This is what the brief
  specified, and it is the **lowest** temperature the gel can reach. **The flow rate does nothing in
  this mode.**
- **none** — no surface cooling at all.

### Feed (the coax braid)

| control | range | notes |
|---|---|---|
| braid arm (mm) | 0 – 30 | how much bare braid behind the feed carries the **return** current. 12 mm (default) makes it the mirror of the element; 0 is the v1/v2 truncated current, which is singular at the feed. **See §4** |

Changing it re-runs the field solve, so it costs about five seconds, not one.

### Catheter

| control | range | notes |
|---|---|---|
| length on the gel (mm) | 0 – 33 | measured back from the antenna tip at x = 12 mm |
| catheter is | — | see below |

**catheter is**

- **cooled by the flow** — it gets its own film coefficient, which is *lower* than the bare surface
  because the wall adds thermal resistance. So a **longer catheter runs slightly hotter**, which is
  the right way round even though it reads backwards at first.
- **ignored** — no effect. **The catheter length does nothing in this mode.**
- **blocking the flow** — the footprint is adiabatic. The pessimistic bracket.

### Advanced

Hidden until you tick **advanced**.

| control | notes |
|---|---|
| grid | 0.5 mm (default, a few seconds) · 1.0 mm (fast, coarse) · 0.25 mm (~90 s, the full baseline grid) |
| SAR scale | log10 or linear |
| T scale | auto, or fixed 20–100 °C for comparing runs side by side. The floor is 20 °C so that a room-temperature run does not clip at the bottom of the colour bar |

Picking the 0.25 mm grid **turns automatic redrawing off by itself** and says so — at 90 s a solve,
a window that redraws on every slider move would look like it had hung. Press **update now** when
you are ready.

---

## 6. How long a change takes

| what you changed | what has to be redone | time |
|---|---|---|
| power, exposure time | conduction only | 1 – 4 s |
| water temperature, flow, catheter, surface model | conduction only | 1 – 4 s |
| permittivity, conductivity, gel temperature | **the field solve as well** | ~5 s |
| grid → 0.25 mm | everything, on 1.8 M nodes | ~90 s |

The electromagnetic solution is cached and reused whenever the gel and the grid are unchanged, which
is why the drive is so much cheaper to move than the gel.

If you would rather it never redrew on its own, untick **update as I change things**. The status
line then turns orange and says the maps are out of date until you press **update now**.

---

## 7. Reading the flag lines

Under the maps, lines beginning `!` are warnings about the run you just did. The ones you will meet:

| flag | what it means |
|---|---|
| *peak N °C is above boiling* | there is no phase change in this model. Past 100 °C the number is a warning, not a temperature. |
| *above 60 °C the constant-property assumption costs 10–15 %* | and it **under**-predicts — the real gel absorbs more as it warms |
| *above ~50 °C low-acyl gellan may be softening* | the phantom may stop being a phantom before the model stops being valid |
| *Dirichlet surface: the flow rate has NO effect* | you are on perfect heat sink; the flow slider is inert |
| *transparent catheter: the catheter length has NO effect* | the catheter is set to "ignored"; the length slider is inert |
| *the gel starts N K from the coolant* | the first seconds are a relaxation as well as a heating transient |
| *braid arm off: the current is truncated at the feed* | you set the braid arm to 0. A 2 mm sphere is blanked and the peak is the value on its edge, not a prediction (§4) |
| *pulsed: N W during a pulse, M W mean* | a reminder to compare at equal mean power |
| *the pulse period is far shorter than the local thermal time constant* | this run is indistinguishable from continuous drive at the mean power |
| *the braid arm reaches x = −N mm but the block stops at −9 mm* | the arm runs out of the modelled block, whose side is pinned at the coolant temperature. Heating from the last few mm of it is outside the model |

---

## 8. Three things you might actually want to do

**"How much power can I use for 60 s?"**
Start from *bench baseline*, then drag the power slider down until the sentence turns green (below
43 °C) or until the peak reads whatever limit you are working to. For a precise answer without
hunting, use `P_for_T(100, 60)` in the Command Window instead — it inverts the stored model directly.

**"My gel measured differently."**
Type the measured permittivity and conductivity straight into their boxes — the tie to the gel
temperature releases itself — and compare the maps against the default pair (74.29 and 2.723 S/m).
The report line will read `(entered by hand)` once it is your numbers driving the solve. This is the
case nothing else in the project can answer.

**"I want to compare against a thermocouple."**
Run *gentle test* (5 W, 60 s) — deliberately chosen to stay in the range where the model is testing
itself rather than the gel's softening behaviour. Read the **probe point** number, not the peak.
Then run it again on *perfect heat sink*: the measurement should fall between the two surface
models. If it falls outside that bracket, the flow geometry is not what the model assumes, and that
is the first thing to check.

For a full bench run with traces, a cool-down curve and a CSV to plot against your measurement, use
`heat_scenario` instead — see `Heat_model_writeup.md` §11.4.

---

## 9. How much to trust it

The GUI runs on a coarser grid and a smaller FFT window than the full model, and only as far in time
as the exposure lasts. That was measured against the stored quarter-hour solves, not assumed:

| check | result |
|---|---|
| vs the stored perfect-heat-sink baseline, probe point at 60 s | **1.1 %** of the temperature rise |
| vs the stored convective baseline, probe point at 60 s | **3.8 %** of the temperature rise |
| grid convergence over the last refinement | 0.5 % |
| scaling with power | exact to machine precision |

Run `verify_gui_model` to reproduce all twelve checks (~4 min); the stored output is in
`results/gui_verification_log.txt`.

**Use the full model instead when** you need the 600 s behaviour, the steady state, thermal dose
maps, or a number you are going to put in a report — `run_heat_model`, then `T_gel`.

---

## 10. Saving what you see

**save PNG** writes the whole window — controls included, so the image records the parameters that
produced the maps — to `Heatmodel/figures/gui_<date>_<time>.png`.

---

## 11. If something goes wrong

| symptom | fix |
|---|---|
| `Undefined function 'heat_gui'` | you are in the wrong folder. `cd` to `Edward\Heatmodel\matlab` first. |
| `Unrecognized function 'saline_stogryn'` | the field model's `matlab\` folder is missing from the path — `heat_gui` normally adds it, so check the folder layout is intact |
| window is bigger than the screen | it sizes itself to the screen on open; reopen it rather than resizing |
| a warning about Stogryn being an extrapolation | gel temperature is outside 0–40 °C, where the saline model is fitted. Untick the box and enter eps′ and σ by hand. |
| the maps never change | check **update as I change things** is ticked, or press **update now** |
| status line says FAILED | the message says why; the commonest is a coolant temperature outside 17–96 °C, where the water property table ends |

---

## Related

| file | what |
|---|---|
| `Heat_model_writeup.md` §11.5 | the GUI's place in the model, and what the fast solve costs |
| `Heat_model_writeup.md` §6.2 | the feed singularity — the grey wedge |
| `matlab/gui_model.m` | the same solve without the window, for scripting a sweep |
| `matlab/verify_gui_model.m` | the twelve validation checks |
| `matlab/heat_scenario.m` | one bench run with traces and a CSV |
