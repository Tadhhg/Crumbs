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

Temperature in °C at the time shown by the slider underneath. White contour lines are drawn at 43,
50, 60, 80 and 100 °C. Drag the time slider to watch the run; it re-reads the stored solution, so
scrubbing is instant and does not re-solve anything.

### Markings on both maps

| what | meaning |
|---|---|
| cyan / white bar along the top | the 12 mm monopole, from x = 0 to x = 12 mm |
| white dotted line | the catheter footprint (only drawn when the catheter is doing something) |
| **grey wedge, dashed arc** | the feed singularity — **not a prediction**, see §4 |

---

## 4. The grey wedge, and why the peak sits on its edge

The field model prescribes the antenna current as a sine that is still at maximum where it is cut
off at the feed. That step puts a fictitious point charge at the origin, so the deposition there
blows up. It is an artefact of the prescribed current, not physics — at the tip the current goes
smoothly to zero and there is no such spike.

So a 2 mm sphere about the feed is **blanked grey and excluded from the reported peak**. The
deposition is still rising towards the feed just outside that sphere, which is why the peak is
almost always reported at around **x = 0.5 mm, 2 mm deep** — sitting on the boundary of the excluded
region. Refining the grid keeps finding a larger value there. That is expected and it is why the
probe point is quoted beside it.

**Do not put a bench probe within 2 mm of the feed.** Put it under the middle of the element, around
x = 6 mm. A real feed junction *is* a hot spot, but this model does not predict its size correctly.

Full diagnosis: `Heat_model_writeup.md` §6.2.

---

## 5. The controls

### Start from

Six starting points. Pick one, then change whatever you like.

There is **no pulsing in this window** — the drive is continuous by definition, and "exposure time"
is how long it stays on. For a duty cycle use `heat_scenario(P, period, t, 'duty', d)` instead.

| preset | what it gives |
|---|---|
| **room-temp water** — 50 W continuous, 60 s | the one this window opens on: gel and water at 25 °C. peak 400 °C, probe 239 °C — past boiling; 100 °C is reached at 4.5 s |
| **bench baseline** — 12 W for 60 s | peak 131 °C, probe 89 °C |
| **gentle test** — 5 W for 60 s | peak 76 °C, probe 59 °C — the one to use for a validation run |
| **high power** — 40 W for 10 s | peak 200 °C, probe 103 °C |
| **cold gel from the bench** — gel at 25 °C | peak 121 °C, probe 83 °C |
| **write-up baseline** — perfect heat sink | peak 98 °C, probe 78 °C |

The last one reproduces the main document's §9.5 to the last degree: **97.6 °C** here against the
97.7 °C from the quarter-hour solve. It is the check that the fast model is telling the truth.

**start over** (bottom left) re-applies whichever preset is currently selected.

### Drive

| control | range | notes |
|---|---|---|
| generator power (W) | 0.1 – 200 | forward power at the generator; the gel receives 75 % of it |
| exposure time (s) | 0.1 – 600 | how long the microwave is on |

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
