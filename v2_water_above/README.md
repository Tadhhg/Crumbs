# v2 — the antenna is in the bath, not in air

Self-contained copy of the model. `../` (v1) is untouched.

---

# v2.1 — the feed, and pulsing

Two changes on top of v2. **Neither touches the bath**, so every v2 number about where the power
goes is unchanged.

## 1. The current returns along the braid, so the grey wedge is gone

**The symptom.** Both maps had a grey semicircle punched out of them at the feed, with a dashed arc
round it, and the reported peak temperature sat on its edge. The photographs of the 50 W run show
nothing of the kind: the gel there is heated smoothly, the hottest part of the blob is a couple of
millimetres along the element and a couple of millimetres down, and the heated region reaches
**back past the junction, under the braid**. No hole, no spike.

**The cause.** v1 and v2 prescribed the current on `0 ≤ x ≤ L` and then stopped it. Under the
quarter-wave hypothesis the current is at its *maximum* where it stops, so charge continuity
`dI/dx + jωρ = 0` puts a point charge `Q = I_m/ω = 2.6e-10 C` at the feed. Its quasi-static field
0.5 mm away in ε′ = 74.3 is 1.3e5 V/m, which is exactly what the solver returned. The deposition —
and so the temperature — was singular there, and blanking a 2 mm sphere was the only defence.

**Why that is an artefact.** The current arriving at the junction does not stop and does not pile up
as charge: **it returns**. There is no ground plane and no balun on this bench, so the return path is
the outside of the bare, immersed coax braid. The structure is not a monopole over a ground plane —
it is an **asymmetric dipole whose second arm is the braid**. Truncating the current at `x = 0` is
exactly the statement *the second arm does not exist*, and the point charge is what that statement
costs. v2's own notes already listed "current on the coax braid, which is bare and immersed — not
modelled" as an unfixed source of error. This is that item.

**What changed.** `model_config` gains `L_sheath` (default 12 mm = L), and the current spectrum
carries a second segment behind the feed:

$$I(x) = \sin\big(k_e (L-x)\big),\ 0\le x\le L; \qquad
  I(x) = I(0)\,\cos\!\Big(\frac{\pi x}{2 L_\text{sheath}}\Big),\ -L_\text{sheath}\le x < 0$$

`I` is continuous at the feed, so `dI/dx` has no delta and there is no point charge. With
`L_sheath = L` the arm is the mirror of the element and `dI/dx` vanishes at `x = 0` from *both*
sides — the charge density is zero at the feed and largest at the two far ends, which is what a
dipole does. The cosine form is used rather than `sin(k_e(L_sheath+x))` because it stays bounded for
every `L_sheath`; the two are identical when `L_sheath = λ_e/4`.

**This is a regularisation with a physical argument, not a solved sheath wave.** The GUI exposes
`L_sheath` ("Feed → braid arm") so it can be swept, and `L_sheath = 0` reproduces v1/v2 **bit for
bit** — the taper segment is then empty and the element quadrature is untouched.

**What it buys.** `r_feed_excl` goes to 0, so nothing is blanked and the reported peak is a genuine
prediction rather than the value on the edge of an arbitrary sphere. On the photographed run
(50 W, 6 min, everything at 25 °C, tap-water bath) the peak falls from **185 °C to 151 °C** and lands
at x = 0.5 mm, 3.0 mm deep — the photographs put the coloured core at roughly x = 2–3 mm, 2–3 mm
deep.

| | truncated (v1/v2) | braid arm (v2.1) |
|---|---|---|
| jump in `I` at the feed, from the spectral tail | **1.0014** (the whole current) | **0.0136** |
| peak deposition, y = 0 plane | 5.3e8 W/m³ | 4.1e6 W/m³ |
| peak temperature at dh = 1 / 0.5 / 0.25 mm | 159 / 117 / 117 °C | 86.6 / 85.6 / 86.1 °C |
| spread across those grids | **32 %** | **1.2 %** |
| blanked | a 2 mm sphere | nothing |

The spectral-tail row is the decisive one. Integrating by parts, a jump `ΔI` at the feed makes
`k_x·Ĩ(k_x) → ΔI` as `k_x → ∞`, and a continuous `I` makes it → 0 (Riemann–Lebesgue). It measures
the point charge directly, with no grid and no threshold to argue about.

> **A first attempt got this wrong and is worth recording.** The obvious fix is to taper the current
> to zero over the *length of the solder joint*, ~2 mm, as though the charge simply lived on the
> junction rather than at a point. That makes it **worse**: squeezing the whole terminal current
> into 2 mm gives `max|dI/dx| = I(0)·π/(2·2 mm)`, **six times** the element's own peak at the tip, so
> the deposition at the feed goes *up* by an order of magnitude and the peak temperature stops
> converging in the other direction (369 → 227 → 191 °C as the grid refines). The junction length is
> not the right scale — the length of the **return path** is, and it is comparable to the element,
> not to the solder.

> **The photographs cannot confirm this.** `run_v2_comparison` asks what temperature each version
> would have to be calling the edge of its blob for that blob to be the 9.5 × 6.3 mm the photograph
> shows, once from the length and once from the depth. Two different answers would mean the
> proportions are wrong. Every version, v1 included, answers within 3 K. The photograph constrains
> the **magnitude** and says nothing about the feed. The case for the braid arm rests on the table
> above, not on the pictures.

## 2. Pulsing

`heat_scenario` could pulse; the GUI could not, and its header said so. It can now: **drive**
(continuous / pulsed), **pulse period** (0.1–60 s) and **duty cycle** (1–100 %).

The gate is applied inside the time loop rather than by superposing step responses, because this
solver is transient already and its initial condition need not be equilibrium. The cost is that the
time step has to land on every pulse edge — otherwise a step straddling an edge silently delivers
the wrong energy and the duty cycle becomes whatever the step size made it. `heat_solve` clips each
step to the next edge, so the delivered energy is exact for any period: **1.9e-14** relative, over
periods from 0.2 s to 30 s.

Two things the window now says that it could not before:

- **Generator power is the power during a pulse**, and the mean is `duty ×` that. Both are quoted.
  Compare pulsed against continuous at equal **mean** power, never at equal generator power.
- **`T_end` is not the worst of the run.** A 30 s period at 50 % duty finishes 15 s into an
  off-period, and the gel has been cooling for all of it — 72 °C at the end against **102 °C at its
  worst**. `heat_solve` tracks the peak at every pulse-off edge, not just at the output times, and
  the headline reports both when they differ.

At a fixed 12 W mean and 60 s, the peak rises monotonically with the pulse period — 59.5 °C
continuous, 59.5 / 59.7 / 61.2 / 76.9 °C at periods of 0.2 / 1 / 5 / 30 s. Slow pulsing at the same
mean power is strictly worse, which is the write-up's §8 conclusion and is now reproducible in the
window.

## 3. The time bar

It was an index into a **log-spaced** list of output times. That is the right grid for a document
figure — it resolves the first millisecond of a run that then lasts ten minutes — and the wrong one
for a scrub bar, because the bar moves in equal steps of *index*: a quarter of its travel covered
the first second of a sixty-second run and the remaining three quarters covered the other
fifty-nine.

The output times are now on a **fixed increment**, chosen from the 1-2-5 ladder to give about 30
frames (2 s for a 60 s run, 10 s for 360 s, 20 s for 600 s), and the bar is linear in **seconds**
with ticks on round numbers and one minor tick per frame. Half way along the bar is half way through
the run.

Anything that genuinely wants the early transient can still ask for it: `gui_model` takes `t_out`
explicitly, and `run_v2_comparison` passes a log-spaced grid so its "time to cross 100 °C" table
still answers in hundredths of a second.

## What v2.1 touched

| file | change |
|---|---|
| `matlab/model_config.m` | **new** `L_sheath` (12 mm), the braid arm behind the feed |
| `matlab/layered_field.m` | the arm in the current spectrum; two Simpson segments meeting at the feed |
| `Heatmodel/matlab/sar_volume.m` | the same, plus `debug_current` so `verify_v2` can read the spectral tail |
| `Heatmodel/matlab/feed_regularise.m` | `r_feed_excl = 0` now disables the mask completely (strict inequality) |
| `Heatmodel/matlab/thermal_config.m` | `r_feed_excl` documented as superseded by the arm |
| `Heatmodel/matlab/heat_solve.m` | pulsing: `period`, `duty`, edge-clipped stepping, `E_in`, `E_loss`, `Th_max_run` |
| `Heatmodel/matlab/gui_model.m` | `L_sheath_mm`, `pulse_period`, `duty`, `t_out`; fixed-increment output times; `T_worst` |
| `Heatmodel/matlab/heat_gui.m` | the Drive and Feed controls, the time bar in seconds, the maps cropped to the data |
| `verify_v2.m` | tests 6 and 7 |
| `run_v2_comparison.m` | all five versions, and the shape test against the photographs |
| `share/.../build_mltbx.m` | **new** — repackages `code/` into the `.mltbx` |

The energy audit changed with it. `P_loss` was integrated by trapezoid **between output times**,
which is fine for a smooth CW run and meaningless under pulsing, where `P_loss` swings with the
drive and the output times sample it aliased. It is now accumulated at every **time step**, using
the same instant of `theta` the forward Euler step itself uses. The audit closes to **5e-15** for
both CW and pulsed runs, against 3e-4 and 5e-2 before.

## Verifying it

`verify_v2` gained tests 6 and 7 — the spectral tail, the grid-convergence of the peak, the budget
with the arm on, and the three energy checks on pulsing. All checks pass in about two minutes.

## Rebuilding the package

`build_mltbx.m` (beside `START HERE.txt`) repackages `code/` into the `.mltbx`. **Run it after
changing the model** — the `.mltbx` is a copy, so without a rebuild the Apps-tab route ships the old
model while the `.bat` route ships the new one.

---

# v2 — what it changed (unchanged by v2.1)

The change v2 itself made: **the medium above the gel**.

## Why

v1 solved a three-layer stack whose top layer was **air** — `layered_field.m` header,
`layer 0 : z > 0 air (source: x-directed wire at z = h)`, with `epsr_air = 1 - 0j`. Against the
gel's ε′ = 74.3 that is an enormous contrast, and it is what drove **99.6 % of the delivered power
into the gel**.

The photographs of the 50 W / 6 min run (`../50W CW 6 min/`) show the antenna **fully immersed**.
Layer 0 is a liquid whose permittivity is close to the gel's own, so the interface nearly
disappears electrically and a large share of the power goes **up into the bath** instead.

## What changed

| file | change |
|---|---|
| `matlab/bath_medium.m` | **new** — bath permittivity, same Stogryn model as the gel |
| `matlab/bath_spectral_E.m` | **new** — spectral E in the bath, above and below the wire |
| `matlab/model_config.m` | `upper` (`'water'`/`'air'`), `bath_NaCl_pct`, `T_bath_C`, `Zw`, `dzw`; resolves `epsr_upper` after overrides |
| `matlab/layered_field.m` | bath absorption in the power budget |
| `Heatmodel/matlab/sar_volume.m` | bath absorption in the power budget; interface self-check |
| `Heatmodel/matlab/gui_model.m` | passes the bath through; reports `frac_bath`, `Pwat`, `epsr_upper` |

### The part that is easy to get wrong

Changing `epsr(1)` alone is **not** enough, and would have made the answer worse rather than better.
v1 normalises by `scale = sqrt(Pdel/Ptot)` with `Ptot = Pgel + Pup + Pdn`. With a *lossless* upper
half-space that is complete: every plane above the source carries the same real flux. With a
*lossy* bath it is not — the bath absorbs, so the flux depends on which plane you measure it on,
and the absorbed power is missing from `Ptot`. Dividing by too small a total would have **inflated**
the field in the gel.

So v2 integrates the dissipation over a bath column `0 ≤ z ≤ Zw` exactly as `Pgel` is integrated
over the slab, and takes `Pup` as the **residual** flux leaving the top of that column. Nothing is
double counted, and the budget closes for any `Zw` — verified below. The bath integral is done in
the spectral domain (Parseval), so it costs no extra FFTs.

## Verification — `verify_v2`

| # | check | result |
|---|---|---|
| 1 | tangential **E** continuous across z = 0 | rel err **1.6e-17** |
| 2 | `upper='air'` reproduces v1 | `Pgel`, `Pup`, `Pdn` and the whole deposition map **identical to machine precision**; bath absorbs exactly 0 |
| 3 | cut height `Zw` does not matter | gel share 0.3917 at Zw = 20/40/80 mm, drift **5e-8**, while bath and escaping flux trade off |
| 4 | budget closes, fractions sum to 1 | exact |
| 5 | uniform medium, wire **on** z = 0, up/down symmetric | **1.000000** |
| 6a | `L_sheath = 0` reproduces the truncated current | max abs dq **exactly 0** |
| 6b | jump in `I` at the feed, from the spectral tail | truncated **1.0014**, braid arm **0.0136** |
| 6c | peak temperature vs grid, dh = 1/0.5/0.25 mm | truncated spread **32 %**, braid arm **1.2 %** |
| 6d | budget still closes with the arm on | rel **1.1e-16** |
| 7a | duty = 1 reproduces continuous drive | max abs dT **exactly 0** |
| 7b | energy delivered = duty x P x t, periods 0.2-30 s | rel **1.9e-14**, peak monotone in the period |
| 7c | continuous at the same mean power is the floor | 59.5 °C against 76.9 °C at a 30 s period |

Test 1 was written before the sign was known and **failed** at first (~1 %), which is what caught it.
Continuity forces `Z0tm*(1-r01b) == Z1tm*(1+r01b)`; substituting `r01b` makes both sides
`2*Z0tm*Z1tm/(Z0tm+Z1tm)` identically, so `S_TM = +1`. Test 5's first version was also wrong — it
left the wire at h = 0.5 mm and expected symmetry, but the wire's near field then lies inside the
upper region (ratio came out 16×). Moving the source onto the interface is what makes the two
sides comparable. **The expectation was wrong, not the code.**

## Result — the photographed run (50 W CW, 360 s, everything at 25 °C)

| | gel % | bath % | P into gel | peak | probe | peak is |
|---|---|---|---|---|---|---|
| v1, air above | 99.6 | 0.0 | 37.4 W | 478 °C | 316 °C | on a blanked sphere |
| v2, deionised/tap bath | 46.4 | 51.3 | 17.4 W | 185 °C | 162 °C | on a blanked sphere |
| v2, 0.9 % saline bath | 39.2 | 60.7 | 14.7 W | 171 °C | 148 °C | on a blanked sphere |
| **v2.1**, deionised/tap bath | 47.9 | 49.3 | 18.0 W | **151 °C** | 105 °C | **a prediction** |
| **v2.1**, 0.9 % saline bath | 43.8 | 56.1 | 16.4 W | **147 °C** | 103 °C | **a prediction** |

The gel share moves by about a point and a half because the braid arm is a second radiator sitting
on the interface, not because anything in the bath calculation changed.

The bath salinity is **not known**, so both ends are run. It matters less than expected: pure water
is already lossy at 2.45 GHz through dipolar relaxation (σ_eff = 1.23 S/m against saline's 2.72),
so the bracket is narrow and the conclusion is robust to not knowing it.

The **shape** barely moves — d50 1.21 → 1.11 mm, d90 7.10 → 6.08 mm, peak still pinned to the feed
exclusion boundary. The bath is essentially a **scalar** on the deposition, which is consistent with
the photographs, where the heated zone's position and proportions looked right and only its
magnitude did not.

## What v2 does NOT fix

**It closes about 2.8× of the gap and no more.** 171 °C still boils, and the photographs show an
intact gel. The remaining factor of 3–6 must come from the causes v2 does not touch:

1. **η = 0.75 is nameplate, not measured** — `Pdel = Pfwd - Pref` = 200 − 50 W, the generator's
   maximum ratings. An unmatched bare wire immersed in water will reflect far more than 25 %.
   **Measure forward and reflected power.** This is still the highest-value measurement.
2. ~~**Current on the coax braid**, which is bare and immersed. Not modelled.~~
   **Modelled in v2.1** -- see the top of this document. It removes the feed singularity and about
   20 % of the remaining peak, and it is a regularisation with a physical argument rather than a
   solved sheath wave: `L_sheath` is asserted, not derived. Sweep it.
3. **`ke_case = 'quarterwave'`** assumes the element is resonant, which maximises the current.
   `model_config` itself calls this the single largest modelling uncertainty.
4. **Standoff** — v1 and v2 both use h = 0.5 mm with a 0.45 mm conductor; the photos show a clear
   gap and a thicker wire.

## Running it

```matlab
cd v2_water_above
verify_v2                          % all the checks above, ~2 min
run_v2_comparison                  % the table above, all five versions
heat_gui                           % the GUI: water bath and braid arm by default

% falling back, for comparison
gui_model('upper','air')                   % v1's medium
gui_model('L_sheath_mm',0)                 % v1/v2's truncated feed, wedge and all
gui_model('upper','air','L_sheath_mm',0)   % v1 exactly

% pulsing
gui_model('Pgen',40,'pulse_period',1,'duty',0.3)   % 40 W peak, 12 W mean
```

In the window the same three fall-backs are **Bath -> above the gel -> air**, **Feed -> braid arm ->
0**, and the **v1/v2 feed** preset.
