# Gel model V4 — the photographed setup by default, every constant editable

Date: 2026-09-16 · Builds on `../Gel_model_V3` (field model unchanged; heating chain re-parameterised).
Why the V3 pattern did not match the photographs, and the measurements behind these defaults:
`../Gel_model_V3/gel_heating_notes.md` §5–6.

```matlab
cd Gel_model_V4
gel_heating_gui          % two tabs: Heating | Constants & boundaries   (first call solves, ~30 s)
run_gel_heating          % the three-panel figure for the photographed run -> gel_heating.png
R = gel_heating_solve(); % the solution struct (cached in solutions/)
```

## What changed from V3

| | V3 | V4 |
|---|---|---|
| constants | scattered through the scripts | **one struct, `model_defaults.m`**, every field with its source; edited on the GUI's second tab |
| antenna | 27 mm bare wire, 27 mm braid arm, $k_e=k_\text{gel}$, axis 0.95 mm up | **12 mm PTFE-covered wire (RG-178 assumed), 4 mm arm, $k_e$ from King's insulated-antenna formula (`king_insulated_ke.m`, 127−30j m⁻¹), axis 1.1 mm up** — as measured on the photographs |
| initial state | gel in equilibrium with 37 °C water | gel in equilibrium with the bath at `T_bath0` (**25 °C**, the photographed day); `T = T_ss + (T_w − T_bath0)·U + (P_set/P_ref)·Θ` |
| block | 30 mm block starting 1.5 mm behind the feed | edges where the photographs put them: **x = −6.5 … 23.5 mm**, 30 mm wide, 30 mm deep |
| top face | Robin, $h$ from the laminar plate | same, or a number, or **`Inf` = pinned at the water temperature** |
| grid | fixed 0.5 mm | **`dh` editable**; node count, time step and an estimated solve time shown before you commit |
| cache | one file per gap | one file per distinct set of constants (`solutions/gel_heating_<md5>.mat`) |
| lossy medium starts at | the conductor radius | the **insulation** radius `b_ins` (set `b_ins = a` for a bare wire) |

`E_field_2d.m` is a verbatim copy of V3's (homogeneous lossy gel; see its notes for the equations). The
water above the gel is still not in the field model — that is the next physics step (V3 notes §7).

## The Constants tab

Every field of `model_defaults` in five groups — antenna, gel, water above the gel, tank walls, block and grid —
as numbers you type (lengths in mm, flow in L/min, frequency in GHz). Three fields accept **empty** to mean
"use the model's own value": `ke` (→ King's formula), `h_top` (→ flat-plate correlation), `h_wall` (→ air and
acrylic in series). `h_top = Inf` pins the surface. Complex numbers are typed as `127-30j`.

The right-hand panel updates as you type: $k_e$ and the element length in wavelengths, the gel's wavelength and
power depth, $h_\text{top}$, $h_\text{wall}$ and their Biot numbers, the water gap, the grid ($n_x n_y n_z$),
the time step, the number of steps and an **estimated solve time**. It also says whether the constants shown
are those of the solution currently on the Heating tab.

**Apply & re-solve** runs `gel_heating_solve` for the new constants (progress dialog), caches the result and
switches to the Heating tab. **Photo defaults** reloads `model_defaults`; **Current solution** reloads the
constants of whatever is displayed. Invalid entries (a non-number, insulation thinner than the conductor,
a block that does not contain the element, …) are reported and nothing is solved.

### Grid fineness

`dh` is the only spatial resolution in the chain (the field is analytic). Cost scales as $dh^{-5}$ — the
node count as $dh^{-3}$ and the stable explicit time step as $dh^{2}$:

| dh | nodes (30 mm block) | steps to 720 s | solve |
|---|---|---|---|
| 1.0 mm | 29 791 | 1 440 | ~1 s |
| 0.5 mm (default) | 226 981 | 3 960 | ~30 s |
| 0.25 mm | 1.77 M | ~15 000 | ~15 min |

The source under the wire goes as $1/\rho^2$, so the power captured by the block is grid-sensitive: 7.1 W at
1 mm vs 6.6 W at 0.5 mm (of 20 W). Explore at 1 mm, quote at 0.5 mm or finer, and check the peak moves by
less than you care about between two spacings (V3 notes §3 lists this as an outstanding verification).

## Result for the photographed run (defaults, 50 W set = 10 W delivered)

| t | peak | at (x, z) | photographs |
|---|---|---|---|
| 29 s | 38.6 °C | (−2.0, −1.0) mm | spot at x ≈ 2–4 mm, ~2.5–3 mm deep |
| 92 s | 44.5 °C | (1.0, −2.5) mm | hottest x ≈ 3.1–3.6 mm, 3.3–3.4 mm deep |
| 332 s | 50.7 °C | (1.5, −3.5) mm | x ≈ 3.3 mm, ~3.5 mm deep; no boiling |

33 % of the delivered power lands in the gel (V3: 16 %). Depth is right; the peak sits 1–2 mm too close to
the feed, and at 29 s it is under the end of the 4 mm braid arm (the cosine arm's charge sits at its end) —
the arm remains an asserted regularisation, not a solved sheath current. The surface node under the wire
still rises 6–7 K where the pictures show it below the pigment's start temperature; `h_top = Inf` shows the
other bracket.

## Still to confirm (each is one field on the tab)

- the cable: conductor radius `a`, PTFE radius `b_ins` — or a bare wire (`b_ins = a`) — and the element length
- the pigment's start temperature and bandwidth (turns the colour bands into isotherms)
- bath temperature and flow cross-section on the day (`T_bath0`, `A_flow`)
- forward/reflected power (`eta`)

## Files

| file | role |
|---|---|
| `model_defaults.m` | the constants, with sources |
| `model_derived.m` | everything that follows from them without a solve (wavenumbers, $h$, grid, dt, cost, cache key) |
| `king_insulated_ke.m` | King, Trembly & Strohbehn's insulated-antenna wavenumber |
| `E_field_2d.m` | field of the wire in the homogeneous lossy gel (V3, unchanged) |
| `gel_heating_solve.m` | 3-D steady state + unit responses for a struct `P`; cached |
| `gel_heating_eval.m` | any power / water temperature / timing from the stored responses; sparing metrics |
| `gel_heating_optimise.m` | least power that makes a ≥ 50 °C lesion under a ≤ 45 °C layer |
| `gel_heating_gui.m` | the two-tab GUI |
| `run_gel_heating.m` | the summary figure |
| `solutions/` | cached solutions, one per distinct `P` (safe to delete) |
