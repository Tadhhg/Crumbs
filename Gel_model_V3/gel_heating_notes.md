# Gel heating — 3-D conduction driven by the E-field model, and why it does not look like the photographs

Date: 2026-09-16 · Code: `gel_heating_solve.m`, `gel_heating_eval.m`, `gel_heating_optimise.m`, `gel_heating_gui.m`,
`run_gel_heating.m` · Tests: `photo_comparison.m` · Source: `E_field_2d.m` (`E_field_2d_notes.md`) ·
Photographs: `../50W CW 6 min/`, `mirrored.jpg`, `photo_series.jpg`

Part A (§1–4) documents the heating chain as it stands. Part B (§5–7) compares it with the 50 W / 6 min
photographs, reports the tests, and ranks the reasons for the mismatch.

---

## Part A — the model

## 0. What you get

```matlab
R = gel_heating_solve();                       % one-off 3-D solve (~1-2 min), cached in gel_heating_solution_gap0.5mm.mat
[T, M] = gel_heating_eval(R, 100, 37, 0, 360); % P_set [W], T_water [degC], pre-cool [s], t [s]  -> plane y = 0, metrics
gel_heating_gui                                % sliders: power, water T, pre-cool, exposure, time; Optimise button
best = gel_heating_optimise(R);                % least power that makes a >= 50 degC lesion under a <= 45 degC layer
run_gel_heating                                % the three-panel figure gel_heating.png

% the antenna of the photographs (added 2026-09-16, s6): 12 mm element, 4 mm braid arm, insulated-wire ke
R12 = gel_heating_solve('L',12e-3, 'L_sheath',4e-3, 'ke',172-26j, 'gap',0.65e-3);   % ~10 s, own cache file
% the photographed day had gel AND bath at 25 degC: then T = 25 + (P/P_ref) Theta, i.e.
f = find(R12.t == 92);  T92 = 25 + (50/100)*double(R12.Theta(:, :, f));               % 50 W set, t = 92 s
% (gel_heating_eval(R12, 50, 25, ...) is a different scenario: gel equilibrated with 37 degC water,
%  water stepped to 25 degC at t = 0.)
```

For the record, that solve puts 5.4 W of the 20 W reference into the gel (the 27 mm default: 3.2 W). With everything
at 25 °C and 50 W set, its peak is 37.9 / 45.0 / 52.4 °C at $z=-1.5/-2.5/-3.5$ mm, $x=-1$…1.5 mm, at 30 / 92 / 332 s,
with the surface node in that column at 30.7–32.2 °C — within a millimetre of the photographs in depth (§5.2), a
millimetre or two too close to the feed in $x$, and the surface still warmer than the pictures suggest (§6.3).

`gel_heating_solve` stores **unit responses** on the antenna plane $y=0$; everything else is instant because the
problem is linear (§1, eq. 4). `'gap'`, `'L'`, `'L_sheath'`, `'ke'` each get their own cache file.

## 1. The chain, five equations

**(1) Source.** From `E_field_2d`: $q=\tfrac12\sigma_\text{eff}|\tilde{\mathbf E}|^2=\sigma_\text{eff}E_\text{rms}^2$ [W/m³],
$\sigma_\text{eff}=\omega\varepsilon_0\varepsilon''=1.84$ S/m. The field of the thin wire is axisymmetric about the wire
axis, so the plane $y=0$ sampled along $\rho=h-z$ gives $q(x,\rho)$ for the whole block: `interp2` on
$\log\rho$, $\rho\in[0.6,70]$ mm. $P_\text{gel}=\sum q\,dh^3$ over the block is reported; the rest of $P_\text{del}$
is dissipated outside the block (the near-field shell round the wire, i.e. in the water) — **16 % in the gel at
0.5 mm gap, 7 % at 3 mm** (§4). Everything downstream scales with $P_\text{gel}$.

**(2) Heat equation** in the $30\times30\times30$ mm block, Phytagel ≈ water: $\rho c\,\partial T/\partial t=k\nabla^2T+q$,
$k=0.60$, $\rho c=4.18\times10^6$, $\alpha=1.435\times10^{-7}$ m²/s. No convection inside a 1 % gel, no perfusion,
no phase change (checked: nothing reaches 100 °C at the powers used).

**(3) Boundaries — every face Robin**, $-k\,\partial T/\partial n=h\,(T-T_\infty)$, as in `gel_heat_2d.m`:

| face | $h$ [W/(m²K)] | $T_\infty$ | where from |
|---|---|---|---|
| top (water) | 1460 | $T_w$ (37 °C) | laminar flat plate $\overline{Nu}_L=0.664\,Re_L^{1/2}Pr^{1/3}$, Incropera & DeWitt ch. 7; $U=Q/(30\times30$ mm$)=0.093$ m/s, $Re_L=3980$. Bi $=hL/k=73$: the face sits within ~0.5 K of the water — *when the source is not right under it* (§6.3) |
| sides, bottom | 7.9 | 25 °C room | free convection + radiation (~10) in series with 5 mm acrylic ($k=0.19$) |

FTCS in 3-D with a ghost-node layer for the Robin faces; $dt=0.9\,dh^2/[\alpha(6+g_\text{top}+2g_\text{wall})]$,
$g=2\,dh\,h/k$; $dh=0.5$ mm; frames every 2 s to 720 s.

**(4) Linear superposition.** With the gel initially in equilibrium with 37 °C water,

$$T(\mathbf r,t)=T_{ss,37}(\mathbf r)+(T_w-37)\,U(\mathbf r,t_\text{pre}+t)+\frac{P_\text{set}}{P_\text{ref}}\,\Theta(\mathbf r,t)$$

$T_{ss,37}$: steady state without power (sparse direct/PCG solve, gel 30.5–37 °C because the walls leak to the room);
$U$: response to a +1 K step of the water; $\Theta$: response to the source at $P_\text{ref}=100$ W ($\eta=0.2$, 20 W
delivered). Superposition is exact because (2)–(3) are linear and $q\propto P_\text{del}$.

**(5) Catheter**, lumped, cylinder in cross-flow (Churchill–Bernstein, Incropera ch. 7):
$h_cA_c(T_c-T_w)=\tfrac{P}{P_\text{ref}}P_\text{film}+G_\text{gap}(T_s-T_w)$, $G_\text{gap}=k_wDL/g$, $P_\text{film}$ the
field model's dissipation in the thermal film round the wire. It puts the catheter within 1 K of the water at every
power — the antenna's own loss is negligible and the cross-flow at 5 L/min is strong ($h_c\approx7800$).

## 2. Assumptions and what they cost

| assumption | note |
|---|---|
| Source from the **homogeneous-gel** E-field model, antenna in gel | The antenna is in water. The near-field shell ($\rho<1$ mm) that takes 70–90 % of $P_\text{del}$ in the model is water, less lossy than the gel ($\sigma_\text{eff}$ 1.2 vs 1.84 S/m, $|E|^2\propto1/|\varepsilon|^2$: ×0.39) — the model probably **under-estimates $P_\text{gel}$ by ~1.5–2×**. Magnitude only; see §6 for the shape |
| Antenna: 27 mm bare wire, $k_e=k_\text{gel}$, 27 mm braid return arm, 0.5 mm gap | From `Experimental_setup_table.xlsx` and v2.1. **None of these four describe the photographed antenna** (§5.1, §6.1) |
| $h_\text{top}$ from a laminar plate with the flow section = the gel face | The channel section is unknown; the antenna in the flow makes the boundary layer turbulent. $h$ could be a few × either way. Matters only in the top ~1 mm (§6.3) |
| Gel starts in equilibrium with 37 °C water | Photographed run: bath and gel at 25 °C (v2 README). Irrelevant to the *rise* (linear), decisive for what the pigment shows |
| 30 × 30 × 30 mm block | The photographs show a block ~30 mm wide, ~19–20 mm tall in the viewing plane; depth unknown. Over 6 min the heat has not reached the walls, so this only affects the steady state |

## 3. Verification — done and not done

Done, in the code: FTCS stability limit enforced; PCG steady state asserted converged (rel. residual $10^{-9}$);
frames land exactly on multiples of 2 s; `'L_sheath',0` in the field model reproduces the truncated current.

**Not yet done** (to do before quoting absolute temperatures): energy audit ($\int q\,dV$ vs. stored energy + boundary
losses), grid convergence of the peak ($dh$ = 1 / 0.5 / 0.25 mm — the v2 README shows the peak of a bare-wire source
is grid-sensitive), and a check of the 3-D solve against the 2-D `gel_heat_2d.m` steady state.

## 4. Baseline results (as run: 27 mm, 0.5 mm gap, 37 °C water, `gel_heating.png`)

| gap | $P_\text{gel}$ of 20 W | peak at 50 W, 360 s | at | peak at 100 W | at 150 W |
|---|---|---|---|---|---|
| 0.5 mm | 3.21 W (16 %) | 50.3 °C | $x=0$, $z=-2$ mm | 64.7 °C | 79 °C |
| 1.0 mm | 2.59 W | 47.2 °C | $x=-0.5$, $z=-2.5$ | | |
| 2.0 mm | 1.89 W | 43.5 °C | $x=-0.5$, $z=-3$ | | |
| 3.0 mm | 1.47 W | 41.4 °C | $x=-0.5$, $z=-3$ | | |

Optimiser (`gel_sparing.png`): 5 °C water, 126 W, 360 s keeps 1.7 mm below 45 °C with a 50 °C lesion under it.
The shape in every case: hottest spot **at the feed** (or under the braid if the grid extends there), **1–3 mm deep,
attached to the surface**, decaying monotonically toward the tip with the standing-wave ripples of a 1.5 λ element.

---

## Part B — the photographs, and why the model disagrees

## 5. What the photographs show

### 5.1 Geometry and scale — the photographed antenna is **12 mm**, not 27 mm

`IMG_7907–7918` (Canon 5D II, 5616 × 3744) are 12 frames, 29–36 s apart, 13:24:43 → 13:30:15 on 2026-08-25:
a **332 s** series (that is v2's "at 332 s"). No colour at 13:24:43, a small spot at +29 s, so power-on is at
0 ± 15 s of the first frame. `mirrored.jpg` is **IMG_7910 = t ≈ 92 s, not the end of the run**. `photo_series.jpg`
is the whole sequence, cropped and mirrored so the feed is on the left as in the model.

Scale from the gel block (30 mm, table): block 3248 px, element 1312 px → **L = 12.1 mm** (12.5–13 with the
antenna 15 mm behind the front face). In `mirrored.jpg`: **105.8 px/mm**, feed (cap edge) at $x=430$ px, gel surface
at $y=355\pm10$ px, wire axis at $y=235$ → **axis 1.1 mm above the gel**, gap ≈ 0.65 mm for a 0.42 mm radius.
Braid OD ≈ 1.6 mm and element ≈ 0.7 mm with a stranded look → an RG-178-class cable with the PTFE left on the inner
conductor (to be confirmed — it matters, §6.2). `Experimental_setup_table.xlsx`'s 27 mm element is a different
antenna from the one in these pictures; this note uses 12 mm throughout.

### 5.2 The heated zone (12 mm scale)

| t (s) | frame | coloured zone $x$ · depth (mm) | hottest (teal) spot $x$ · depth | cool top layer under the element |
|---|---|---|---|---|
| 29 | 7908 | ~3 mm spot, $x$ 2–4 · 1.5–3.5 | (pink only) · ~2.5–3 | **1.1 mm** |
| 62 | 7909 | −1…8 · 0.9–6 | 4.2 · 4.5 (green just appearing) | 0.8 |
| 92 | 7910 | −0.3…7.6 · 0.9–6.3, green 0.8–6.3 · 1.8–5.1 | **3.1–3.6 · 3.3–3.4** | 0.5 |
| 182 | 7913 | −3…9.5 · 0.7–8.4+ | 3.3 · 3.6 ; second lobe $x$ 6–9 · 5–6 | 0.4 |
| 332 | 7918 | −4…12 · 0.6 to below the crop (>8.6) | 3.3 · ~4.5 (region, not a point) | 0.3 |

Colour sequence going **down** through the core ($x=2.6$–3.5 mm, frame 7910): clear 0–0.7 mm → pink 0.8–1.6 →
yellow ~2.0 → **green 2.3–4.2** → yellow 4.6 → pink 4.9–6.5 → clear. That is the monotone red→green play of a
cholesteric (TLC) pigment with no blue and no clear "hole", so (i) the green/teal is the true maximum, (ii) the top
layer is genuinely **below the pigment's start temperature**, not above its clearing point, (iii) the whole blob
spans only the few kelvin of the colour play. The pigment's start temperature and bandwidth are not on file.

Features the model must reproduce, in order of how robust they are:

1. **Depth**: maximum ~2.5–3 mm deep at 29 s (diffusion length $\sqrt{\alpha t}=2$ mm — so this is where the
   *source* peaks, not a thermal artefact), 3.3–3.6 mm from 92 s on; a cool layer of ~1 mm on top that **thins**
   to 0.3 mm by 332 s (an isotherm creeping toward a cooled surface — a pigment-free skin would stay constant).
2. **Position**: $x\approx3$–3.6 mm along the element, never at the feed; a **second, deeper and more diffuse lobe** at
   $x\approx6$–9 mm, 5–6 mm deep, from ~130 s; only a **faint** lobe under the braid.
3. **Growth**: the outer isotherm grows ~$\sqrt t$ (3 → 9 → 18 mm), the core stays put.
4. **Magnitude**: no boiling, no bubbles, gel intact at 50 W for 332 s.

## 6. Why the model looks different — tests and ranking

All tests: `photo_comparison.m` (50 W set, rise above the no-power steady state, antenna plane; "core" = region
where the rise ≥ half the peak rise; "s/p" = rise at the surface node ÷ peak rise in the same column). Photograph
targets: peak 2.5–3 / 3.4 / ~3.5 mm deep at 29 / 92 / 332 s, $x\approx3$–3.6 mm, s/p ≪ 1.

| case (12 mm unless stated) | $P_\text{gel}$ | 29 s: peak K @ ($x,z$) | 92 s | 332 s | s/p 92 s |
|---|---|---|---|---|---|
| **V3 as run**: 27 mm, arm 27 mm, $k_e=k_\text{gel}$, gap 0.5 | 2.1 W | 8.9 @ (0.5, −1.0) | 12.1 @ (0, −2.0) | 15.8 @ (**−3.0**, −3.5) | 0.46 |
| 27 mm, no braid arm | 2.7 | 18.9 @ (0.5, −0.5) | 20.4 @ (0.5, −0.5) | 21.2 @ (0.5, −0.5) | 0.80 |
| 12 mm, arm 12 mm, $k_e=k_\text{gel}$, gap 0.65 | 2.5 | 14.9 @ (2.5, −1.0) | 18.4 @ (2.5, −1.5) | 21.1 @ (2.0, −2.0) | 0.53 |
| 12 mm, arm 4 mm, $k_e=k_\text{gel}$ | 2.5 | 20.3 @ (3.0, −1.0) | 23.7 @ (3.0, −1.0) | 25.7 @ (3.0, −1.5) | 0.57 |
| 12 mm, arm 4 mm, **$k_e=k_\text{gel}/2$** | 3.1 | 12.5 @ (0.5, −1.5) | **19.3 @ (2.0, −2.5)** | **25.9 @ (2.5, −3.5)** | **0.32** |
| 12 mm, arm 4 mm, $k_e=k_\text{gel}/2$, **surface pinned** | 3.1 | 10.9 @ (1.5, −2.0) | **17.2 @ (2.5, −3.0)** | **23.3 @ (2.5, −4.0)** | 0 |
| 12 mm, arm 4 mm, quarter-wave $k_e=\pi/2L$ | 3.1 | 13.1 @ (−2.0, −1.0) | 19.4 @ (1.5, −2.5) | 26.0 @ (1.5, −3.5) | 0.31 |
| 12 mm, arm 4 mm, $k_e=k_\text{gel}$, **$E_x$ only** (diagnostic) | 1.1 | 6.0 @ (0, −1.0) | 7.1 @ (0, −1.5) | 8.3 @ (**8.0, −4.0**) | 0.54 |

Two further tests on the 27 mm geometry (scratch runs, 360 s): a **plane-wave-like source** $e^{-2\alpha z}$ with the gel's
9.7 mm power depth, uniform under the element, still peaks at 4.5 mm with the core's top at 1 mm — a surface-entering
source plus a pinned surface *cannot* put the maximum deeper than ~5 mm in 6 min; and **switching the power off** sinks
the core from 3.2 to 4.4 mm in 60 s but the gel is back within 2 K of the water by then. Neither is needed once the
antenna is right.

### 6.1 Inputs that do not describe the photographed run (explains the *position*)

- **27 mm element with a 27 mm cosine return arm on the braid.** The arm carries the full feed current 27 mm back
  along the coax, so the current — and the charge at its ends — is largest around the junction; the model's hottest
  spot lands *at* the feed, and once the grid extends under the braid, *under the braid* ($x=-3$ mm at 332 s).
  Removing the arm on the 27 mm element moves the peak from $x=-0.5$ to **7.5 mm** (scratch run). The photographs
  show only a faint lobe under the braid, so the arm is over-weighted; 4 mm reproduces "faint".
- **With the actual 12 mm element** the model's source already peaks at $x=2.5$–3.5 mm for every $k_e$ tried
  (rows 3–7): for $k_\text{gel}L=4.1$ rad the current has a node ~3 mm from the feed where the charge, and hence the
  $1/\rho$ near field, is largest; the current antinode at ~7.4 mm gives the $E_x$ lobe — which is exactly where the
  $E_x$-only diagnostic puts its peak at 332 s (**$x=8$ mm, 4 mm deep**) and where the photographs grow their
  **second, deeper lobe**. The along-element pattern is therefore *right in kind*; it was the inputs that were wrong.
- Gap 0.5 vs 0.65 mm, bath 37 vs 25 °C: negligible for the shape.

### 6.2 The depth — the charge-driven near field is too strong (physics)

The field of a thin wire is $\mathbf E=-j\omega\mathbf A-\nabla\phi$. Near the wire the scalar-potential part is the
quasi-static field of the line charge, $E_\rho\simeq q_l/(2\pi\varepsilon_c\rho)$, so $q\propto1/\rho^2$ — **it peaks
at the gel surface directly under the wire** ($\rho=h$) and falls 17× by 3.4 mm depth; the vector-potential part
($E_x$, from the current) is only logarithmic in $\rho$ and spreads smoothly. Charge continuity gives
$|q_l|=|dI/dx'|/\omega\sim|k_e||I|/\omega$: **the surface spike scales with $|k_e|^2$, the smooth part does not.**

V3 takes $k_e=k_\text{gel}=344-52j$ rad/m (bare wire in gel, King & Smith). The element is not in gel and, by the look
of it, not bare: it is a PTFE-covered conductor in **water**. King's insulated-antenna formula
[King, Trembly & Strohbehn, IEEE Trans. MTT-31 (1983) 574; King & Smith, *Antennas in Matter*, ch. 8]

$$k_L=k_2\Big[1+\frac{H_0^{(1)}(k_3b)}{k_3b\,H_1^{(1)}(k_3b)\,\ln(b/a)}\Big]^{1/2},\qquad
k_2=k_0\sqrt{\varepsilon_\text{PTFE}},\ k_3=k_0\sqrt{\varepsilon_\text{water}}$$

gives, in the small-argument limit, **$k_L\approx129-34j$** for an RG-178 conductor ($a=0.155$, $b=0.42$ mm,
$\varepsilon_\text{PTFE}=2.1$, water $73-8.8j$) — $\lambda_L=49$ mm, so the 12 mm element is a **quarter wave**
(the v1 write-up's hypothesis H1, now with a derivation) — and ≈ $212-71j$ for a thin 0.15 mm skin on a 0.45 mm
conductor. $k_\text{gel}/2=172-26j$ sits in that bracket; a *bare* wire in water would be $440-26j$ and would make the
disagreement worse. Halving $k_e$ (rows 5–7) cuts the surface share s/p from 0.57 to 0.31 and moves the peak to
**2.5 / 3.5 mm at 92 / 332 s**, against 3.4 / ~3.5 in the photographs, at $x=2$–2.5 mm; the $E_x$ lobe then
carries a comparable share of the power, so the second lobe becomes visible. This is the single physical change that
the depth series demands.

Two interface effects, both *ignored* by the homogeneous model, push the same way: (i) the normal component obeys
$\varepsilon_wE_{z,w}=\varepsilon_gE_{z,g}$, so $E_z$ is **1.66× larger** in the gel just below the surface than in
the water above it — the surface spike is, if anything, under-estimated, which makes a *small* $k_e$ even more
necessary; (ii) in water the near-field shell dissipates ~4× less than the model's gel-filled shell, so more of
$P_\text{del}$ reaches the gel (magnitude, §6.4).

### 6.3 The cool top layer — the surface is better held than $h=1460$ with a spike under it

Even at Bi = 73 the model's surface node rises to 50–60 % of the peak (s/p, rows 1–4): the $1/\rho^2$ source sits in
the first half-millimetre and no finite $h$ holds it down. Pinning the surface (row 6) gives a detached core with its
top 0.5–1 mm down, as photographed, and the isotherm then creeps up as the bulk warms — the observed 1.1 → 0.3 mm
thinning. With the smaller $k_e$ the surface rise is already only 0.3 of the peak, so a moderately larger $h$
(turbulent wake of the antenna, pump jet) or a bath a few K below the pigment's start temperature does the rest. The
photographed bath was 25 °C; the pigment start is unknown; both are needed to say which.

### 6.4 Magnitude — plausible, no longer the problem it was

At 12 mm, 10 W delivered, 2.5–3.1 W in the gel, the peak rise at 332 s is 21–26 K → **46–51 °C from a 25 °C bath**: an
intact gel and a pigment playing in the 30–45 °C range are both consistent. v2's 147–185 °C came from $\eta=0.75$
(37 W delivered) and the surface spike; $\eta=0.2$ and the shallow-gel share fix most of it. Remaining uncertainty:
$\eta$ is stated, not measured; $P_\text{gel}$ is probably 1.5–2× low (§6.2 ii); the pigment band is unknown.

### 6.5 Ruled out

| candidate | why not |
|---|---|
| TLC above its clearing point (hot core looks clear) | colour order clear→pink→yellow→green going down: the layer next to the clear zone is the *cold* end of the band |
| pigment-free / settled surface layer | its thickness would be constant; the observed layer thins 1.1 → 0.3 mm |
| photograph taken after power-off | frames every 30 s during the run; the blob grows monotonically |
| lower gel loss ($\varepsilon''=8.8$), plane-wave illumination, longer exposure | each moves the peak ≤ 1–2 mm and leaves it attached to the surface |
| total internal reflection trapping in the water layer ($\theta_c=51°$) | flattens the lateral pattern, does not create a sub-surface maximum |
| thermal convection in the gel, perfusion, boiling | solid 1 % gel; nothing near 100 °C |

## 7. What to do next

Measure / confirm (one line each decides a section above):

1. **The element in the photographs**: length (12 vs the table's 27 mm), cable type, and whether the PTFE was left on
   the inner conductor — §6.2 hinges on it. If the 27 mm element is the *new* antenna, the photographs cannot be used
   to validate the 27 mm model at all.
2. **Pigment datasheet**: start temperature, bandwidth, clearing point. Turns the colour bands into isotherms.
3. **Bath temperature and flow geometry** on the day (25 °C tap water per v2; channel section for $h_\text{top}$).
4. **Power-on time** relative to IMG_7907 (±15 s now), and forward/reflected power ($\eta$).
5. A thermocouple 3 mm deep under $x=3$ mm would calibrate the whole chain.

Model changes, in order:

1. Done today: `gel_heating_solve` accepts `'L'`, `'L_sheath'`, `'ke'` (own cache files). Solve the photographed
   configuration — `gel_heating_solve('L',12e-3,'L_sheath',4e-3,'ke',172-26j,'gap',0.65e-3)` — and compare the
   isotherm of the pigment's start temperature with `photo_series.jpg` frame by frame (linear in $P$, so $\eta$ and
   the start temperature trade off; the *shape* does not).
2. Replace the asserted $k_e$ by King's insulated-antenna $k_L$ computed from the cable geometry (a 10-line function;
   the bracket above is the small-argument limit).
3. Bring the water back above the gel (the v2 layered field, `../v2_water_above/matlab/layered_field.m`): it fixes
   the near-field shell's loss, the 1.66× $E_z$ jump, and the power split — magnitude and the top millimetre.
4. Verification items of §3 before quoting absolute temperatures.

## References

R. W. P. King, B. S. Trembly, J. W. Strohbehn, "The electromagnetic field of an insulated antenna in a conducting
or dielectric medium," *IEEE Trans. Microwave Theory Tech.* MTT-31 (7), 574–583, 1983.
R. W. P. King, G. S. Smith, *Antennas in Matter*, MIT Press 1981, ch. 8 (insulated antennas) —
[archive.org](https://archive.org/details/antennasinmatter00rono).
F. P. Incropera, D. P. DeWitt, *Fundamentals of Heat and Mass Transfer*: ch. 7 (flat plate, Churchill–Bernstein), ch. 9.
Internal: `E_field_2d_notes.md`; `../v2_water_above/README.md` (braid arm, bath, the 332 s series);
`../E_field_model_writeup.md` §6.3 (hypothesis H1); `../Heatmodel/Heat_model_writeup.md` §7 (why 3-D).
