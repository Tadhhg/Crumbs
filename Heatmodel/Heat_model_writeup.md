# Heating of the Phytagel slab: generator power in, temperature out

**Working document — model derivation, sources, numerical implementation and results**
Date: 2026-09-06 · Code: `Heatmodel/matlab/` · Figures: `Heatmodel/figures/` · Numerical output: `Heatmodel/results/`

Companion to [`../E_field_model_writeup.md`](../E_field_model_writeup.md), which supplies the
electromagnetic source. Read that first; this document assumes its geometry, its coordinates and
its $\tilde{\mathbf E}$.

---

## 0. The answer, up front

Cycle-averaged microwave deposition heats the gel; conduction carries the heat to the blood, which
holds the surface at 37 °C. The governing equation is

$$\rho c\,\frac{\partial T}{\partial t} \;=\; k\,\nabla^{2}T \;+\; q(\mathbf r),
\qquad q(\mathbf r)=\tfrac12\,\sigma_\text{eff}\big|\tilde{\mathbf E}(\mathbf r)\big|^{2}
=\sigma_\text{eff}E_\text{rms}^{2}$$

with $T=37$ °C on the gel/blood surface, $T(\mathbf r,0)=37$ °C, and **no perfusion term** — the gel
is a phantom, so the blood only appears as a boundary condition (§2.2).

Everything in that problem is linear and the source is proportional to the delivered power. So one
solve answers every power:

$$\boxed{\;T(\mathbf r,t)\;=\;37\ \text{°C}\;+\;\frac{\eta\,P_\text{gen}}{150\ \text{W}}\;\Theta(\mathbf r,t)\;}
\qquad \eta=\frac{P_\text{del}}{P_\text{gen}}=0.75$$

$\Theta(\mathbf r,t)$ is the temperature rise computed once, for 150 W delivered, and stored in
`results/heat_model.mat`. **The scaling in $P$ is exact**, not a fit — within the constant-property
model, whose limits are quantified in §10F. Note the contrast with the field model, where
$|\mathbf E|\propto\sqrt{P}$: temperature rise goes with $P$ because heating goes with $E^2$.

For an arbitrary generator waveform $P(t)$ — pulsed, ramped, duty-cycled — the same linearity gives
the answer by Duhamel superposition of step responses (§11.2); no new solve is needed.

**The deliverable is `T_gel.m`:**

```matlab
addpath('Heatmodel/matlab');
run_heat_model                       % once: builds results/heat_model.mat (~15 min)

T_gel(12, 60)                        % peak gel temperature at 12 W generator, 60 s  [°C]  -> 97.7
T_gel(12, 60, 6e-3, -2e-3)           % at x = 6 mm, 2 mm deep, on y = 0
T_gel(12, 60, 6e-3, 1e-3, -2e-3)     % anywhere in the block
T_gel(12, Inf)                       % steady state
P_for_T(100, 60)                     % generator power that just boils it at 60 s  [W]  -> 12.46
time_to_T(50, 100)                   % seconds until 100 °C at 50 W               [s]  ->  3.39

heat_scenario(12, [], 60)            % one bench run: power, pulse period, stimulation
heat_scenario(40, 1, 60, 'duty',0.3) % time.  Draws the figures and writes a CSV (Sec. 11.5)
```

### The headline

At **200 W generator power** — the setting the field model was normalised to — the gel reaches
43 °C in **0.04 s** and 100 °C in **0.72 s**, and would settle, if it could, at about 1280 °C.
**The phantom cannot absorb anything like full power.**

| what you want | generator power |
|---|---|
| 100 °C after a 60 s exposure | **12.4 W** (8.6 W if the surface is convective rather than pinned) |
| 100 °C, run indefinitely | **10.2 W** (7.1 W convective) |
| stay below 43 °C indefinitely | **1.0 W** (0.7 W convective) |

Per generator watt the peak rises **0.39 K/W at 1 s, 2.79 K/W at 10 s, 5.07 K/W at 60 s** and
6.19 K/W at steady state. The peak is **82 % of its final value by one minute**.

The hottest gel is **2–4 mm below the surface, never at it**: the blood pins the surface at 37 °C
and the catheter, cooled by the same flow, adds no hot spot of its own (§4, §9.3). At 12.4 W for
60 s the zone that would ablate tissue (CEM43 ≥ 240 min) is **9.0 mm deep and 15 mm long**.

Two warnings that the sensitivity study makes unavoidable. **The 37 °C surface you asked for is an
optimistic bound** — it is the $h\to\infty$ end of a range that costs a factor 2.6 in allowed power
(§10B). And **pulsing does not buy headroom**: at the same mean power, 30 % duty at 1 Hz raises the
peak by 2 °C and at 0.02 Hz by 79 °C (§11.2).

---

## 1. What was specified, and what I had to assume

**Given:**

| item | value |
|---|---|
| heat sink | blood, modelled as water at **37 °C**, flowing at **5 L/min** |
| surface treatment | "you can basically model the surface as just a constant 37 °C" |
| the catheter | lies on the gel and **is cooled by the same water** |
| region of interest | 30 mm × 30 mm below the antenna; the rest is "taken care of by homeostasis" |
| the source | the 2.45 GHz field model, unchanged |

**Assumed, because it was not specified.** Each is an explicit parameter and each is swept in §10.

| | assumption | baseline | why | effect |
|---|---|---|---|---|
| T1 | gel thermal conductivity $k$ | 0.60 W m⁻¹K⁻¹ | >98 % water; measured agar phantoms bracket it (§5) | §10F |
| T2 | gel specific heat $c$ | 4180 J kg⁻¹K⁻¹ | water at 37 °C | §10F |
| T3 | gel density $\rho$ | 1000 kg m⁻³ | **held identical to the field model**, so SAR and $q$ agree exactly | — |
| T4 | flow geometry above the gel | 30 × 30 mm channel, patch at the leading edge | needed to turn 5 L/min into a velocity; nothing else was given | **the largest**, §10B |
| T5 | what is under the gel | adiabatic (air-backed dish) | matches assumption A5 of the field model | small, §10D |
| T6 | walls of the region of interest | held at 37 °C | this is exactly what "homeostasis" means as a boundary condition | small, §10E |
| T7 | catheter outer diameter / wall | 1.5 mm, 0.25 mm of PTFE-like polymer | typical ablation catheter | §10C |
| T8 | delivery efficiency $\eta$ | 0.75 | 200 W forward − 50 W max reflected, the field model's worst case | exact linear scaling |
| T9 | gel permittivity temperature | 25 °C | the field model's baseline, kept for exact consistency | 4 %, §10H |

**"30 mm × 30 mm" — what I took it to mean.** The square is the vertical plane $y=0$: 30 mm along
$x$ (centred on the antenna, so $-9 \le x \le 21$ mm) by the full 30 mm depth of the gel. That is
precisely the plane the field model reports $|\mathbf E|$ on, so the two documents show the same
picture.

**But the conduction is solved in three dimensions**, on the 30 × 30 × 30 mm block containing that
plane. This is not optional. The thermal diffusion length is

$$\sqrt{4\alpha t}=2.4\ \text{mm at 10 s},\qquad 5.9\ \text{mm at 60 s},\qquad \alpha=\frac{k}{\rho c}=1.435\times10^{-7}\ \text{m}^2\text{/s}$$

which is comparable to the transverse ($y$) half-width of the deposition lobe. A 2-D solve in the
$y=0$ plane would forbid heat from spreading sideways and would overstate the temperature
substantially. The block captures 94.7 % of the power absorbed in the whole gel.

---

## 2. The equation, and the two terms that are *not* in it

### 2.1 Why the microwave ripple disappears

The field magnitude oscillates at $2f = 4.9$ GHz (field-model §10.1). The thermal response cannot
see it: the ratio of the thermal time constant to the 204 ps ripple period is about $10^{10}$. Only
the cycle average survives, and it is exact:

$$\langle \mathbf J\!\cdot\!\mathbf E\rangle_t=\tfrac12\sigma_\text{eff}\big|\tilde{\mathbf E}\big|^{2}
=\sigma_\text{eff}E_\text{rms}^{2}\;\equiv\;q(\mathbf r)\quad[\text{W/m}^3]$$

$q = \rho\cdot\text{SAR}$, and $\rho$ is 1000 kg/m³ in both models, so $q$ and the SAR map of the
field model are the same object in different units.

### 2.2 There is no perfusion term

The Pennes bioheat equation used throughout the ablation literature,

$$\rho c\,\partial_t T = k\nabla^2T+\omega_b\rho_bc_b\,(T_a-T)+q+q_\text{met},$$

carries a volumetric perfusion sink $\omega_b\rho_bc_b(T_a-T)$ because *tissue* is threaded with
capillaries. **A Phytagel phantom is not perfused.** Putting a perfusion term in would be inventing
a cooling mechanism the bench experiment does not have. The blood in this problem is entirely
outside the gel, and it enters only through the surface boundary condition — which is exactly how
the brief describes it. The metabolic term $q_\text{met}$ (≈ 1 kW/m³ in muscle) is likewise absent,
and would in any case be six to eight orders of magnitude below $q$.

This is worth saying explicitly because it is the main reason phantom temperatures run *hotter*
than the tissue they mimic at the same power. Muscle perfusion is about 37 mL min⁻¹kg⁻¹ [T4], i.e.
$\omega_b\approx6.7\times10^{-4}$ s⁻¹, so $\omega_b\rho_bc_b\approx2.5$ kW m⁻³K⁻¹ and the perfusion
time constant is $\rho c/(\omega_b\rho_bc_b)\approx1700$ s. That is the same order as the conduction
time constant of this 30 mm block ($1/[\alpha(\pi/2d)^2]=2540$ s, §7.2), so **in real tissue the two
sinks carry comparable weight, and the phantom is missing one of them.** Over the first minute,
where conduction has barely reached 6 mm, the omission changes little; over minutes it does not.

### 2.3 The equation actually solved

Writing $\theta = T - T_\text{blood}$ makes every boundary condition homogeneous:

$$\frac{\partial\theta}{\partial t}=\alpha\nabla^{2}\theta+\frac{q(\mathbf r)}{\rho c},
\qquad \theta(\mathbf r,0)=0$$

| boundary | condition | why |
|---|---|---|
| $z=0$ (gel/blood) | $\theta=0$ **(baseline)** or $-k\,\partial_z\theta = h\theta$ | as specified; the convective case is run alongside, §3 |
| $z=-30$ mm (floor) | $\partial_z\theta=0$ | air-backed dish (field-model A5) |
| four side walls | $\theta=0$ | "minor temperature rises taken care of by homeostasis" |

Because the problem is linear, homogeneous in its boundary conditions and zero in its initial
condition, and because $q\propto P_\text{del}$, the solution is *exactly* proportional to
$P_\text{del}$. That single fact is what turns one 3-D transient solve into the function the brief
asked for.

---

## 3. The heat sink: 5 L/min of 37 °C water

### 3.1 What the brief gives, and what it does not

The flow rate is given; the cross-section it flows through is not, and without that there is no
velocity and therefore no heat-transfer coefficient. I take the channel above the gel to have the
same 30 × 30 mm section as the region of interest, with the patch starting at the leading edge of
the flow (assumption T4). `matlab/convection_h.m` then gives, for water at 37 °C
($\rho=992.99$, $c_p=4178$, $k=0.6282$, $\mu=6.929\times10^{-4}$, $\Pr=4.609$; Incropera
Table A.6):

$$U=\frac{\dot V}{A}=\frac{8.333\times10^{-5}}{9\times10^{-4}}=0.0926\ \text{m/s},
\qquad \mathrm{Re}_L=\frac{UL}{\nu}=3981,\qquad \mathrm{Re}_{D_h}=3981$$

**Which correlation?** The laminar entrance lengths are $L_\text{hydro}=0.05\,\mathrm{Re}\,D_h=5.97$ m
and $L_\text{therm}=0.05\,\mathrm{Re}\,\Pr\,D_h=27.5$ m. Both are enormous compared with the 30 mm
patch, so over the gel the flow is an **undeveloped boundary layer on a flat plate**, not duct flow.
That fixes the baseline:

| correlation | $\mathrm{Nu}$ | $h$ [W m⁻²K⁻¹] | comment |
|---|---|---|---|
| **flat plate, laminar** $\;0.664\,\mathrm{Re}_L^{1/2}\Pr^{1/3}$ | 69.8 | **1460** | **baseline** — the entrance-length argument above |
| flat plate, turbulent $(0.037\mathrm{Re}_L^{4/5}-871)\Pr^{1/3}$ | 0 | 0 | $\mathrm{Re}_L$ is far below transition; listed to show it does not apply |
| duct, Dittus–Boelter $\;0.023\,\mathrm{Re}_D^{4/5}\Pr^{0.4}$ | 32.2 | 673 | if the patch sat far downstream in a developed duct |
| duct, fully developed laminar, one wall heated | 4.9 | 103 | the pessimistic bound |

So $h$ is somewhere between about 100 and 1500 W m⁻²K⁻¹ depending on a geometry that has not been
specified. **This is the single largest uncertainty in the thermal model**, the counterpart of the
standoff gap in the field model, and §10B sweeps it.

### 3.2 Is "constant 37 °C" defensible?

The instruction was to model the surface as held at 37 °C, i.e. $h\to\infty$. The relevant
comparison is the **cooling length** $k_\text{gel}/h$ — the thickness of gel whose conduction
resistance equals the surface film resistance:

$$\frac{k_\text{gel}}{h}=\frac{0.60}{1460}=0.41\ \text{mm}\quad\text{(baseline }h\text{)},
\qquad \frac{0.60}{103}=5.8\ \text{mm}\quad\text{(worst case)}$$

Against the ~1.3 mm $1/e$ depth of the power deposition, 0.41 mm is a modest correction and 5.8 mm
is not. So: **the Dirichlet surface is a lower bound on the gel temperature**, tight if the flow is
as assumed and loose if it is not. Both are computed — see §10B for the numbers and
`results/heat_model_robin.mat` for the full convective solution.

### 3.3 Ignoring the coolant temperature rise is safe

The brief says to ignore it. Checking:
$\dot m = \rho\dot V = 0.0828$ kg/s, so the bulk rise per watt removed is
$1/(\dot m c_p) = 2.9\times10^{-3}$ K/W. Even at the extreme of the model removing 110 W, the water
leaves 0.32 K warmer. Ignoring it is right.

---

## 4. The catheter

The reminder in the brief is the important physical point of the whole model, so it is worth
stating in full.

The catheter lies on the gel along $0\le x\le 12$ mm at $y=0$, and it is bathed in the same 37 °C
flow. Three things follow.

**It is not a heat source.** The metal and the polymer wall have no dielectric loss worth the name
at 2.45 GHz. Nothing in the applicator heats except by conduction from the gel.

**It is not an insulating blanket either.** The catheter outer surface is a small cylinder in cross
flow: $\mathrm{Re}_D = UD/\nu = 199$, and Churchill–Bernstein gives $\mathrm{Nu}_D=14.2$, hence
$h_\text{out}=5940$ W m⁻²K⁻¹ — *higher* than the flat surface, because the diameter is small. In
series with the wall conduction resistance $t/k = 0.25\text{ mm}/0.25\text{ W m}^{-1}\text{K}^{-1}
= 1.0\times10^{-3}$ m²K/W, the effective coefficient of the footprint is

$$h_\text{cath}=\Big(\tfrac{1}{5940}+1.0\times10^{-3}\Big)^{-1}=857\ \text{W m}^{-2}\text{K}^{-1}$$

which is the same order as the 1460 of the bare surface. The catheter footprint is a *comparably
cooled* patch, not a hole in the heat sink.

**Therefore the hottest gel is not at the surface.** The deposition is largest right under the
conductor, at $z=0$; the cooling is also strongest there. The maximum is pushed to a few millimetres
depth, and the surface itself stays at 37 °C. This is the same mechanism that irrigated microwave
ablation catheters use deliberately — including the Qian et al. work [R13 of the field model] that
this phantom recipe comes from — to put the lesion below the endocardial surface instead of
scorching it.

In the baseline (surface pinned at 37 °C) the catheter needs no separate treatment: a cooled object
resting on a plane already held at 37 °C changes nothing. `cath_bc` therefore defaults to
`'transparent'`, and the two informative alternatives — `'cooled'` at $h_\text{cath}$, and
`'insulating'` (an adiabatic strip, the pessimistic bracket) — are run in §10C.

---

## 5. Material properties — every number and its source

### 5.1 The gel

>98 % water by mass, so the thermal properties are water's, pulled down slightly by the polymer.

| quantity | value | source and reasoning |
|---|---|---|
| $\rho$ | 1000 kg m⁻³ | **held equal to the field model's SAR density** so $q=\rho\,\mathrm{SAR}$ is exact. Water at 37 °C is 993. |
| $c$ | 4180 J kg⁻¹K⁻¹ | water at 37 °C is 4178 (Incropera Table A.6 [T1]). Measured agar tissue phantoms: 4090–4130 [T2]. |
| $k$ | 0.60 W m⁻¹K⁻¹ | water at 37 °C is 0.628 [T1]. Measured water–agar gels fall a few percent below water and decrease with polymer fraction [T3]; measured agar tissue phantoms give 0.52 ± 0.06 (brain-equivalent) and 0.57 ± 0.10 (muscle-equivalent) [T2]. At 1 % w/v the value should sit near the top of that range. Swept 0.50–0.628 in §10F. |
| $\alpha=k/\rho c$ | $1.435\times10^{-7}$ m² s⁻¹ | derived |
| $\sigma_\text{eff}$ | 2.723 S m⁻¹ | field model §4.3, unchanged |

For comparison, muscle is $k=0.49$, $c=3421$, $\rho=1090$ (IT'IS database [T4]) — run as a
sensitivity case in §10F.

### 5.2 The coolant

Water at 37 °C, from Incropera Table A.6 [T1], implemented in `matlab/water_properties.m` as a
table at 5 K intervals from 290 K to 370 K with `pchip` interpolation (the bottom two rows exist so
that room-temperature coolant is interpolated rather than extrapolated). The table stores $\Pr$
independently of $c_p$, $\mu$ and $k$, and the code asserts $\Pr = c_p\mu/k$ at every lookup — a
transcription check that passes to better than 0.2 %.

| $\rho$ | $c_p$ | $k$ | $\mu$ | $\nu$ | $\Pr$ |
|---|---|---|---|---|---|
| 992.99 kg m⁻³ | 4178 J kg⁻¹K⁻¹ | 0.6282 W m⁻¹K⁻¹ | $6.929\times10^{-4}$ Pa s | $6.978\times10^{-7}$ m² s⁻¹ | 4.609 |

Real blood ($\rho=1050$, $c_p=3617$, $k=0.52$, $\mu=3.5\times10^{-3}$ Pa s [T4]) is available as
`water_properties('blood')`. Its $\nu=3.33\times10^{-6}$ m² s⁻¹ is 4.8× water's and its
$\Pr=24.4$ is 5.3× water's, so at fixed volumetric flow
$\mathrm{Nu}\sim\mathrm{Re}^{1/2}\Pr^{1/3}$ falls by 20 % and $h=\mathrm{Nu}\,k/L$ by a further
17 % on the lower conductivity — **$h$ drops about a third**, to roughly 960 W m⁻²K⁻¹. The brief
specified water, so water is the baseline, and the $h$ sweep of §10B covers the difference; but if
you ever run this with real blood, use the lower coefficient.

---

## 6. From the field model to the source term

### 6.1 What had to be recomputed, and why

`matlab/layered_field.m` in the field model keeps the complex phasor field on the $y=0$ plane and
throws the rest away. The heat equation needs $q$ in **three** dimensions. Holding the full 3-D
phasor field over the 256 mm FFT window at 0.25 mm would be about 200 GB, so
`Heatmodel/matlab/sar_volume.m` runs the same solver and keeps
$q=\tfrac12\sigma_\text{eff}|\tilde{\mathbf E}|^2$ on the 30 × 30 × 30 mm block instead.

It is the same physics, the same Weyl expansion, the same TE/TM transmission-line recursion, the
same branch handling, the same FFT normalisation and the same delivered-power normalisation. Two
differences of substance:

* the field is evaluated on **nodes** ($z=0,-0.25,\dots,-30$ mm) rather than cell centres, so that
  it lands exactly on the thermal grid and no interpolation is needed;
* the depth integral for the absorbed power uses **Simpson's rule** on those 121 nodes rather than
  a midpoint rule on 120 cells, because the integrand varies on a 1.3 mm scale near the surface.

Duplicating a solver is a risk, so it is **checked rather than asserted** — see §8, test 5.

### 6.2 The feed singularity — a real problem this model exposes

The field model prescribes the antenna current as $I(x)=I_m\sin[k_e(L-x)]$ on $0\le x\le L$ and zero
outside. With $L=\lambda_e/4$ that current is **maximum at the feed** and is then truncated. Charge
continuity, $\partial_x I + j\omega\rho=0$, converts that step into a point charge at the feed:

$$Q=\frac{I_m}{\omega}=\frac{4.02}{1.539\times10^{10}}=2.6\times10^{-10}\ \text{C}$$

Its quasi-static field 0.5 mm away in a medium with $\varepsilon'=74.3$ is
$Q/(4\pi\varepsilon_0\varepsilon' r^2)=1.3\times10^{5}$ V/m, against the
$2.4\times10^{5}$ V/m amplitude the solver actually returns at $(0,0,0)$ — the same size to within
a factor of two, which is all a bare point-charge estimate can claim next to an interface. What
identifies the mechanism beyond doubt is the asymmetry: at the **tip** the current goes smoothly to
zero, there is no truncation, no point charge, and the solver shows **no spike at all** —
$6.3\times10^{8}$ W/m³ at the tip against $7.7\times10^{10}$ at the feed, a factor of 120, for two
ends of the same wire at the same height.

Is it physical? Partly. In a real applicator the feed current does not stop — it continues onto the
outside of the coaxial outer conductor. The field model does not include that sheath (its §11,
limitation 3, which already predicted it would misrepresent the region behind the feed). A real feed
junction *is* a hot spot — the "backward heating along the shaft" of the catheter-ablation
literature — but neither its magnitude nor its shape is predicted correctly here.

This never mattered for the field model, which quoted values at $x=6$ mm. It matters a great deal
for a temperature model, whose natural output is a maximum. The treatment adopted is:

1. **The source is left exactly as the field model gives it.** No hidden edits; `T_gel` and `E_gel`
   describe the same field everywhere.
2. **The reported peak excludes a 2 mm sphere about the feed** (`thermal_config.r_feed_excl`).
3. **The alternative is computed and reported**: `feed_fix='cap'` limits $q$ inside that sphere to
   the largest value on its surface, which removes the point charge, and the difference is given in
   §9.2 and §10I.

The evidence that this is the right call, and that it does not contaminate the rest of the model:

| $\Delta h$ | power in a 2 mm sphere at the feed | peak $\Theta_\text{ss}$, feed excluded | peak $\Theta_\text{ss}$, unrestricted | $\Theta_\text{ss}$ at (6, 0, −2 mm) |
|---|---|---|---|---|
| 1.00 mm | 50.3 W | 1118 K | 1118 K | 884.5 K |
| 0.50 mm | 24.6 W | 1196 K | 1325 K | 928.6 K |
| 0.25 mm | 21.5 W | 1239 K | **1535 K** | **939.9 K** |

The unrestricted peak **diverges** under refinement — the signature of a genuine singularity. The
value at a point 6 mm away **converges**, and it is within 0.8 % of the value obtained with the
singularity removed (932.4 K capped vs 939.9 K raw). **More than about 3 mm from the feed, the
artefact changes nothing.** Everything quoted in this document is either far from the feed or
explicitly flagged.

---

## 7. Numerical implementation

`matlab/heat_grid.m`, `heat_lap.m`, `heat_solve.m`, `heat_steady.m`.

### 7.1 Discretisation

Node-centred second-order finite differences on a uniform 0.25 mm grid — the same spacing as the
field model's $\Delta x$, so the EM solution is sampled directly onto the thermal nodes.
121 × 121 × 121 = 1.77 M nodes. Nodes sit **on** the boundaries, so a boundary node owns half a
control volume in the normal direction. Boundary conditions use ghost nodes:

| condition | ghost |
|---|---|
| Dirichlet | node removed from the unknowns |
| adiabatic | ghost = mirror of the first interior node |
| Robin at $z=0$ | with depth $s=-z$: $k\,\partial_s\theta|_0 = h\theta_0$, so $\theta_g=\theta_2-2\Delta h\,\theta_1/k$ |

The Robin form deserves a note. Substituted into the ordinary three-point stencil it gives, at the
surface node,

$$\left.\nabla^2_z\theta\right|_1=\frac{2(\theta_2-\theta_1)}{\Delta^{2}}-\frac{2h\theta_1}{k\Delta}$$

which is **identical** to the finite-volume energy balance on the half cell (the factor $\frac12$
divides out of both sides). The scheme is therefore conservative as well as second-order, which is
what makes the energy audit in §7.3 meaningful.

### 7.2 Time integration and the steady state

Forward Euler, with $\Delta t=\min\big(\text{CFL}\cdot\Delta t_\text{stab},\ t_\text{next}-t\big)$ and

$$\Delta t_\text{stab}=\frac{\Delta^{2}}{\alpha\,(6+2\Delta h/k)}=72.6\ \text{ms (Dirichlet)},\quad 60.3\ \text{ms }(h=1460)$$

Explicit stepping needs no linear solve and is trivially verifiable; at this grid the stability
limit is not restrictive. Output is on 241 log-spaced times from 1 ms to 600 s, and the first
interval is integrated exactly because $\theta$ is linear in $t$ before diffusion acts.

The steady state is a separate solve, because the slowest mode of a 30 mm block has time constant
$1/[\alpha(\pi/2d)^2]=2540$ s — a ten-minute transient is nowhere near it, and the steady state is
what sets the maximum continuous power. `heat_steady.m` solves

$$A\theta=b,\qquad A=-W\odot\nabla^2_h,\qquad b=W\odot q/k$$

by matrix-free preconditioned conjugate gradients. Multiplying by the control-volume weight $W$
(½ on a face, ¼ on an edge, ⅛ at a corner) is what makes $A$ **symmetric** — the bare ghost-node
Laplacian is not, because a boundary node owns half a cell. The preconditioner is the exact
diagonal, $W(6+2\Delta h/k)/\Delta^2$.

### 7.3 The energy audit

$$\frac{dU}{dt}=P_\text{in}-P_\text{loss},\qquad
U=\sum_i W_i\,\rho c\,\theta_i\,\Delta^{3},\qquad
P_\text{loss}=-k\Delta^{3}\sum_i W_i(\nabla^{2}_h\theta)_i$$

The second identity is exact for this discretisation: the weighted sum of the discrete Laplacian
telescopes to the boundary fluxes. It is therefore a real test of the boundary implementation
rather than a tautology, and it is evaluated at every output time.

One bookkeeping point worth stating: the power deposited **at** a node pinned to 37 °C never enters
the solution — the perfect heat sink removes it instantly. Both totals are reported.

---

## 8. Validation

`matlab/verify_heat_model.m`. Nine checks, in `results/full_heat_log.txt`.

**8.1 A three-dimensional closed form.** A source $q_0\sin\frac{\pi x}{L}\sin\frac{\pi y}{L}\sin\frac{\pi z}{L}$
in a box with all six faces Dirichlet is an eigenfunction of the Laplacian, so

$$\theta(\mathbf r,t)=\frac{q_0}{k\Lambda}\Big(1-e^{-\alpha\Lambda t}\Big)\sin\tfrac{\pi x}{L}\sin\tfrac{\pi y}{L}\sin\tfrac{\pi z}{L},
\qquad \Lambda=3\Big(\frac{\pi}{L}\Big)^{2}$$

exactly. This exercises the whole operator, all six boundaries and the time integration at once.

| $\Delta h$ | vs the continuum eigenvalue | vs the discrete eigenvalue | steady solver vs discrete |
|---|---|---|---|
| 1.00 mm | $1.45\times10^{-2}$ | $1.42\times10^{-2}$ | $7.6\times10^{-16}$ |
| 0.50 mm | $3.68\times10^{-3}$ | $3.59\times10^{-3}$ | $3.7\times10^{-15}$ |
| 0.25 mm | $9.27\times10^{-4}$ | $9.06\times10^{-4}$ | $1.4\times10^{-14}$ |

The continuum column falls by a factor 3.95 per halving — **second order**, as it must be. The
"discrete" column replaces $\Lambda$ by its exact discrete counterpart
$\Lambda_d=6(1-\cos\pi\Delta h/L)/\Delta h^2$, which removes the spatial error and leaves only
forward-Euler error; it happens to fall at the same rate because $\Delta t\propto\Delta h^2$ through
the stability limit. Cutting the time step at fixed grid separates them: at $\Delta h=1$ mm,
CFL 0.9 → 0.1 (a factor 9 in $\Delta t$) reduces it from $1.42\times10^{-2}$ to
$1.60\times10^{-3}$, a factor 8.9. **The steady solver reproduces the discrete solution to
machine precision.**

**8.2 The Robin boundary condition, exactly.** Uniform source, adiabatic sides and floor,
convective surface. The exact solution is a quadratic,
$\theta(s)=q_0d/h+(q_0/k)(ds-s^2/2)$, and the three-point stencil plus the ghost-node Robin form
reproduce a quadratic **exactly** — so any error here is a coding error, not truncation.

| | exact | computed |
|---|---|---|
| $\theta$ at the surface | 3.333 K | 3.333 K |
| $\theta$ at the floor | 70.000 K | 70.000 K |
| max relative error over the profile | | $1.7\times10^{-14}$ |
| $P_\text{in}$ vs $P_\text{out}$ through the surface | 0.100000 W | 0.100000 W, error $1.4\times10^{-14}$ |

**8.3 The adiabatic box.** With every boundary insulated, no heat can leave and $\theta$ must rise
uniformly as $q_0t/\rho c$. Computed 3.588517 K against exactly 3.588517 K, relative error
$7.4\times10^{-16}$; boundary loss identically zero; energy audit $4.7\times10^{-16}$. This is the
sharp test of the mirror ghosts and of the energy bookkeeping.

**8.4 The steady operator.** Symmetry $\langle Au,v\rangle=\langle u,Av\rangle$ to
$2.3\times10^{-16}$ relative on random vectors, and $\langle u,Au\rangle>0$ — so PCG is legitimate.
The steady solution agrees with the explicit transient run out to 20 000 s to **0.00 %**, and the
steady-state heat balance closes: $P_\text{in}$ = 110.2039 W against $P_\text{loss}$ = 110.2039 W, relative error
$1.8\times10^{-13}$. (The block receives 169.0309 W in total on that coarse grid; 58.8 W of it is
deposited at nodes pinned to 37 °C and is removed by the sink without ever entering the solution,
which is why the balance is taken over the free nodes.)

**8.5 Electromagnetic cross-validation.** `sar_volume.m` duplicates the field-model solver, so it is
checked against it rather than trusted.

| | this model | field model | difference |
|---|---|---|---|
| delivered power | 150.000 W | 150.000 W | — |
| absorbed in the gel | 149.586 W (99.72 %) | 149.584 W (99.72 %) | $1\times10^{-5}$ |
| antenna current scale $I_m$ | 4.0082 A | 4.0186 A | 0.26 % |
| $\sigma_\text{eff}$ | 2.72283 S/m | 2.72283 S/m | 0 |
| $q$ on $y=0$, 174 points, 1–29 mm deep | | | max 1.6 %, **rms 0.53 %** |
| `E_gel(6 mm, −5 mm)` power density | $5.0654\times10^{7}$ W/m³ | $5.0812\times10^{7}$ W/m³ | 0.31 % |

The 0.26 % in $I_m$ is the depth quadrature — Simpson on 121 nodes here against a midpoint rule on
120 cell centres there — and it propagates as $2\times0.26$ % into $q$, which accounts for most of
the 0.53 % rms. **The two models describe the same field.**

94.71 % of the power absorbed in the whole gel falls inside the 30 × 30 × 30 mm block.

**8.6 Grid convergence.** Two columns of this table behave completely differently, and that is the point.

| $\Delta h$ | nodes | peak, unrestricted | peak, feed excluded | $\Theta$(6, 0, −2 mm) | $P_\text{in}$ |
|---|---|---|---|---|---|
| 1.00 mm | 31³ | 896.65 K | 693.71 K | 506.637 K | 110.204 W |
| 0.50 mm | 61³ | 1228.82 K | 820.89 K | 532.706 K | 124.267 W |
| 0.25 mm | 121³ | **1437.14 K** | 865.02 K | **539.420 K** | 132.269 W |
| *with the feed source capped:* | | | | | |
| 1.00 mm | 31³ | 851.14 K | 687.82 K | 506.575 K | 110.073 W |
| 0.50 mm | 61³ | 771.08 K | 722.86 K | 531.900 K | 120.746 W |
| 0.25 mm | 121³ | 731.64 K | 701.87 K | 538.152 K | 124.635 W |

A probe 6 mm from the feed converges at **second order** (differences 26.1 then 6.7 K, ratio 3.9).
The unrestricted peak **diverges**, because it sits on the singularity of §6.2. Excluding the 2 mm
sphere tames it but does not cure it — the masked peak still creeps up, because it hugs the
exclusion boundary. Capping the source removes the divergence entirely and gives a peak stable to
about 5 % across a 4× change in grid spacing. This table is the evidence behind the treatment
adopted in §6.2, and behind the ±15 % assigned to it in §10I.

**8.7 Time-step convergence.** | CFL | $\Delta t$ | steps | $\Theta_\text{max}$(30 s) | energy-audit closure |
|---|---|---|---|---|
| 0.900 | 0.2612 s | 279 | 820.7132 K | $1.9\times10^{-3}$ |
| 0.450 | 0.1306 s | 388 | 820.1738 K | $1.4\times10^{-3}$ |
| 0.225 | 0.0653 s | 611 | 819.8769 K | $9.3\times10^{-4}$ |

Converged to **0.1 %** across a 4× change in time step, so the CFL 0.9 used in production is not a
limiting error. The residual in the energy audit is the trapezoid applied to $P_\text{loss}$ between
output times, not the solver: it falls when output times are added, not when the time step is cut.

**8.8 Linearity in power.** $3.7\,\theta(P)$ against $\theta(3.7P)$: relative difference 0 at 5 s
and $1.3\times10^{-16}$ at 60 s. Exact by construction — the check is that no power-dependent term
has crept into the implementation, since this identity is what `T_gel` relies on.

**8.9 Thermal dose.** `cem43.m` against hand calculation, exploiting the fact that each degree above 43 °C doubles
the dose and each degree below quarters it:

| held for | at | computed | exact |
|---|---|---|---|
| 10 min | 43 °C | 10.0000 | 10 |
| 10 min | 33 °C | $9.5367\times10^{-6}$ | $10\times0.25^{10}$ |
| 1 min | 53 °C | 1024.00 | $0.5^{-10}$ |
| 1 min | 60 °C | $1.3107\times10^{5}$ | $0.5^{-17}$ |

Relative error $\le 4\times10^{-15}$ throughout. The integrand spans 17 orders of magnitude over
this range, which is why each interval is integrated as an exponential rather than by trapezoid.

**8.10 The energy audit on the production run.** Over the full 600 s baseline, the closure of
$dU/dt=P_\text{in}-P_\text{loss}$ is $1.2\times10^{-3}$ — and that residual is dominated by the
trapezoid used on $P_\text{loss}$ between output times, not by the solver.

---

## 9. Results

All temperature rises below are $\Theta$, i.e. **per 150 W delivered = 200 W generator**. Divide by
200 to get kelvin per generator watt, or use the boxed formula of §0.

### 9.1 The peak temperature, and the fact that 200 W is impossible

| $t$ [s] | $\Theta_\text{max}$ [K] | $T_\text{max}$ at 200 W [°C] | hot spot $(x,\,y,\,\text{depth})$ mm |
|---|---|---|---|
| 0.01 | 1.53 | 38.5 | (7.7, 0, 0.25) |
| 0.10 | 13.77 | 50.8 | (7.5, 0, 0.25) |
| 0.51 | 48.22 | 85.2 | (6.5, 0, 0.50) |
| 1.00 | 78.82 | 115.8 | (6.5, 0, 0.50) |
| 2.05 | 152.81 | 189.8 | (0.5, −1.2, 1.50) |
| 5.00 | 356.49 | 393.5 | (0.25, 0, 2.00) |
| 9.75 | 557.56 | 594.6 | (0.25, 0, 2.00) |
| 29.7 | 862.56 | 899.6 | (0.25, 0, 2.00) |
| 61.2 | 1014.75 | 1051.7 | (0.25, 0, 2.00) |
| 119 | 1119.64 | 1156.6 | (0.25, 0, 2.00) |
| 308 | 1207.69 | 1244.7 | (0.25, 0, 2.00) |
| 600 | 1231.94 | 1268.9 | (0.25, 0, 2.00) |
| ∞ | **1238.78** | 1275.8 | (0.25, 0, 2.00) |

Read the third column as a *reductio*: **at the generator setting the field model was normalised to,
the gel boils in under a second.** The phantom is not a load that can absorb 150 W.

The peak reaches 82 % of its steady value by 60 s and 90 % by 120 s. The near-surface region
equilibrates in seconds because it is only millimetres from a fixed-temperature boundary; the
remaining slow creep is the rest of the 30 mm block warming up.

### 9.2 Where the hot spot is, and where it moves

Two regimes, and the boundary between them is the feed artefact of §6.2:

* **$t \lesssim 1$ s** — the maximum sits under the *middle of the radiating element*, at
  $x=6.5$–$7.7$ mm and only 0.25–0.5 mm deep. This is the physically trustworthy answer: it is
  6–8 mm from the feed, where §6.2 shows the artefact changes nothing.
* **$t \gtrsim 2$ s** — the maximum migrates to $(0.25,\,0,\,-2.0)$ mm, i.e. onto the surface of the
  2 mm feed-exclusion sphere. Heat from the (unphysical) feed charge has by then diffused out far
  enough to beat the element's own contribution. **This part of the answer is set by the artefact,
  and is reported for completeness rather than believed.**

Two consequences worth acting on. First, the feed-regularised model (§10I) moves the late-time peak
to $(3.0\text{–}3.5,\,0,\,-4.0)$ mm — under the element, at a sensible depth — and lowers it by only
3–6 % depending on the grid, so the *numbers* are robust even though the *location* is not. Second,
a real applicator does have a feed-junction hot spot; this model cannot size it, so **instrument the
feed on the bench.**

The occasional off-plane report ($y=-1.2$ mm at 2.05 s) is not physics: the deposition is symmetric
in $y$, and that entry is simply the hottest point on the exclusion sphere, which cuts the $y=0$
plane there.

### 9.3 The cooled surface, and why the catheter does not make a hot spot

On the vertical through the middle of the element, $x=6$ mm:

* surface ($z=0$): $\Theta = 0$ by construction — the blood holds it at 37 °C;
* peak: $\Theta = 1101$ K at **4.25 mm depth**;
* the profile is monotone from 0 to 4.25 mm and falls slowly beyond.

The deposition is *largest* at the surface and the temperature is *smallest* there. That inversion
is the whole physical content of the catheter reminder in the brief: the applicator sits in the
coolant, so the most intensely irradiated gel is also the best cooled, and the thermal maximum is
driven a few millimetres inward. It is the mechanism irrigated microwave ablation uses on purpose.

The catheter's own footprint changes nothing in the baseline, and it cannot: a cooled body resting
on a plane already held at 37 °C leaves it at 37 °C. But §10C shows that this is only harmless
*because the catheter really is cooled*. Give the footprint the coefficient the flow actually
provides ($h_\text{cath}=857$ W m⁻²K⁻¹) and the penalty is 8 % in allowed power. Make the same
1.5 mm strip perfectly insulating and $\Theta$ rises by **62 %**, because that strip lies exactly
over the most intensely deposited gel in the whole block.

So the reminder in the brief is load-bearing, not a detail. **If the catheter did not shed heat to
the blood, it would be the most important surface feature in the model** — more important than the
choice of gel conductivity, and comparable to the flow geometry itself. Two practical corollaries:
keep the catheter wall thin and thermally conductive where it contacts the gel, and do not let an
air gap or a dry contact form under it.

### 9.4 Where the power goes

Deposited in the 30 mm block, at 150 W delivered: **142.81 W**, of which 10.54 W lands on nodes held
at 37 °C and is removed instantly by the sink, leaving **132.27 W** to heat the gel.

| $t$ [s] | stored [J] | leaving [W] | as % of input |
|---|---|---|---|
| 1.0 | 124 | 12.6 | 9.5 % |
| 9.8 | 1 051 | 34.9 | 26.4 % |
| 61.2 | 5 054 | 67.7 | 51.2 % |
| 600 | 17 448 | 126.2 | 95.4 % |
| ∞ | — | 132.27 | 100.00 % |

The last row is the steady-state closure check: everything deposited leaves. Even at that extreme
the coolant bulk temperature rises only **0.38 K**, which confirms the brief's instruction to ignore
it.

### 9.5 The operating chart — the number you actually need

Generator forward power at which the **peak** temperature just reaches a limit, at that exposure:

| exposure | 43 °C | 50 °C | 60 °C | 100 °C |
|---|---|---|---|---|
| 1 s | 15.22 W | 32.99 W | 58.36 W | 159.9 W |
| 5 s | 3.37 | 7.29 | 12.90 | 35.34 |
| 10 s | 2.15 | 4.66 | 8.25 | 22.60 |
| 30 s | 1.39 | 3.01 | 5.33 | 14.61 |
| 61 s | 1.18 | 2.56 | 4.53 | **12.42** |
| 119 s | 1.07 | 2.32 | 4.11 | 11.25 |
| 308 s | 0.99 | 2.15 | 3.81 | 10.43 |
| 600 s | 0.97 | 2.11 | 3.73 | 10.23 |
| indefinitely | **0.97** | 2.10 | 3.71 | **10.17** |

(The exposures are the nearest stored output times; `P_for_T` interpolates to whatever time you ask
for, and returns 12.46 W at exactly 60 s.)

This is `figures/fig4_safe_power.png`, and it is the single most useful output of the model.
Two readings of it:

* **You have about 12 W to play with for a one-minute exposure, and about 10 W for ever.** Not 200.
* The curve flattens after roughly 100 s, because by then the deposition and the conduction losses
  are nearly in balance. Longer exposures do not cost much more power headroom; they cost
  thermal dose (§9.6).

**Caveat on the "indefinitely" row.** The steady state of this model is the steady state of a 30 mm
block whose outer walls are held at 37 °C — which is the boundary condition the brief specified. By
600 s the thermal diffusion length is 18.6 mm and those walls are doing real work, so the last two
rows depend on the size of the region of interest in a way the short-exposure rows do not (§10E, J).
Read them as "a 30 mm neighbourhood in thermal equilibrium with tissue at 37 °C", not as a universal
constant.

### 9.6 Thermal dose and the size of the lesion

`cem43.m` converts the temperature history into cumulative equivalent minutes at 43 °C
(Sapareto & Dewey [T6]). A Phytagel phantom has no biology, so this is a **translation of the
phantom result into the tissue it stands for**, not a property of the gel. The conventional
ablation threshold is CEM43 = 240 min.

Dimensions of the dosed zone after a 60 s exposure, on the $y=0$ plane with the feed sphere masked
out — so these are the zone the **radiating element** produces:

| generator power | peak $T$ | CEM43 ≥ 240 min | CEM43 ≥ 60 min | $T\ge43$ °C |
|---|---|---|---|---|
| 6.21 W | 68.5 °C | 5.75 mm deep, 10.7 mm long | 6.50 mm, 11.7 mm | 10.50 mm deep |
| 9.31 W | 84.3 °C | 7.75 mm deep, 13.5 mm long | 8.25 mm, 14.2 mm | 12.25 mm deep |
| 12.42 W | 100.0 °C | **9.00 mm deep, 15.0 mm long** | 9.50 mm, 15.5 mm | 13.50 mm deep |

Three things follow.

**The dose contour is sharp; the temperature contour is not.** Between CEM43 = 60 and 240 min there
is only 0.5 mm of depth, because the dose is exponential in temperature. The 43 °C isotherm, by
contrast, reaches 4.5 mm deeper than the ablation contour. If you are looking for a lesion boundary
on the bench, the 240 min contour is the one that will look like an edge.

**The heated zone is longer than the antenna.** At 12.42 W it runs from $x=-3.5$ to $+11.5$ mm — it
extends about 3.5 mm *behind* the feed and stops just short of the tip. That backward reach is
partly the feed artefact of §6.2 and partly real; the field model's own limitation 3 notes that it
*under*-predicts heating behind the feed because coax sheath currents are absent. **Expect more
backward heating on the bench than this model shows, not less.**

**Depth scales weakly with power.** Doubling the power from 6.2 to 12.4 W deepens the ablation
contour only from 5.75 to 9.00 mm. Depth is set by the deposition profile and by conduction to the
blood, and neither of those changes with power. On the vertical through the middle of the element,
at 60 s:

| generator power | peak at $x$ = 6 mm | 43 °C reaches | 50 °C reaches |
|---|---|---|---|
| 5 W | 54.8 °C | 9.52 mm | 5.70 mm |
| 10 W | 72.7 °C | 12.56 mm | 9.16 mm |
| 12.42 W | 81.3 °C | 13.56 mm | 10.12 mm |
| 20 W | 108.4 °C | 15.86 mm | 12.20 mm |
| 50 W | 215.5 °C | 20.22 mm | 16.54 mm |

Eight times the power buys about twice the 43 °C depth, at the cost of boiling everything above it.
**If you want a deeper lesion, you need a longer exposure or a different antenna, not more watts.**

### 9.7 Figures

| file | content |
|---|---|
| `fig1_rise_per_watt.png` | $\Theta/P_\text{gen}$ on the $y=0$ plane at 10 s and at steady state — the power-independent core result |
| `fig2_T_maps.png` | temperature maps at 1, 10, 60 and 600 s at 12.42 W |
| `fig3_Tmax_vs_time.png` | peak temperature vs time for 1–200 W generator |
| `fig4_safe_power.png` | **the safe-operating chart** — allowed power vs exposure time, four limits |
| `fig5_depth_profiles.png` | depth profiles, and the sub-surface peak in close-up |
| `fig6_thermal_dose.png` | CEM43 map and lesion dimensions vs power |
| `fig7_energy.png` | the energy audit — where the deposited power goes with time |
| `fig8_pulsed.png` | pulsed vs CW at the same mean power |
| `fig9_sensitivity.png` | the parameter study |

---

## 10. Sensitivity analysis

`matlab/run_heat_sensitivity.m` → `results/heat_sensitivity_log.txt`, `figures/fig9_sensitivity.png`.

Every case is run on the same 0.5 mm grid to 120 s plus its exact steady state, so the rows are
comparable with each other. That grid reads the baseline $P_{100}$ 5 % high relative to the 0.25 mm
production grid (13.03 W against 12.42 W); the offset is common to every row and does not affect
the comparisons.

$\Theta$ is in kelvin per 150 W delivered. $P_{100}$ is the generator forward power at which the
peak reaches 100 °C after 60 s — the number a bench operator actually needs.

| case | $\Theta$(1 s) | $\Theta$(10 s) | $\Theta$(120 s) | $\Theta_\text{ss}$ | hot spot (x, y, depth) mm | $P_{100}$ [W] |
|---|---|---|---|---|---|---|
| **A baseline (surface pinned at 37 °C)** | **73.8** | **526.9** | **1076.8** | **1195.8** | (0.5, 0, 2.00) | **13.03** |
| B Robin, $h$ = 5000 (near-Dirichlet) | 88.7 | 598.0 | 1228.7 | 1359.4 | (0.5, 0, 2.00) | 11.39 |
| B Robin, $h$ = 1460 (flat plate, laminar) | 127.1 | 711.1 | 1521.6 | 1681.8 | (0.5, 0, 2.00) | 9.20 |
| B Robin, $h$ = 673 (duct, Dittus–Boelter) | 159.7 | 805.0 | 1869.4 | 2081.5 | (0.5, 0, 2.00) | 7.56 |
| B Robin, $h$ = 300 | 179.9 | 890.1 | 2362.5 | 2698.6 | (0.5, 0, 2.00) | 6.13 |
| B Robin, $h$ = 103 (laminar duct, worst case) | 192.1 | 983.3 | 3031.8 | 3720.1 | (0.5, 0, 2.00) | 5.01 |
| C catheter cooled, $h_\text{cath}$ = 857 (with $h$ = 1460) | 150.8 | 762.4 | 1648.6 | 1815.9 | (0.5, 0, 2.00) | 8.47 |
| C catheter **insulating** (with $h$ = 1460) | 198.3 | 924.5 | 2168.8 | 2370.6 | (0.5, 0, 2.00) | 6.43 |
| C catheter **insulating** (pinned surface) | 196.0 | 848.5 | 1742.0 | 1884.7 | (0.5, 0, 2.00) | 7.88 |
| D floor held at 37 °C (water-cooled dish) | 73.8 | 526.9 | 1076.8 | 1195.0 | (0.5, 0, 2.00) | 13.03 |
| E ROI side walls adiabatic | 73.8 | 526.9 | 1077.6 | 1532.3 | (3.5, 0, 6.00) | 13.03 |
| F $k$ = 0.50 W m⁻¹K⁻¹ | 77.8 | 568.7 | 1262.0 | 1435.0 | (0.5, 0, 2.00) | 11.27 |
| F $k$ = 0.55 | 75.8 | 547.2 | 1162.0 | 1304.6 | (0.5, 0, 2.00) | 12.15 |
| F $k$ = 0.628 (pure water at 37 °C) | 72.8 | 516.0 | 1034.4 | 1142.5 | (0.5, 0, 2.00) | 13.52 |
| G $c$ = 4000 J kg⁻¹K⁻¹ | 76.1 | 539.6 | 1082.5 | 1195.8 | (0.5, 0, 2.00) | 12.92 |
| G muscle thermal properties | 85.0 | 613.8 | 1304.0 | 1464.3 | (0.5, 0, 2.00) | 10.83 |
| H gel permittivity at 37 °C, not 25 °C | 76.6 | 554.9 | 1118.1 | 1235.7 | (0.5, 0, 2.00) | 12.50 |
| I feed exclusion radius 1 mm | 263.6 | 862.0 | 1159.5 | 1246.3 | (0.0, 0, 1.50) | 11.41 |
| I feed exclusion radius 3 mm | 73.8 | 326.1 | 984.8 | 1173.6 | (3.0, 0, 4.00) | 15.20 |
| I feed source regularised (capped) | 73.8 | 453.5 | 971.2 | 1158.5 | (3.0, 0, 4.00) | 14.59 |
| J ROI 20 mm wide | 73.8 | 526.8 | 1043.6 | 1099.0 | (0.5, 0, 2.00) | 13.20 |
| J ROI 40 mm wide | 73.8 | 526.9 | 1077.3 | 1220.4 | (0.5, 0, 2.00) | 13.03 |

### The ranking

| rank | parameter | effect on $P_{100}$ | what to do about it |
|---|---|---|---|
| 1 | **surface film coefficient $h$** | **2.6×** (5.0 → 13.0 W) | measure or control the flow geometry over the phantom; §10B |
| 2 | **the catheter, if it insulates** | 1.65× (7.9 → 13.0 W) | it should not, but a thick or dry-contacting wall would; §10C |
| 3 | how the feed singularity is handled | ±15 % | intrinsic to the field model's current distribution; §6.2 |
| 4 | gel thermal conductivity | ±9 % | measure it, or accept the range |
| 5 | temperature-dependent properties | +12–15 % on $\Theta$ | breaks exact linearity; §10K |
| 6 | everything else | < 5 % | ignore |

Note what is **not** in the list: the delivery efficiency $\eta$. It multiplies the answer exactly
and is therefore not a modelling uncertainty at all — it is a measurement you should make.
$P_{100}$ scales as $1/\eta$: at a perfect match ($\eta=1$) the 13.03 W becomes 9.77 W.

### A. The baseline

$\Theta$(1 s) = 73.8 K, $\Theta_\text{ss}$ = 1196 K per 150 W delivered; $P_{100}$ = 13.03 W. Every
comparison below is against this row.

### B. The surface film coefficient — the dominant sensitivity

Going from the specified Dirichlet condition to the most pessimistic defensible film coefficient
costs a factor **2.6** in allowed power. This is the thermal counterpart of the standoff gap in the
field model, and it has the same character: it depends on a geometric detail nobody has measured.

The internal ordering is exactly as it should be — $h=5000$ is already within 13 % of the Dirichlet
limit, confirming that "held at 37 °C" is the $h\to\infty$ end of a continuum, and that the
Dirichlet condition is an **optimistic** bound, not a neutral one.

At the baseline $h=1460$ the full production run (`results/heat_model_robin.mat`) gives
$P_{100}$ = 8.63 W at 60 s and 7.08 W indefinitely, against 12.42 and 10.17 W for the pinned
surface. **The convective surface also stops the surface being at 37 °C**: at $x=6$ mm the steady
surface rise is 571 K per 150 W delivered instead of zero, which is the film resistance the
Dirichlet condition throws away.

**Recommendation: design the flow cell so that the gel patch sits at the leading edge of a fast,
thin channel** — that maximises $h$, pushes the answer towards the Dirichlet limit, and makes the
model you were given the right one.

### C. The catheter

This is the case the brief specifically flagged, so all three treatments are run.

* **Cooled** (the physically expected case, $h_\text{cath}=857$ from Churchill–Bernstein in series
  with the wall): $P_{100}$ falls from 9.20 to 8.47 W, i.e. **8 %**. The catheter is very nearly
  thermally transparent, because its own convective coefficient (5940 W m⁻²K⁻¹, high because the
  diameter is small) is much better than the flat surface's and only the polymer wall degrades it.
* **Insulating** (a catheter that blocks the coolant and conducts nothing): $P_{100}$ falls to
  6.43 W with a convective surface, or 7.88 W with a pinned one — **a 30–40 % penalty.**
* In every case the strip is 1.5 mm wide, which is smaller than the thermal diffusion length after
  the first second, so the gel underneath is cooled sideways even when the strip itself is not.
  That is why even the pessimistic bracket costs less than a factor of two.

**The catheter never becomes a heat source and never creates a surface hot spot.** The maximum stays
2–4 mm deep in all three treatments.

### D and E. The boundaries of the region of interest

The floor makes **no difference at all** to three decimal places — the same conclusion the field
model reached for $|\mathbf E|$ below 20 mm, for the same reason: almost nothing gets that far.

The side walls make no difference to 120 s ($\Theta$ 1077.6 adiabatic vs 1076.8 pinned) and a large
one at steady state (+28 %, and the hot spot moves out to (3.5, 0, −6.0) mm). Combined with case J
— 20 / 30 / 40 mm wide regions give $\Theta_\text{ss}$ of 1099 / 1196 / 1220 K — this says:

> **Everything up to a couple of minutes is independent of the region of interest. The steady state
> is not, and should be read as "a 30 mm neighbourhood surrounded by tissue at 37 °C" — which is
> exactly the homeostasis assumption in the brief.**

### F and G. Thermal properties

$k$ across its plausible range 0.50–0.628 W m⁻¹K⁻¹ moves $P_{100}$ by ±9 %. Specific heat is
negligible (0.8 % for a 4 % change) because at these timescales the answer is governed by conduction
away, not by heat capacity. Substituting muscle properties wholesale costs 17 % — worth knowing if
you want the phantom result to stand for tissue, and a reminder that the phantom is more conductive
and higher heat capacity than what it mimics.

### H. Permittivity temperature

Evaluating the gel permittivity at 37 °C rather than the field model's 25 °C changes $P_{100}$ by
−4 %. Small, and it is kept at 25 °C in the baseline only so that the cross-validation of §8.5
against `E_gel` is exact.

### I. The feed treatment

$P_{100}$ = 11.41 / 13.03 / 15.20 W for exclusion radii of 1 / 2 / 3 mm, and 14.59 W with the
source regularised. So **the feed artefact is worth about ±15 %** on the headline number, and the
2 mm baseline is the conservative choice among them. Note also that the regularised case puts the
hot spot at (3.0, 0, −4.0) mm — under the radiating element, which is where physics says it should
be — and gives essentially the same answer as a 3 mm exclusion. That agreement is the reason the
±15 % can be quoted as an uncertainty rather than a failure.

### K. Temperature-dependent properties — how far the exact power scaling can be trusted

The Stogryn model gives, for this gel at 2.45 GHz, $\sigma_\text{eff}=2.782$ S/m at 37 °C rising to
3.179 at 60 °C (+0.62 %/K) and $\varepsilon'$ falling at −0.38 %/K. Re-running the transient with
$q$ scaled pointwise by $\sigma(T)/\sigma(37)$, holding $|\tilde{\mathbf E}|^2$ fixed:

| generator power | $T_\text{max}$(60 s), linear | with $\sigma(T)$ | difference in the rise |
|---|---|---|---|
| 5 W | 61.2 °C | 64.0 °C | **+11.7 %** |
| 10 W | 85.4 °C | 92.6 °C | **+14.9 %** |
| 20 W | 133.7 °C | 135.5 °C | +1.8 % |

So the exact power scaling of §0 is good to about **10–15 % once the gel is tens of kelvin above the
coolant, and it errs low** — the real gel absorbs more as it warms. (The 20 W row is smaller only
because $\sigma$ turns over above 80 °C in the Stogryn model, which is past the point where the
phantom has stopped being a phantom anyway.)

This is a **bound, not a coupled solution**: the field distribution is held fixed, so it captures
the change in absorbed power but not the redistribution caused by $\varepsilon'$ falling. A properly
coupled electromagnetic–thermal iteration is the honest next step if you need better than 15 %.

---

## 11. The deliverable functions

### 11.1 `T_gel` — power and time in, temperature out

```matlab
T = T_gel(P, t)                  % peak temperature in the block          [°C]
T = T_gel(P, t, x, z)            % at a point on the y = 0 plane
T = T_gel(P, t, x, y, z)         % anywhere in the block
[T, S] = T_gel(...)              % + rise, CEM43, hot-spot location, validity flags
T = T_gel(P, t, ..., 'delivered')% P is delivered power, not generator forward power
T = T_gel(P, Inf, ...)           % steady state
```

Arguments are SI — metres and seconds — and $z$ is **negative** inside the gel, exactly as in
`E_gel`. `P` is by default the **generator forward power**, because that is the number on the front
panel; the gel sees $\eta P$.

Verified output (`results/full_heat_log.txt`):

| $P_\text{gen}$ | $t$ | peak [°C] | at (6, 0, −2 mm) | at (6, 0, −10 mm) |
|---|---|---|---|---|
| 5 W | 10 s | 51.13 | 44.39 | 37.86 |
| 5 W | 60 s | 62.27 | 54.22 | 42.27 |
| 5 W | ∞ | 67.97 | 60.50 | 56.27 |
| 10 W | 10 s | 65.25 | 51.78 | 38.72 |
| 10 W | 60 s | 87.55 | 71.43 | 47.54 |
| 10 W | ∞ | 98.94 | 84.00 | 75.54 |
| 20 W | 60 s | 138.09 | 105.87 | 58.09 |
| 50 W | 60 s | 289.74 | 209.17 | 89.72 |

and off the symmetry plane, at 10 W and 60 s, $T(6\ \text{mm},\,y,\,-2\ \text{mm})$ = 71.36, 69.63,
65.52, 52.12 °C at $y$ = 0, 1, 2, 5 mm. The temperature *rise* is therefore half its centre value at
$|y|\approx4.5$ mm: the heated zone is about **9 mm wide transversely** after one minute, against
the 15 mm it spans along the antenna.

**Two behaviours worth knowing before you use it.** The peak returned by `T_gel(P,t)` excludes the
2 mm feed sphere (§6.2) — point queries do not. And `S.flags` reports every validity boundary the
result has crossed; at 12 W for 60 s it returns three, including *"T exceeds 60 °C: permittivity and
conductivity drift ~10 %, power scaling no longer exact"*.

### 11.2 Pulsed and ramped drive, without another solve

The system is linear and time-invariant, so the response to any generator history $P(t)$ is a
Duhamel superposition of step responses:

$$\theta(\mathbf r,t)=\frac{1}{P_\text{ref}}\int_0^{t}P_\text{del}(\tau)\,
\frac{\partial\Theta}{\partial t}(\mathbf r,\,t-\tau)\,d\tau$$

Pass a function handle and `T_gel` evaluates it by discrete superposition, which is exact for a
piecewise-constant power history:

```matlab
Pw = @(tt) 41.55*(mod(tt,1) < 0.3);        % 30 % duty at 1 Hz, same mean as 12.46 W CW
T  = T_gel(Pw, linspace(0,60,4000));
```

**The result is the useful one: duty cycling does not buy headroom, and slow duty cycling is
actively dangerous.** All of the following deliver the same 12.42 W mean over 60 s:

| drive | peak power | peak temperature after 60 s |
|---|---|---|
| CW | 12.42 W | **99.8 °C** |
| 50 % duty, 1 Hz | 24.83 W | 100.6 °C |
| 50 % duty, 0.1 Hz | 24.83 W | 112.9 °C |
| 50 % duty, 0.02 Hz | 24.83 W | 138.9 °C |
| 30 % duty, 1 Hz | 41.39 W | 101.9 °C |
| 30 % duty, 0.1 Hz | 41.39 W | 117.1 °C |
| 30 % duty, 0.02 Hz | 41.39 W | **178.5 °C** |
| 10 % duty, 1 Hz | 124.17 W | 102.1 °C |
| 10 % duty, 0.1 Hz | 124.17 W | 118.8 °C |
| 10 % duty, 0.02 Hz | 124.17 W | **269.7 °C** |

The pattern is set by one number: the thermal time constant of the hot region, which is of order
seconds. A 1 Hz cycle is far faster than that, so the gel simply averages it and the peak lands
within 2 °C of CW. A 0.02 Hz cycle (50 s period) is far slower, so the gel follows each pulse and
the peak is set by the **peak power**, not the mean — at 10 % duty that is 124 W, and the answer is
the 200 W disaster of §9.1 in miniature.

**Rule of thumb: pulse at 1 Hz or faster, or do not pulse.** If you want to lower the peak at fixed
energy, lower the power and lengthen the exposure — the operating chart of §9.5 tells you how.
`figures/fig8_pulsed.png` shows the waveforms.

### 11.3 The inverses

```matlab
P_for_T(100, 60)        % = 12.46 W  : the power that just boils it in a minute
P_for_T(100, Inf)       % = 10.17 W  : the power it can take for ever
time_to_T(20, 100)      % = 12.44 s  : how long 20 W has before it boils
time_to_T(10, 100)      % = Inf      : 10 W never reaches 100 °C
```

Both accept vectors and broadcast them, so the whole safe-operating chart is one line:

```matlab
te = logspace(-1, 2, 60);  plot(te, P_for_T([43 100], te));
```

`P_for_T` is not a search — it inverts the exact linear relation of §0 algebraically.

### 11.4 `heat_scenario` — one bench run, drawn and tabulated

For validating the model against a measurement, `T_gel` is the wrong shape of tool: you want to type
in what the generator was set to and get back the figures and a table you can lay a thermocouple
trace on. That is `heat_scenario`:

```matlab
heat_scenario(12, [],  60)                      % 12 W continuous for 60 s
heat_scenario(40,  1,  60, 'duty', 0.3)         % 40 W, 30 % duty at 1 Hz, for 60 s
heat_scenario(12, [],  60, 'model', 'robin')    % the same, convective surface
heat_scenario(12, [],  60, 'probes', [6e-3 0 -2e-3; 6e-3 0 -6e-3])
heat_scenario(12, [],  60, 'measured', log)     % overlay a bench record
```

The three positional arguments are **power, pulse period, stimulation time**, *in that order and
before any options*. Pass `[]` for the period to get continuous wave. A period without a duty cycle
defaults to 50 % and says so.

```matlab
heat_scenario(20, [], 'model','robin', 360)     % WRONG -- 'model' lands in t_stim
heat_scenario(20, [], 360, 'model','robin')     % right
```

That is the easy mistake to make, so it is checked: the first form now reports which argument is
wrong and prints the corrected call, rather than failing later with an unrelated message. Option
names are checked against the known list too, so a typo is caught instead of silently ignored.

Each call writes three files:

| file | what |
|---|---|
| `figures/scenario_<name>_traces.png` | the drive waveform, the peak and every probe against time (with the cool-down after the stimulation ends), depth profiles, and a summary block |
| `figures/scenario_<name>_maps.png` | temperature on the $y=0$ plane at five times through the run, plus the thermal dose over the whole run |
| `results/scenario_<name>.csv` | time, drive power, peak and each probe trace — the file to load next to a measurement |

and prints a console summary: peak and where it is, the temperature at every probe at the end of
stimulation and its maximum, the dosed zone, and every validity flag the run has crossed.

**No new solve happens.** The response to any piecewise-constant power history is a superposition of
the stored step response,

$$\theta(\mathbf r,t)=\frac{1}{P_\text{ref}}\sum_j \Delta P_j\,\Theta(\mathbf r,\,t-\tau_j)$$

where the drive steps by $\Delta P_j$ at $\tau_j$. Continuous wave is exactly two terms — on at
$t=0$, off at $t_\text{stim}$ — so the **cool-down after the exposure is exact and free**. A pulse
train is the same sum evaluated as a discrete convolution. A run takes under a second.

Two grids are used, for two different jobs: the maps live on a coarse time grid because a
121 × 121 plane at many times is a large dense object, and the traces live on a grid fine enough to
resolve each individual pulse (20 samples inside the shortest on- or off-phase). Getting that second
grid wrong smears the pulse peaks — with 5 samples per pulse the reported peak came out 4 % low.

**Three things to know before comparing against a measurement.**

1. **The baseline is a lower bound.** The surface is held at exactly 37 °C, as specified. Run
   `'model','robin'` as well and expect the measurement to fall between the two — at 12 W for 60 s
   they bracket the peak at 97.7 and 124.3 °C, and a probe at (6, 0, −2) mm at 78.3 and 88.3 °C.
   If the measurement falls outside that bracket, the flow geometry is not what §3 assumes and that
   is the first thing to check.
2. **Do not put a probe near the feed.** Inside 2 mm of $(0,0,0)$ the field model is singular
   (§6.2). Those cells are blanked grey on the maps, and the reported peak excludes them. Put
   probes under the middle of the element, around $x=6$ mm.
3. **Check the flags.** Above 60 °C the constant-property assumption starts to cost 10–15 %, and it
   under-predicts (§10K); above ~50 °C the phantom itself may be softening (§12.5). A validation run
   that stays below 50 °C tests the model rather than the gel — at 12 W a 10 s exposure, or 5 W for
   60 s, will do that.

Seven worked examples ship with the model, in `Heatmodel/results/`: `scenario_P12W_cw_60s` and its
`_robin` twin (the same drive, convective surface), `scenario_P40W_1Hz_30pc_60s` (the same mean
power, pulsed), `scenario_P50W_cw_60s` and `scenario_P50W_0.0333333Hz_50pc_60s` (the slow-pulsing
case of §9.7), and `scenario_P20W_cw_360s` with its `_robin` twin (a six-minute exposure, where the
two surface models are furthest apart).

### 11.5 `heat_gui` — sweeping the physics, not just the drive

Everything above takes the physics as fixed. `T_gel` and `heat_scenario` are fast because they reuse
one stored $\Theta$, and that is only legitimate while the gel, the coolant and the catheter are the
ones this document assumed. Change the permittivity and $\Theta$ itself is wrong.

`heat_gui` exists for the other question — *what if the gel is not what we assumed?* It re-solves
the whole chain, electromagnetics included, in a few seconds, and draws the two maps that answer it:

```matlab
heat_gui                              % the window
R = gui_model('sigma',1.8,'T_water',30,'Qflow_Lmin',20,'Pgen',25);   % headless
```

| what you can change | what it reaches |
|---|---|
| gel permittivity $\varepsilon'$, conductivity $\sigma$ | the source term $q$ — a new EM solve |
| gel temperature | the initial condition, **and** $\varepsilon',\sigma$ via Stogryn when the box is ticked |
| water temperature | the sink temperature, and the water properties that set $h$ |
| water flow rate | $h$, through `convection_h` |
| catheter length | the footprint of the second cooled surface |
| generator power, exposure time, grid | the drive and the numerics |

**How it is used.** Pick one of five starting points from the top drop-down — *bench baseline*,
*gentle test*, *high power*, *cold gel from the bench*, *write-up baseline* — then move whatever you
like. The maps redraw by themselves: about a second when only the power or the time changed, because
the electromagnetic solve is cached and only the conduction is redone, and about five seconds when
the gel itself changed. Above the maps is the answer in one sentence — peak temperature, the
temperature at the probe point, and whether that is past 43 °C, past 60 °C, or past boiling. The
numerical detail sits underneath, and the grid and colour-scale controls stay hidden behind
*advanced* until they are wanted. Choosing the 0.25 mm grid turns the automatic redraw off by
itself, because a solve there takes 90 s and a window that appears to hang is worse than a button.

*Selecting the write-up baseline preset reproduces §9.5 to the last degree: 97.6 °C against the
97.7 °C computed on the quarter-hour grid.*

**Two defaults had to change, and it matters.** This document's baseline is a Dirichlet surface with
a transparent catheter. Under it, *the flow rate and the catheter length do nothing whatsoever* — a
perfect heat sink does not care how fast the coolant moves, and a cooled body resting on a surface
already pinned at 37 °C changes nothing (§3.2, §4). Two of the six parameters would be dead
controls. The GUI therefore opens on the physical pair instead — convective surface, cooled catheter
— and keeps both baselines in the drop-downs. The Dirichlet answer is still the lower bound.

**The tick box is a trap worth understanding.** Linked, the gel temperature moves $\varepsilon'$ and
$\sigma$ *and* the initial condition together, so a sweep of it is not a clean experiment. Two tests
in `verify_gui_model` were wrong on first writing for exactly this reason. Untick it to vary the
temperature and the permittivity independently — which is also what you want once the gel has been
measured with a coaxial probe rather than inferred (§4, §11).

**A gel that does not start at 37 °C.** A phantom taken from the bench sits at room temperature, so
$\theta(\mathbf r,0)\neq0$ and the exact power scaling of §11.1 no longer holds: the answer is a
driven part, linear in $P$, plus a relaxation towards the coolant that does not scale with $P$ at
all. `heat_solve` grew an `opt.theta0` for this; it defaults to zero, so nothing else in the model
changed. Test 6 of `verify_gui_model` confirms the two parts superpose to single precision.

**What it costs.** The GUI runs at 0.5 mm with a 512-point FFT window instead of 0.25 mm and 1024,
and only as far in time as the exposure lasts. `verify_gui_model` measures the price against the
stored quarter-hour baseline rather than asserting it:

| check | result |
|---|---|
| against `heat_model.mat` (Dirichlet), bench point at 60 s | −0.44 K, **1.1 %** of the rise |
| against `heat_model_robin.mat`, bench point at 60 s | −1.98 K, **3.8 %** of the rise |
| grid convergence 1.0 → 0.5 → 0.25 mm, bench point | 0.5 % over the last refinement |
| FFT window 512 vs 1024 | 0.06 % on the bench point, 0.8 % on the SAR peak |
| exact scaling in delivered power | to machine precision |
| initial-condition superposition | 3 × 10⁻⁸ (single-precision storage floor) |
| energy balance at the end of the run | 4 × 10⁻⁴, and first order in $\Delta t$ |

A few percent, for a solve that takes 5 s instead of 15 min. The **peak** is deliberately not
asserted on: it sits on the feed-exclusion sphere (§6.2), where refining the grid keeps finding a
larger value, so the GUI reports the bench point (6, 0, −2) mm beside it and that is the number to
compare with a thermocouple.

### 11.6 What is stored

`results/heat_model.mat` (68 MB, `-v7.3`) holds $\Theta$ on the $y=0$ plane at full 0.25 mm
resolution and all 241 output times, the 3-D history decimated to 0.5 mm at 61 times, the
full-resolution steady state, the peak history and its location, the deposition $q$, the energy
audit, and the configuration struct. `results/heat_model_robin.mat` is the same for the convective
surface.

---

## 12. Limitations — read before quoting numbers

1. **The perfect heat sink is optimistic, and it is the largest uncertainty.** The brief asked for
   a 37 °C surface, and that is the baseline; but it is the $h\to\infty$ end of a range that runs
   down to $h\approx100$ W m⁻²K⁻¹ for the same 5 L/min through a different channel. The allowed
   power spans a factor 2.6 across that range (§10B). **Fix the flow geometry, and the model
   becomes quantitative.** The same caveat applies to the catheter footprint: it is harmless in the
   model only because the coolant reaches it. An insulating strip in the same place costs another
   30–40 % (§10C), which makes the thermal contact under the catheter the second thing to get right.

2. **The feed singularity.** The field model's prescribed current is truncated at its maximum,
   which puts a fictitious point charge at the feed and makes the deposition singular there
   (§6.2). Everything more than about 3 mm away is unaffected to better than 1 %, and the headline
   numbers move by ±15 % depending on how the feed region is excluded — but the *location* of the
   late-time hot spot is set by the artefact, not by physics. A real feed junction is a genuine hot
   spot that this model cannot size. **Instrument it on the bench.**

3. **Constant properties.** $\sigma_\text{eff}$ rises with temperature and $\varepsilon'$ falls, so
   the exact power scaling is good to 10–15 % once the gel is well above 37 °C, and it
   under-predicts (§10K). Removing this needs a coupled electromagnetic–thermal iteration.

4. **Above 100 °C nothing in this model is true.** The gel boils, water vapour changes both the
   permittivity and the geometry, and in a real ablation this is the regime of steam pops. The model
   will happily report 1275 °C. `S.flags` tells you when you have left the physical range.

5. **The gel may fail before it boils.** Low-acyl gellan gels are reported to undergo gel–sol
   transitions in the 30–50 °C range at low ionic strength, shifting upward with added cations
   [T8]. The specific 0.9 % saline / 1 % Phytagel formulation has **not** been characterised here.
   It is entirely possible that the phantom stops behaving as a solid well below the temperatures in
   §9. **Check the melting behaviour before trusting any result above ~50 °C.**

6. **The steady state belongs to the region of interest, not to the world.** By 600 s the diffusion
   length is 18.6 mm and the 37 °C walls of the 30 mm block are carrying real heat flux (§10E, J).
   The "indefinitely" row of §9.5 means "in a 30 mm neighbourhood held at 37 °C at its edge", which
   is precisely the homeostasis assumption the brief supplied — but it is an assumption.

7. **No perfusion, because there is none.** The model is of the *phantom*. Real tissue in the same
   field would run cooler, because capillary perfusion adds a distributed sink of comparable weight
   to conduction over these timescales (§2.2). Do not read the phantom temperatures as tissue
   temperatures without adding that term.

8. **Everything inherits the field model's assumptions**, above all the 0.5 mm standoff between the
   conductor and the gel surface, which the field model identifies as its own dominant sensitivity.
   $\Theta\propto q\propto|\mathbf E|^2$, so a ±20 % error in $|\mathbf E|$ is a ±44 % error in
   temperature rise. **The standoff matters twice as much here as it does there.**

9. **Thermal dose thresholds are biological and the phantom is not.** CEM43 is reported as a
   translation of the phantom result into the tissue equivalent it stands for, not as a property of
   Phytagel.

10. **Ideal geometry.** Flat surface, uniform flow, perfect thermal contact between catheter and
    gel, no bubbles, no container walls, no irrigation ports.

---

## 13. Running the code

```
Heatmodel/matlab/
  water_properties.m      water/blood properties, Incropera Table A.6        [T1]
  convection_h.m          h from the 5 L/min flow, four correlations         [T5]
  thermal_config.m        every constant + boundary condition, each sourced
  sar_volume.m            3-D power deposition q(x,y,z) from the field model
  feed_regularise.m       the feed singularity: diagnosis and optional cap
  heat_grid.m             grid, boundary masks, control-volume weights
  heat_lap.m              the discrete Laplacian with its ghost nodes
  heat_solve.m            transient solver + energy audit
  heat_steady.m           steady state, matrix-free PCG
  cem43.m                 thermal dose                                       [T6]
  verify_heat_model.m     the validation suite (9 tests)
  run_heat_model.m        baseline run + 8 figures + results/heat_model.mat
  run_heat_sensitivity.m  22-case parameter study
  make_heat_figures.m     figures 1-8
  make_sensitivity_figure.m  figure 9, redrawable from the saved study
  extra_diagnostics.m     lesion sizes, isotherm depths, pulsing, Dirichlet vs Robin
  heat_scenario.m         *** ONE BENCH RUN: power, period, time -> figures + CSV ***
  make_scenario_figures.m the two figures a scenario produces
  gui_model.m             one fast forward solve: parameters in, SAR + T maps out
  heat_gui.m              *** THE PARAMETER SWEEP: the window ***
  verify_gui_model.m      validation of the fast model against the stored baseline
  run_all_heat.m          all of the above, logged
  T_gel.m                 *** THE DELIVERABLE FUNCTION ***
  P_for_T.m               inverse: power for a temperature limit
  time_to_T.m             inverse: time to a temperature limit
```

```matlab
addpath('matlab'); addpath('Heatmodel/matlab');
verify_heat_model        % ~10 min  - validation (run this first)
run_heat_model           % ~15 min  - baseline + figures + heat_model.mat
run_heat_sensitivity     % ~10 min  - 22-case parameter study
run_all_heat             % all of it, ~40 min, logged to results/full_heat_log.txt

heat_gui                 % ~25 s to open, then a few seconds per solve
verify_gui_model         % ~4 min   - what the fast model costs (12 checks)
```

`heat_gui` and `gui_model` need none of the stored `.mat` files — they re-solve from
`model_config` and `thermal_config` every time. `verify_gui_model` does need
`results/heat_model.mat` and `results/heat_model_robin.mat`, since comparing against them is the
point; it says so and skips those two tests if they are absent.
`results/gui_verification_log.txt` is the console output of a full run.

`results/full_heat_log.txt` is the console output of a complete `run_all_heat`, and every number
quoted in this document can be found in it. Note the annotation partway down: three reporting bugs
in tests 4, 6 and 7 were found and fixed after that run, and the corrected tests are appended at the
end of the same file. The quoted values are the corrected ones.

`run_heat_model` requires `results/field_model.mat` from the field model only for the cross-check in
`verify_heat_model` test 5; the source itself is recomputed from `model_config`, so the two models
cannot silently drift apart.

Tested on MATLAB R2025b. No toolboxes beyond base MATLAB.

---

## References

Numbering continues from the field model's [R1]–[R17]; thermal sources are [T*].

**[T1]** F. P. Incropera and D. P. DeWitt, *Fundamentals of Heat and Mass Transfer*. Wiley.
Table A.6 (thermophysical properties of saturated water) and Ch. 7 (external flow correlations:
laminar and turbulent flat plate, Churchill–Bernstein cylinder in cross flow), Ch. 8
(Dittus–Boelter). Values cross-checked against the open
[Purdue ME300 property tables](https://engineering.purdue.edu/ME300/Spring%202020/ME%20300%20Property%20Table%209ed.pdf).

**[T2]** M. Choi *et al.*, "Acoustic and thermal characterization of agar based phantoms used for
evaluating focused ultrasound exposures," *J. Therapeutic Ultrasound*, vol. 5, art. 14, 2017.
doi:[10.1186/s40349-017-0093-z](https://doi.org/10.1186/s40349-017-0093-z) ·
[PMC5452295](https://pmc.ncbi.nlm.nih.gov/articles/PMC5452295/) — measured
$k = 0.52\pm0.06$ and $0.57\pm0.10$ W m⁻¹K⁻¹, $c = 4130$ and $4090$ J kg⁻¹K⁻¹,
$\rho = 1050$ and $1070$ kg m⁻³ for 2 % agar brain- and muscle-equivalent phantoms.

**[T3]** H. Zhang, L. Che, X. Chen, F. Zhao, Y. Yang, X. Zhong and Z. Lu, "Experimental
determination of thermal conductivity of water–agar gel at different concentrations and
temperatures," *J. Chem. Eng. Data*, vol. 56, no. 4, pp. 859–864, 2011.
doi:[10.1021/je100570h](https://doi.org/10.1021/je100570h) — line-heat-source measurements,
mass fractions 0.001–0.1, 278–323 K; $k$ falls with polymer concentration and rises with
temperature.

**[T4]** IT'IS Foundation, *Tissue Properties Database — Thermal Properties*.
[itis.swiss](https://itis.swiss/virtual-population/tissue-properties/database) — muscle
$k=0.49$ W m⁻¹K⁻¹, $c=3421$ J kg⁻¹K⁻¹, $\rho=1090$ kg m⁻³; blood $\rho=1050$, $c=3617$,
$k=0.52$ W m⁻¹K⁻¹.

**[T5]** The four convection correlations are all standard textbook results and are quoted in
`convection_h.m` in the form used: laminar flat plate
$\overline{\mathrm{Nu}}_L=0.664\,\mathrm{Re}_L^{1/2}\Pr^{1/3}$; mixed turbulent flat plate
$\overline{\mathrm{Nu}}_L=(0.037\,\mathrm{Re}_L^{4/5}-871)\Pr^{1/3}$; Dittus–Boelter
$\mathrm{Nu}_D=0.023\,\mathrm{Re}_D^{4/5}\Pr^{0.4}$; Churchill–Bernstein for a cylinder in cross
flow. Source [T1], Ch. 7–8.

**[T6]** S. A. Sapareto and W. C. Dewey, "Thermal dose determination in cancer therapy,"
*Int. J. Radiat. Oncol. Biol. Phys.*, vol. 10, no. 6, pp. 787–800, 1984.
[PubMed 6547421](https://pubmed.ncbi.nlm.nih.gov/6547421/) — the CEM43 definition,
$R=0.5$ above 43 °C and $R=0.25$ below.

**[T7]** P. Keangin, P. Rattanadecho *et al.* and the review by C. Rossmann and D. Haemmerich,
"Review of temperature dependence of thermal properties, dielectric properties, and perfusion of
biological tissues at hyperthermic and ablation temperatures," *Crit. Rev. Biomed. Eng.*, 2014.
[PMC4859435](https://pmc.ncbi.nlm.nih.gov/articles/PMC4859435/) — tissue $\sigma$ and
$\varepsilon'$ at 2.45 GHz drift at −0.13 to −0.20 %/°C and −0.15 to −0.18 %/°C respectively;
thermal conductivity rises at +0.004 to +0.78 %/°C.

**[T8]** M. Milas, K. Shi and M. Rinaudo, and subsequent DSC/rheology studies of gellan sol–gel
transitions, e.g. *Carbohydrate Polymers* and *Food Hydrocolloids* series — low-acyl gellan gels
exhibit gel–sol transitions in the 30–50 °C range at low ionic strength, shifting upward with added
cations. Cited only as a caution (§12, limitation 5): the thermal stability of the specific
0.9 % saline / 1 % Phytagel formulation has **not** been established here.
