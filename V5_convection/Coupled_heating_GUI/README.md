# V5 coupled heating GUI

This GUI models dielectric heating plus catheter surface heating of the thermochromic Phytagel under the insulated 2.45 GHz monopole.
It is built from:

- `../Experimental_setup_table.xlsx`: the constants and starting conditions.
- `../Catheter_heating/Convection model.docx`, sections 1–6: the heat equation, the catheter boundary condition, the superposition split and the first-guess constants.
- `../Insulated_monopole_Erms/`: the E-field.
- `../../50W CW 6 min/`: the photographs, used for the colouring and for comparison.

```matlab
>> cd V5_convection/Coupled_heating_GUI
>> v5_gui            % setup table: 60 W, 27 mm antenna, 2 mm gap, 37 °C water
>> v5_gui('photo')   % photographed run: 50 W, 12 mm, 1.1 mm, 25 °C
>> v5_selftest       % checks against closed forms (all pass)
```

## What is on screen

**Model tab.** The layout has three columns.

| Where | What |
|---|---|
| Left | The controls you will use most. Values can be typed in the box beside each slider. **↻** means the stored solution is recomputed; everything else updates instantly. |
| Centre | The antenna plane (y = 0), coloured like the photographs. Water, cable, feed block and insulated element are drawn to scale. The dashed line is the isotherm at the colour-start temperature, and ○ marks the hottest point. Click the map to place a probe. Below the map are the time slider, Play and speed, and the history plot: catheter tip, peak in the plane, surface maximum, probe, and the pigment band shaded. |
| Right | The **power envelope** at the antenna: pulses, the per-frame average the heating uses, and the mean over the treatment. Tick **show photographs** to swap in the photograph nearest the current time, registered to the same mm axes (feed at x = 0, gel surface at z = 0). Under it, the depth profile through the hottest column (or the probe), split into T_b + θ_d (dielectric) and T_b + θ_s (catheter diffusion). At the bottom, a readout of the derived quantities. |

The left-column controls are:

| Group | Controls | Why they are on the front |
|---|---|---|
| Starting conditions | input power, antenna length ↻, antenna–gel distance ↻ | the table's own "Starting conditions" rows |
| Experiment | power drop-off, water temperature, initial gel temperature, circulation flow | what changes from run to run |
| Power envelope | treatment time, duty cycle, pulse period, burst length, rest between bursts | the power timing; all instant |
| Pigment | colour start, band width | decides what the picture shows |
| View | pigment (as photographed / vivid), T, rise, θ_d, θ_s, SAR | |
| Resolution and runtime | 3-D block or 2-D slice, grid spacing Δh 1–0.15 mm, Solve | shows the last solve time (E-field, modes, frames), the node count, and the evaluate and draw times per update |

**Advanced constants tab.** Every other constant is here, each tagged *instant* or *re-solve*:

- reflected power and γ
- antenna radii, ε_ins, feed block
- gel and water permittivities, gel ρ, c, k
- h_s (Inf pins the surface) and h_w
- room temperature
- β and L_c
- catheter m, c, k_c
- gel block size and position
- stored run length and frame interval

A derived panel shows σ, the diffusivity, √(at), k_L, λ_L/4, Z_in, α, P_c, ΔT_c and τ_c. It also shows the node count and the **estimated solve time** before you press Apply. The estimate is refitted to your machine after every solve.

## Model (as in the docx)

$$\rho c\,\partial_t T = k\nabla^2 T + Q_d,\qquad Q_d = \sigma E_\text{rms}^2\,\frac{P_\text{del}}{1\text{ W}},\qquad P_\text{del} = (P_\text{in}-P_\text{refl})(1-P_\text{dropoff})\ \text{while on}$$

- **Surface:** $-k\,\partial_z T = h_s\,(T - T_\text{diff})$ with $T_\text{diff} = T_b + \alpha(x)\,(T_\text{cath}-T_b)$.
  Here $\alpha = \min(1, \beta/(d\sqrt f))$ on $0\le x\le L_c$ and 0 elsewhere.
- **Other faces:** Robin with $h_w$ to room air. The plane $y=0$ is a symmetry plane.
- **Power envelope:** $P_\text{del}\,p(t)$, where $p=1$ in the on-phase of each pulse (period $T_\text{pulse}$, on for duty × period), inside bursts of $t_\text{burst}$ separated by rests of $t_\text{rest}$, until the treatment time $t_\text{on}$. The pulse train restarts with each burst.
  The model is linear, so $\theta_d(t_n)=P_\text{del}\sum_m \bar p_m\,[D(t_n-t_m)-D(t_n-t_{m+1})]$, with $\bar p_m$ the mean of $p$ over frame $m$. This is exact for pulses longer than a frame. Shorter pulses enter through their average (the diffusion length in one 2 s frame is about 0.5 mm). Lower the frame interval in the Advanced tab to resolve faster pulses.
- **Catheter:** lumped, $mc\,\dot T = P_c\,p(t) - k_c(T-T_b)$ with $P_c=\gamma(P_\text{in}-P_\text{refl})P_\text{dropoff}$, stepped exactly per frame. It relaxes after switch-off and between pulses.
- **Superposition:** $T = T_b + P_\text{del}[D(t)-D(t-t_\text{on})] + \sum_k \Delta\Theta_k\,U(t-t_k) + (T_\text{room}-T_b)W + (T_0-T_b)Z$.
  The first three terms are the docx §3 split. W and Z are the third part the docx mentions for $T_0 \ne T_\text{room} \ne T_b$.
- **Stored responses:** D, U, W and Z are what `v5_solve` stores. Everything else is instant (`v5_eval`, 10–60 ms).
- **E_rms:** exactly the `insulated_monopole_Erms.m` physics: King's $k_L$, a standing wave on the element, return charge on the feed block, and the thin-wire field with the gel Green's function.
  The field depends only on (x, distance from the axis), so it is computed on that 2-D table and interpolated onto the 3-D grid. `v5_selftest` reproduces the original script's printed peaks.

**Solver.** Finite volumes on a box grid. The operator separates into three 1-D problems, so their generalised eigenvectors diagonalise it. Each mode is then integrated **exactly in time**: there is no time step, and only Δh sets the accuracy.

Typical runtimes on this machine:

| Grid | Nodes | Solve time |
|---|---|---|
| 3-D, Δh = 0.5 mm | 213 k | 1.5 s |
| 3-D, Δh = 0.35 mm | 606 k | 3.1 s |
| 3-D, Δh = 0.25 mm | 1.66 M | 7.2 s |
| 2-D | — | < 0.3 s |

Times are after MATLAB warm-up (the first solve of a session is slower, up to about 5 s at 0.5 mm) and depend on what else the machine is doing. The block was enlarged on 2026-09-28, so node counts are about 1.9× the first version's.
The answer is converged by 0.5 mm: on the first, 30 mm block the 332 s peak for the table preset was 74.8 / 75.0 / 75.0 °C at Δh = 1 / 0.5 / 0.35 mm. The enlarged block gives 74.81 / 74.81 / 74.84 °C at 0.5 / 0.35 / 0.25 mm.
The 2-D slice runs about 12 K hot, because heat cannot spread sideways (the docx §5 caveat).

## Colouring (the pigment)

`v5_tlc.m` maps T to the photographed colours: hazy blue-grey gel, then pink, yellow, green and teal, then blue at the top of the band, then clear above it. The colours are sampled from `photos/photo_332s.jpg`. "Pigment (vivid)" uses the same hues saturated, which is easier to read.

There is no pigment datasheet, so the band was **fitted**. The photographed bloom edges (behind, in front and below) sit at 40–42 °C in this model at 29, 91, 181 and 332 s. The core is pink at 29 s, green at 91 s, and teal with blue hints at 332 s, which gives colour start **40.5 °C** and width **23 K**.
The start trades off against the power drop-off (both scale the rise), so treat the pair as one calibration. Both presets use the same band, because it is a property of the gel.
Consequence for the table preset: at 60 W with 37 °C water the peak passes the top of the band (63.5 °C) at about 92 s and reaches 75.5 °C by switch-off. The middle should then look **blue, then clear**, not teal.

## Photo comparison (photo preset, 3-D, Δh = 0.5 mm)

| t | Photo: bloom extent | Model: colour-start isotherm | Photo core | Model peak |
|---|---|---|---|---|
| 91 s | x −2.5…9, z −0.8…−7 mm | x −2.5…10, z −0.5…−7.0 mm | green, (1–5, −1.5…−3) | 53.5 °C at (3.5, −2.5) |
| 332 s | x −5.5…12, z −0.3…−12.2 mm | x −6.5 (the block edge)…12, z 0…−12.0 mm | teal, (2–6, −2…−4) | 61.6 °C at (3.5, −3.0) |

**Residual:** the model's top edge is 3–8 K warmer than the other edges. By 332 s parts of the surface reach the colour start, where the photographs still show a 0.3 mm clear layer. A larger h_s, or a smaller β or γ, moves it.

## Assumptions to confirm (each is one field)

- **27 mm antenna radii.** The table's 1.8 mm is read as the insulated element diameter (b = 0.9 mm). The conductor a = b/3.35 is the 50 Ω PTFE ratio, **assumed**.
- **Solved block.** It is the table's 30 × 30 mm gel enlarged 20 % in length and width (36 mm) and 30 % in depth (39 mm), as requested.
  The 27 mm element is **assumed** centred on it (x_start = −4.5 mm).
  The photo preset grows the photographed block (−6.5…23.5 mm) by 3 mm on each side (−9.5…26.5 mm). The gel behind x = −6.5 mm was water or tank in the photographs, so it is modelled but not physical. It hardly changes the result (the 332 s peak is 61.3 vs 61.6 °C).
- **β, m, c_cat, k_c, γ, h_s.** These are the docx §6 first guesses. A thermocouple on the tip fixes only ΔT_c = P_c/k_c and τ_c.
- **Irrigation flow.** The table's 0.03 L/min is not in the model; only the 5 L/min circulation enters α.
- **Antenna–gel distance.** d is a re-solve control here, not instant as the docx lists it, because it moves the E-field as well as α.
- **Constant properties.** Properties are held constant with temperature, since the split needs linearity; they actually drift by 1–2 %/K. There is no perfusion.

## Files

| File | Role |
|---|---|
| `v5_gui.m` | the GUI |
| `v5_defaults.m` | every constant, both presets, with sources |
| `v5_erms.m` | insulated-monopole E_rms on the (x, ρ) table, at 1 W |
| `v5_solve.m` | the four stored responses (exact modal solution), with timings |
| `v5_eval.m` | instant settings → T(x, z, t), catheter, α |
| `v5_envelope.m` | the power envelope: duty cycle, pulse period, bursts and rests |
| `v5_key.m` | which constants force a re-solve |
| `v5_tlc.m` | thermochromic colouring |
| `v5_selftest.m` | erfc, early-time adiabatic rise, adiabatic Z, linearity, switch-off, 2-D = 3-D, the E-field vs the original script, and the envelope (pulses, rests, mean power, catheter vs analytic) |
| `photos/` | mirrored, registered 25 px/mm crops of IMG_7907–7918 with EXIF times (`photos.csv`, `README.txt`) |
| `exports/` | written by "Export frame": a PNG plus a .mat with the settings and T plane |
