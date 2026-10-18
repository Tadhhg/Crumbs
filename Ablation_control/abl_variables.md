# `abl.m`: equations and variables

## Governing equations

| | |
|---|---|
| Power split | $P_{del} = (1-\delta)\,P_{gen}, \qquad P_{lost} = \delta\,P_{gen}$ |
| Dielectric source | $Q(\mathbf r) = \sigma\,E_{rms}^2(\mathbf r)\,P_{del}, \qquad \sigma = \omega\,\varepsilon_0\,\varepsilon_g''$ |
| Gel | $\rho c\,\dfrac{\partial T}{\partial t} = k\,\nabla^2 T + Q$ |
| Surface ($z=0$) | $-k\,\dfrac{\partial T}{\partial z} = h_s\,(T - T_w), \qquad T_w = T_b + \alpha\,(T_{cat} - T_b)\ \text{on}\ 0 \le x \le L_c$ |
| Walls | $-k\,\dfrac{\partial T}{\partial n} = h_w\,(T - T_b)$ (surroundings held at body temperature) |
| Catheter tip | $m_{cat}\,c_{cat}\,\dfrac{dT_{cat}}{dt} = \gamma\,P_{lost} - k_{cat}\,(T_{cat} - T_b)$ |
| Surface coefficient | $h_s = h_{ref}\,(f_w / f_{ref})^{n}$ (laminar flat plate: $Nu \propto Re^{1/2}$, so $n = 1/2$) |
| Surface transfer | $\alpha = \min\!\left(1,\ \dfrac{\beta}{d\,\sqrt{f_w}}\right)$ |
| Ablation | $T_{max}(\mathbf r) = \max_t T(\mathbf r,t), \qquad \text{ablated where } T_{max} \ge T_{abl}$ |

Coordinates: $x$ along the element (feed at $x=0$), $y$ across it, $z$ up (gel surface at $z=0$, gel at $z<0$). The model is symmetric about $y=0$, the antenna plane.

## Parameters (`P = abl.params`)

### Antenna
| Code | Symbol | Value | Meaning | Units |
|---|---|---|---|---|
| `P.f` | $f$ | $2.45\times10^9$ | microwave frequency | Hz |
| `P.L` | $L$ | $27\times10^{-3}$ | element length, feed to tip | m |
| `P.a` | $a$ | $0.27\times10^{-3}$ | inner conductor radius | m |
| `P.b` | $b$ | $0.9\times10^{-3}$ | outer radius of the insulation | m |
| `P.eps_ins` | $\varepsilon_{ins}$ | $2.1$ | insulation (PTFE) relative permittivity | – |
| `P.l_f` | $l_f$ | $1.5\times10^{-3}$ | exposed feed block behind the feed | m |
| `P.d` | $d$ | $2\times10^{-3}$ | element axis to gel surface | m |

### Media
| Code | Symbol | Value | Meaning | Units |
|---|---|---|---|---|
| `P.eps_g` | $\varepsilon_g = \varepsilon_g' - j\varepsilon_g''$ | $43.98 - 13.5j$ | gel complex permittivity | – |
| `P.eps_w` | $\varepsilon_w = \varepsilon_w' - j\varepsilon_w''$ | $73 - 8.8j$ | water complex permittivity | – |
| `P.rho` | $\rho$ | $1000$ | gel density | kg/m³ |
| `P.c` | $c$ | $4180$ | gel specific heat | J/(kg K) |
| `P.k` | $k$ | $0.6$ | gel thermal conductivity | W/(m K) |

### Temperatures and boundaries
| Code | Symbol | Value | Meaning | Units |
|---|---|---|---|---|
| `P.T_b` | $T_b$ | $37$ | circulating water (blood) temperature | °C |
| `P.T_0` | $T_0$ | $37$ | gel temperature at $t=0$ | °C |
| `P.h_s` | $h_{ref}$ | $1400$ | water-to-gel surface heat transfer coefficient at $f_{ref}$ (laminar flat plate, V4) | W/(m² K) |
| `P.flow_ref` | $f_{ref}$ | $5$ | flow $h_{ref}$ was derived at | L/min |
| `P.n_h` | $n$ | $0.5$ | flow exponent, $h_s \propto f_w^n$: 0.5 laminar plate, 0.8 turbulent, 0 developed pipe flow | – |
| `P.h_w` | $h_w$ | $8$ | walls to surroundings at $T_b$ heat transfer coefficient | W/(m² K) |

### Surface transfer from the catheter
| Code | Symbol | Value | Meaning | Units |
|---|---|---|---|---|
| `P.L_c` | $L_c$ | $27\times10^{-3}$ | heated catheter footprint on the surface | m |
| `P.beta` | $\beta$ | $9\times10^{-4}$ | transfer constant in $\alpha$ | m (L/min)$^{1/2}$ |
| `P.flow` | $f_w$ | $5$ | circulation flow rate | L/min |

### Power split and catheter
| Code | Symbol | Value | Meaning | Units |
|---|---|---|---|---|
| `P.drop` | $\delta$ | $0.7$ | fraction of $P_{gen}$ lost before the antenna | – |
| `P.gamma` | $\gamma$ | $0.043$ | share of $P_{lost}$ that heats the catheter tip | – |
| `P.m_cat` | $m_{cat}$ | $2\times10^{-3}$ | heated tip mass | kg |
| `P.c_cat` | $c_{cat}$ | $1500$ | tip specific heat | J/(kg K) |
| `P.k_cat` | $k_{cat}$ | $0.1$ | tip-to-water heat loss conductance | W/K |

### Ablation
| Code | Symbol | Value | Meaning | Units |
|---|---|---|---|---|
| `P.T_abl` | $T_{abl}$ | $50$ | ablation cutoff temperature | °C |

### Grid
| Code | Symbol | Value | Meaning | Units |
|---|---|---|---|---|
| `P.gel` | $[L_x,\ W_y/2,\ H_z]$ | $[36,\ 18,\ 39]\times10^{-3}$ | gel block length, half-width, depth | m |
| `P.x0` | $x_0$ | $-4.5\times10^{-3}$ | block edge behind the feed | m |
| `P.dh` | $\Delta h$ | $0.5\times10^{-3}$ | grid spacing | m |

### Graded grid (`abl_forvibes.m` only)

`abl_forvibes.m` is a copy of `abl.m` that also accepts `P.grid` (point positions `.x`, `.y`, `.z` in m), which replaces `P.dh`. `abl_forvibes.graded(lo, hi, a, b, h_fine, h_coarse, growth)` spaces points `h_fine` apart in $[a, b]$ and grows each gap by `growth` outward, up to `h_coarse`. `abl_forvibes.modes` then uses each point's own length $w_i = (\Delta_{i-1} + \Delta_i)/2$ and gap conductances $1/\Delta_i$. `S.cell` holds the plane area of each point, for the ablated area. `abl_pid_viberuntimereduction.m` runs the controller on such a grid (0.25 mm in the top 6 mm and within 5 mm of the antenna plane, 0.5 mm along x, up to 2 mm far away).

## Solver (`abl.init`, `abl.step`)

The gel temperature rise is expanded in the eigenmodes $\phi_m$ of $k\nabla^2$ with the convection conditions on the faces:

$$\theta(\mathbf r,t) = T - T_b = \sum_m c_m(t)\,\phi_m(\mathbf r), \qquad \dot c_m = -\lambda_m c_m + g_m$$

$$g_m = P_{del}\,g_{Q,m} + \Theta\,g_{U,m}, \qquad \Theta = \alpha\,(T_{cat} - T_b)$$

The walls exchange heat with surroundings at $T_b$, so they add no forcing term; they only enter through the Robin condition in the modes.

| Code | Symbol | Meaning | Units |
|---|---|---|---|
| `S.gQ` | $g_Q$ | dielectric source per W delivered, projected on the modes | K m$^{3/2}$/(s W) |
| `S.gU` | $g_U$ | surface water 1 K warmer over $0 \le x \le L_c$ | m$^{3/2}$/s |
| `S.c` | $c_m$ | modal coefficients of $\theta$; $c_m(0)$ from $\theta(0) = T_0 - T_b$ | K m$^{3/2}$ |
| `S.Lam` | $\lambda_m = \kappa\,(\lambda_x + \lambda_y + \lambda_z)$ | modal decay rates, $\kappa = k/\rho c$ | 1/s |
| `S.h_s` | $h_s$ | surface heat transfer coefficient at the current flow | W/(m² K) |
| `S.alpha` | $\alpha$ | catheter-to-surface-water transfer | – |
| `S.tau` | $\tau_{cat} = m_{cat}c_{cat}/k_{cat}$ | catheter tip time constant | s |
| `Th` | $\Theta$ | surface water rise, mean over the step | K |

Each step of length $\Delta t$ at constant $P_{gen}$ is exact:

$$T_{cat}^+ = T_b + (T_{cat} - T_b)\,e^{-\Delta t/\tau_{cat}} + \frac{\gamma\,\delta\,P_{gen}}{k_{cat}}\left(1 - e^{-\Delta t/\tau_{cat}}\right)$$

$$c_m^+ = e^{-\lambda_m \Delta t}\,c_m + \frac{1 - e^{-\lambda_m \Delta t}}{\lambda_m}\,g_m$$

$$T_{max}^+ = \max(T_{max},\ T^+)$$

## Controller states (`S`)

| Code | Symbol | Meaning | Units |
|---|---|---|---|
| `S.T` | $T(x,z)$ | temperature on the antenna plane $y=0$ ($n_x \times n_z$) | °C |
| `S.Tc` | $T_{cat}$ | catheter tip temperature | °C |
| `S.Tmax` | $T_{max}(x,z)$ | peak temperature so far | °C |
| `S.t` | $t$ | elapsed time | s |
| `abl.depth(S)` | $d_{abl}$ | deepest point with $T_{max} \ge T_{abl}$, interpolated between nodes | m |

## Changing the flow mid-run (`abl.reflow`)

`S = abl.reflow(S, flow)` switches to a new flow without losing the heat in the gel. Flow changes $h_s = h_{ref}(f_w/f_{ref})^n$ and $\alpha$; only the top surface feels $h_s$, so only the z-shapes change. The mode amounts are carried across exactly:

$$M = V_{z,new}^\top\,\mathrm{diag}(w_z)\,V_{z,old}, \qquad c \leftarrow M c,\quad g_Q \leftarrow M g_Q \ \text{(along z)}$$

| Code | Meaning | Units |
|---|---|---|
| `S.wz` | length of gel each z point stands for | m |
| `S.Lxy` | $\kappa(\lambda_x + \lambda_y)$, x and y part of the fading rates | 1/s |
| `S.gUxy` | x and y part of $g_U$ | – |

## Closed-loop control (`abl_pid.m`)

Every `DT` = 0.5 s, on a perfect temperature map:

| | |
|---|---|
| PID 1 (priority) | $e_1 = (T_{abl} - \text{MARGIN}) - \hat y_1$, $\hat y_1$ = hottest the spared layer ($z$ shallower than `ABLATION_DEPTH`) gets if the last power runs one more step, then stops for `T_LOOK` |
| PID 2 | $e_2 = R_{target} - \dot A$, $A$ = area with $T_{max} \ge T_{abl}$ on the antenna plane, $\dot A$ smoothed over `TAU_R` |
| Applied power | $P_{gen} = \min(P_1, P_2)$, $P_2 \in [P_{MIN}, P_{MAX}]$ |
| Safety | if "$P_{gen}$ for one step, then off for `T_LOOK`" brings the layer within MARGIN/2 of $T_{abl}$, halve $P_{gen}$ until safe |
| Anti-windup | the loop not in charge has its integral set to the applied power (external reset) |

| Setting | Default | Meaning | Units |
|---|---|---|---|
| `ABLATION_DEPTH` | 1.5e-3 | spared surface layer (no ablation shallower); input 1 of `abl_pid`, mm or m | m |
| `MAX_WATTS` | 60 | highest power PID 2 may ask for | W |
| `P_MIN` | 10 | lowest power PID 2 may ask for (PID 1 can go lower) | W |
| `MODE` | `'pid'` | `'hold'` = fixed `P_HOLD`, PID 1 still guards | – |
| `MARGIN` | 2 | PID 1 set point below $T_{abl}$ | K |
| `T_LOOK` | 10 | look-ahead with the power off | s |
| `R_TARGET` | 0.5 | target growth rate of the ablated area | mm²/s |
| `Kp1, Ki1` | 8, 0.4 | PID 1 gains | W/K, W/(K s) |
| `Kp2, Ki2` | 10, 1 | PID 2 gains | W/(mm²/s), W/mm² |
