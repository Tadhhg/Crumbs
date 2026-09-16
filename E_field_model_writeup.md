# E-field of a 12 mm 2.45 GHz monopole lying on a 30 mm Phytagel slab

**Working document — model derivation, sources, numerical implementation and results**
Date: 2026-09-06 · Code: `matlab/` · Figures: `figures/` · Numerical output: `results/`

---

## 0. The answer, up front

The excitation is time-harmonic at $f=2.45$ GHz, so the field at every point is a phasor
$\tilde{\mathbf E}(x,z)=(\tilde E_x,\tilde E_y,\tilde E_z)$ and the physical field is
$\mathbf E(x,z,t)=\mathrm{Re}\{\tilde{\mathbf E}e^{j\omega t}\}$. Each Cartesian component is a
pure cosine, so the **magnitude** has an exact closed form:

$$\boxed{\;|\mathbf E(x,z,t)| \;=\; \sqrt{\tfrac12\Big[\,\big|\tilde{\mathbf E}(x,z)\big|^{2} \;+\; \big|Q(x,z)\big|\cos\!\big(2\omega t + \arg Q(x,z)\big)\Big]}\;}$$

$$\big|\tilde{\mathbf E}\big|^{2}=\tilde{\mathbf E}\!\cdot\!\tilde{\mathbf E}^{*}=|\tilde E_x|^2+|\tilde E_y|^2+|\tilde E_z|^2 ,
\qquad Q=\tilde{\mathbf E}\!\cdot\!\tilde{\mathbf E}=\tilde E_x^2+\tilde E_y^2+\tilde E_z^2$$

(note: $Q$ uses the **un-conjugated** dot product and is complex; derivation in §10.1).

Three things follow immediately, and they are the physically important ones:

| quantity | expression | meaning |
|---|---|---|
| time RMS | $E_\text{rms}=\big|\tilde{\mathbf E}\big|/\sqrt2$ | the only thing that sets heating/SAR |
| peak over time | $E_\text{peak}=\sqrt{(\big|\tilde{\mathbf E}\big|^{2}+|Q|)/2}$ | major semi-axis of the polarisation ellipse |
| minimum over time | $E_\text{min}=\sqrt{(\big|\tilde{\mathbf E}\big|^{2}-|Q|)/2}$ | minor semi-axis |

**The magnitude ripples at $2f = 4.9$ GHz**, not at $f$ — because $|\mathbf E|$ is insensitive to the
sign of the field. It only reaches zero where the polarisation is linear ($|Q|=\big|\tilde{\mathbf E}\big|^2$); elsewhere
the field vector traces an ellipse and $|\mathbf E|$ never vanishes.

$\tilde{\mathbf E}(x,z)$ itself has **no elementary closed form** — that is the whole difficulty of the
problem (§5). It is computed by an exact spectral (Sommerfeld) solution of the three-layer
problem, validated against an independent analytic solution to $10^{-4}$–$10^{-3}$ (§7), and delivered two ways:

* **Exact:** `E_gel(x,z,t)` in `matlab/E_gel.m`, interpolating the computed phasor field.
* **Compact analytic approximation** (§10.3): a single line integral, accurate to **1.4 dB rms /
  2 dB for 90 % of the volume** over $4\le R\le 28$ mm.

Absolute scale is set by the delivered power, not by an assumed antenna field strength — see §6.4.
At $P_\text{del}=150$ W (200 W forward − 50 W max reflected) the model gives, on the plane $y=0$:

| position (x, depth) | $E_\text{rms}$ (V/m) | $E_\text{peak}$ (V/m) | SAR (W/kg) |
|---|---|---|---|
| 6 mm, 0.5 mm | 13 970 | 17 340 | 531 000 |
| 6 mm, 2 mm | 7 640 | 9 460 | 158 900 |
| 6 mm, 5 mm | 4 400 | 5 970 | 52 800 |
| 6 mm, 10 mm | 2 340 | 3 280 | 14 870 |
| 6 mm, 20 mm | 690 | 975 | 1 300 |
| 6 mm, 29 mm | 496 | 702 | 671 |

Probe depths are snapped to the nearest 0.25 mm cell centre; `E_gel` interpolates exactly
(e.g. 4 320 V/m at precisely 5.000 mm depth). Everything scales as $\sqrt{P_\text{del}}$ —
multiply by 1.155 for 200 W delivered.

---

## 1. What was specified, and what I had to assume

**Given by you:**

| symbol | value | note |
|---|---|---|
| $f$ | 2.45 GHz | ISM band |
| $L$ | 12 mm | monopole length |
| $d$ | 30 mm | gel depth |
| $P_\text{fwd}$ | ≤ 200 W | generator setting |
| $P_\text{ref}$ | ≤ 50 W | generator max reflection |
| medium | Phytagel, lossy | |
| geometry | antenna **parallel** to the gel surface | |

**Assumed, because it was not specified.** Each is carried as an explicit parameter and swept in §9.

| assumption | baseline | why | sensitivity of $|E|$ |
|---|---|---|---|
| A1 Delivered power | $P_\text{del}=P_\text{fwd}-P_\text{ref}=150$ W | worst case; you may be better matched | exact $\sqrt{P}$ scaling |
| A2 Gel composition | 1 % w/v Phytagel in 0.9 % w/v NaCl saline, 25 °C | the standard formulation (§4) | ±15 % (§9C) |
| A3 Standoff $h$ | 0.5 mm (wire axis above surface, i.e. resting on it) | conductor radius + sheath | **largest**: ±20 % on $|E|$ at depth, **12×** on the near-surface $1/e$ distance (§9B) |
| A4 Current wave $k_e$ | $L=\lambda_e/4$ (resonant) | consistent with the ≤ 50 W reflection (§6.3) | ±8 % (§9A) |
| A5 What is under the gel | air (open dish) | | ±35 % but only below 25 mm (§9E) |
| A6 Gel is homogeneous, non-magnetic, at constant temperature | | | see §11 |
| A7 CW drive | | if pulsed, scale by duty cycle | |

**The "2D plane".** I take it to be the vertical plane containing the antenna axis, $y=0$, spanning
$-30\le z\le 0$ mm. This is the plane in which the deposition pattern is asymmetric and therefore the
informative one. The solver produces the full 3D field, so any other plane is available at no extra
cost (`fig4`, `results/field_model.mat` also stores a horizontal slice).

---

## 2. Literature search and screening

You asked me to reject work on **conducting media** and on **different plane geometries**. That
screen is physically well-founded here, and it is worth saying exactly why.

**Is this medium "conducting" or a "lossy dielectric"?** The discriminant is the loss tangent
$\tan\delta = \varepsilon''/\varepsilon'$. For this gel at 2.45 GHz, $\tan\delta = 0.269$ (§4), i.e.
displacement current beats conduction current by roughly 4:1. It is a **lossy dielectric**, not a
conductor. The classical half-space antenna literature (Sommerfeld, Baños, King) is written for
$\tan\delta \gg 1$ — earth at VLF, seawater — where quasi-static approximations to the Sommerfeld
kernel are legitimate. Those approximations are *not* valid here. Your screen is correct.

### 2.1 Papers kept

| # | Reference | Why it applies |
|---|---|---|
| K1 | Michalski & Mosig 1997, *IEEE TAP* **45**, 508–519 [R1] | The transmission-line / spectral-domain formulation for **plane-stratified** media. This is the machinery I use (§5.3). Not restricted to conductors. |
| K2 | Chew 1995, *Waves and Fields in Inhomogeneous Media* [R2] | Weyl plane-wave (Sommerfeld) representation and layered-media reflection recursion. Planar, arbitrary complex $\varepsilon$. |
| K3 | King & Smith 1981, *Antennas in Matter* [R3] | The reference text for antennas embedded in / adjacent to dissipative **and dielectric** media; source of the insulated-antenna argument in §6.3. |
| K4 | Gao, Šipuš & Skrivervik 2022, arXiv:2205.04971 [R4] | Uses exactly the Green's-function-for-multilayered-planar-media approach at **2.45 GHz** in a lossy **planar** phantom; shows the same normalised E-field map style. *Caveat: their source is buried at 3 cm depth, not lying on the surface — so it is a method/benchmark reference, not a geometry match.* |
| K5 | Bernardi, Cavagnaro & Lin 2004, *IEEE T-MTT* [R5] | 2450 MHz **cardiac ablation catheter antenna** in a homogeneous dielectric phantom; SAR + temperature; FDTD vs measurement. Closest published application. |
| K6 | Lin et al. 2004, *J. Interv. Card. Electrophysiol.* [R6] | 2.45 GHz **monopole** catheter, phantom with myocardium-equivalent properties, antenna **in contact with the tissue surface**, and explicitly discusses partial contact. The nearest geometry match I found. |
| K7 | Balanis 2005, *Antenna Theory* 3rd ed. [R7] | Sinusoidal current distribution Eq. (4-56); dipole fields. |
| K8 | Stogryn 1971 [R8] / Klein & Swift 1977 [R9] / Gadani et al. 2012 [R10] | Complex permittivity of NaCl solutions (§4). |
| K9 | Gabriel et al. 1996 [R11,R12] | Tissue permittivity, for the muscle-equivalence comparison only. |
| K10 | Qian et al. 2020 [R13] / Qian et al. [R14] | Source of the Phytagel phantom formulation. **See the caveat in §4.4.** |

### 2.2 Papers rejected, and the reason

| Reference | Reason for rejection |
|---|---|
| Michalski & Mosig 2018, *Radio Science*, HED over an imperfectly **conducting** half-space, with plasmonics extensions | Conducting lower medium **and** far-zone asymptotics only. We need the reactive/radiating near field within 2 wavelengths. |
| Baños/King, *The Horizontal Electric Dipole in a Conducting Half-Space* | Conducting medium. |
| NIST, *Dipole radiation in a conducting half space* | Conducting medium. |
| Parise, *Horizontal linear antennas above a lossy half-space: a new model for the Sommerfeld kernel* | Geometry matches, but the kernel model is developed for and validated on conduction-dominated ground at low frequency. Noted, not used. |
| Multi-slot / dual-slot / helical-slot / cap-choke **interstitial** coaxial ablation antennas (numerous) | **Different geometry** — applicator inserted *into* the medium along its axis, nearly cylindrically symmetric, effectively an unbounded medium. No planar air/gel interface parallel to the antenna. |
| Monopole patch antennas for head imaging; wearable/textile monopoles | Different geometry (patch on a substrate radiating toward a curved body), and the object of interest is a link budget, not an internal E-field map. |
| Waveguide, horn and loop **hyperthermia applicators** | Aperture sources, not a wire parallel to the surface. |
| *Effect of non-parallel applicator insertion on ablation zone* | Inserted applicator. |

### 2.3 The gap

**No published paper solves this exact problem**: a bare/thinly-insulated 12 mm wire lying *on* the
planar surface of a lossy-dielectric slab of finite depth at 2.45 GHz, with the internal E-field map
as the output. K5 and K6 are the closest in application; K1–K2 give the correct machinery. That gap
is why the model below is built from first principles rather than taken from a paper.

---

## 3. Coordinate system and geometry

```
                    z
                    ^            air  (layer 0,  eps_r = 1)
   feed             |    tip
     x=0 ===========|=========== x = L = 12 mm      wire axis at z = +h
  ---------------------------------------------  z = 0   gel surface
                    |
                    |            GEL (layer 1)  eps_r = 74.29 - j19.98
                    |            d = 30 mm
                    |
  ---------------------------------------------  z = -d
                                 air  (layer 2)   [or metal dish, §9E]
```

* $x$ — along the antenna axis. The 12 mm element occupies $0\le x\le L$. $x=0$ is the **feed**,
  $x=L$ the **open tip**.
* $y$ — transverse horizontal. The modelled plane is $y=0$.
* $z$ — height. Gel surface at $z=0$; **$z$ is negative inside the gel**; floor at $z=-30$ mm.

Time convention throughout: $e^{+j\omega t}$, so a lossy medium has
$\varepsilon_r = \varepsilon' - j\varepsilon''$ with $\varepsilon''\ge0$, and a wave travelling in
$+z$ carries $e^{-jkz}$ with $\mathrm{Im}(k)\le 0$.

---

## 4. Material properties — every number and its source

### 4.1 Physical constants

| constant | value | source |
|---|---|---|
| $\varepsilon_0$ | $8.8541878128\times10^{-12}$ F/m | CODATA 2018 recommended values |
| $\mu_0$ | $1.25663706212\times10^{-6}$ H/m | CODATA 2018 |
| $c_0$ | 299 792 458 m/s (exact) | SI definition of the metre |
| $M$(NaCl) | 58.443 g/mol | IUPAC standard atomic weights |
| $\rho_\text{gel}$ | 1000 kg/m³ (assumed = water) | used only to convert power density to SAR |

### 4.2 The gel: model and justification

At 1 % w/v the gellan polymer occupies well under 1 % by volume; the material is >98 % water by
mass. Its complex permittivity is therefore that of the **saline** to well within the model's own
uncertainty. (The one measured study of 1 % gellan gels with added NaCl, Zhang et al. 2015 [R15],
confirms that salt content — not the polymer — is what moves $\varepsilon''$; its numeric tables are
paywalled and were not obtained.)

I therefore compute $\varepsilon_r$ from the **Stogryn (1971)** regression model for NaCl solutions
[R8]. Its coefficients are transcribed from the open-access reproduction by Gadani et al. (2012)
[R10], Eqs. (1)–(10), which is freely readable and was checked line by line. Implementation:
`matlab/saline_stogryn.m`.

$$\varepsilon_r(\omega) \;=\; \varepsilon_\infty \;+\; \frac{\varepsilon_s-\varepsilon_\infty}{1+j\omega\tau} \;-\; j\,\frac{\sigma_i}{\omega\varepsilon_0}
\qquad\text{[R10, Eq. (1)]}$$

$$\varepsilon' = \varepsilon_\infty + \frac{\varepsilon_s-\varepsilon_\infty}{1+(\omega\tau)^2},
\qquad
\varepsilon'' = \frac{(\varepsilon_s-\varepsilon_\infty)\,\omega\tau}{1+(\omega\tau)^2} + \frac{\sigma_i}{\omega\varepsilon_0}$$

with, for temperature $T$ in °C and normality $N$ in mol/L:

| quantity | expression | source |
|---|---|---|
| $\varepsilon_s(T,0)$ | $87.74-0.40008T+9.398\times10^{-4}T^2+1.410\times10^{-6}T^3$ | [R10] Eq. (6) |
| $a(N)$ | $1-0.2551N+5.151\times10^{-2}N^2-6.889\times10^{-3}N^3$ | [R10] Eq. (4) |
| $\varepsilon_s(T,N)$ | $\varepsilon_s(T,0)\,a(N)$ | [R10] Eq. (2) |
| $2\pi\tau(T,0)$ | $1.1109\times10^{-10}-3.824\times10^{-12}T+6.938\times10^{-14}T^2-5.096\times10^{-16}T^3$ | [R10] Eq. (7) |
| $b(N,T)$ | $1+0.1463\times10^{-2}NT-0.04896N-0.02967N^2+5.644\times10^{-3}N^3$ | [R10] Eq. (5) |
| $2\pi\tau(T,N)$ | $2\pi\tau(T,0)\,b(N,T)$ | [R10] Eq. (3) |
| $\sigma_i(25,N)$ | $N\,[10.394-2.3776N+0.68258N^2-0.13538N^3+1.0086\times10^{-2}N^4]$ | [R10] Eq. (10) |
| $\sigma_i(T,N)$ | $\sigma_i(25,N)\,[1-1.962\times10^{-2}\Delta+8.08\times10^{-5}\Delta^2-\Delta N(3.020\times10^{-5}+3.922\times10^{-5}\Delta+N(1.721\times10^{-5}-6.584\times10^{-6}\Delta))]$, $\Delta=25-T$ | [R10] Eq. (9) |
| $\varepsilon_\infty$ | 4.9 | Klein & Swift 1977 [R9] — *"best fit to available high frequency data is achieved if $\varepsilon_\infty$ is assumed a constant value of 4.9"* |

**Normality.** 0.9 % w/v NaCl = 9 g/L ⇒ $N = 9/58.443 = 0.154$ mol/L. Cross-check with Stogryn's own
salinity conversion $N=S(1.707\times10^{-2}+1.205\times10^{-5}S+4.058\times10^{-9}S^2)$ [R10, Eq. (8)]
at $S\approx8.96$ ppt gives $N=0.1539$ — agreement to 0.1 %, confirming the unit conversion.

### 4.3 Resulting baseline values (25 °C, 2.45 GHz)

| quantity | value |
|---|---|
| $\varepsilon_s$ | 75.363 |
| $\tau$ | 8.078 ps |
| $\sigma_\text{ionic}$ | 1.547 S/m |
| $\varepsilon_r$ | $\mathbf{74.29 - j\,19.98}$ |
| $\tan\delta$ | 0.2689 |
| $\sigma_\text{eff}=\omega\varepsilon_0\varepsilon''$ | **2.723 S/m** |
| $\alpha = -\mathrm{Im}\,k_1$ | **58.98 Np/m** → field $1/e$ depth **16.95 mm**, power $1/e$ depth 8.48 mm |
| $\beta = \mathrm{Re}\,k_1$ | 446.49 rad/m → $\lambda_\text{gel}$ = **14.07 mm** |
| $L/\lambda_\text{gel}$ | 0.853 |

Note $\sigma_\text{eff}$ (2.72 S/m) is the **total** effective conductivity — ionic (1.55 S/m) plus
dipolar relaxation (1.18 S/m). All of $\varepsilon''$ heats the gel, so $\sigma_\text{eff}$, not
$\sigma_\text{ionic}$, is what belongs in the SAR expression.

### 4.4 The glycerol variant — an open question in the sourcing

The prior notes in `OLD_dont_use/References.md` identify the published Phytagel formulation as
**70 % 0.9 %-saline + 30 % glycerol + 1 % w/v Phytagel**, "muscle equivalent", from Qian et al.
[R13,R14]. Two things must be said plainly:

1. **That quotation is still unverified.** Both papers are paywalled; the recipe text was taken from
   a search index, not from the PDF. It has not been possible to verify it in this session either.
   Check it against the published PDF through the library before relying on it.
2. **If** glycerol is present, $\varepsilon'$ drops substantially and the medium moves toward muscle.
   Rather than guess at a glycerol–saline mixture permittivity, I run the sensitivity case at the
   *stated design target* — muscle — computed from Gabriel's four-pole Cole-Cole model (§4.5).

If your gel is plain saline Phytagel with no glycerol, §4.3 is the right column and no caveat
applies.

### 4.5 Muscle, for comparison only

Computed, not asserted: `matlab/tissue_cole_cole.m` implements

$$\varepsilon_r(\omega)=\varepsilon_\infty+\sum_{n=1}^{4}\frac{\Delta\varepsilon_n}{1+(j\omega\tau_n)^{1-\alpha_n}}+\frac{\sigma_i}{j\omega\varepsilon_0}
\qquad\text{[R12]}$$

with the muscle parameter column $\varepsilon_\infty=4.0$, $\sigma_i=0.2$ S/m,
$\Delta\varepsilon=(50.0,\;7000,\;1.2\times10^{6},\;2.5\times10^{7})$,
$\tau=(7.23\ \text{ps},\;353.68\ \text{ns},\;318.31\ \mu\text{s},\;2.274\ \text{ms})$,
$\alpha=(0.10,\,0.10,\,0.10,\,0.00)$ from Gabriel et al. 1996 [R11,R12].

**Verification:** at 2.45 GHz this returns $\varepsilon'=52.73$, $\sigma_\text{eff}=1.738$ S/m —
which is exactly the value the IT'IS database serves for muscle (52.7, 1.74 S/m). The same code
returns liver 43.03 / 1.686 S/m and fat 5.28 / 0.105 S/m, both consistent with published values.
That round-trip is the evidence that the parameters were transcribed correctly.

**So the plain-saline gel is roughly 1.6× lossier than muscle** ($\sigma_\text{eff}$ 2.72 vs
1.74 S/m) and has a 40 % higher $\varepsilon'$. It is *not* muscle-equivalent without the glycerol.

---

## 5. Theory

### 5.1 Why there is no closed form

The problem is a horizontal electric dipole (strictly, a line of them) beside a planar interface
between two half-spaces of very different permittivity, with a third interface 30 mm below. This is
the Sommerfeld half-space problem. It has **no solution in elementary functions**; the field is a
Sommerfeld integral. All the classical closed forms in the literature are asymptotic — valid either
in the far zone ($k_1\rho\gg1$) or for a conduction-dominated lower medium. Here:

* the region of interest ($\rho \lesssim 30$ mm) is **2 wavelengths** of $\lambda_\text{gel}=14.07$ mm — near field;
* the medium is **not** conduction-dominated ($\tan\delta=0.27$).

So neither family of approximations applies, and the Sommerfeld integral is evaluated directly.

### 5.2 Plane-wave (Weyl) representation of the source

The $x$-directed current filament is described by a magnetic vector potential
$\mathbf A=\hat{\mathbf x}A_x$ with $\nabla^2\mathbf A+k^2\mathbf A=-\mu\mathbf J$, and

$$\mathbf E=\frac{1}{j\omega\mu\varepsilon}\big[\nabla(\nabla\!\cdot\!\mathbf A)+k^2\mathbf A\big]$$

(standard Lorenz-gauge result; e.g. Balanis [R7] §3.6, [R16] Ch. 6). The free-space Green function is
expanded into plane waves by the **Weyl identity** [R2, §2.2; originally Weyl 1919 / Sommerfeld 1909]:

$$\frac{e^{-jkR}}{4\pi R}=\frac{1}{8\pi^{2}}\iint \frac{e^{-j(k_xx+k_yy+k_z|z|)}}{j k_z}\,dk_x\,dk_y ,
\qquad k_z=\sqrt{k^{2}-k_x^{2}-k_y^{2}},\;\;\mathrm{Im}\,k_z\le0$$

Fourier convention used throughout:
$f(x,y)=\frac{1}{4\pi^2}\iint \tilde f\,e^{-j(k_xx+k_yy)}dk_xdk_y$,
$\tilde f=\iint f\,e^{+j(k_xx+k_yy)}dxdy$.

With $\tilde I(k_x)=\int I(x)e^{jk_xx}dx$ the spectral vector potential in the source layer is
$\tilde A_x = \mu_0\tilde I\,e^{-jk_{z0}|z-h|}/(2jk_{z0})$.

### 5.3 TE/TM decomposition and the three-layer stack

Define $k_\rho=\sqrt{k_x^2+k_y^2}$, the horizontal unit vectors
$\hat{\mathbf s}=(k_y,-k_x,0)/k_\rho$ (TE direction) and $\hat{\boldsymbol\rho}=(k_x,k_y,0)/k_\rho$,
and the **modal wave impedances** [R1]:

$$Z_n^{\text{TE}}=\frac{\omega\mu_0}{k_{zn}},\qquad Z_n^{\text{TM}}=\frac{k_{zn}}{\omega\varepsilon_n}$$

Writing the TE wave by its $E$-amplitude $a$ along $\hat{\mathbf s}$ and the TM wave by its
$H$-amplitude $b$ along $\hat{\mathbf s}$ (this choice removes all up/down sign ambiguity), matching
$E_t$ and $H_t$ at an interface gives

$$r^{01}_{a}=\frac{Z_1^{\text{TE}}-Z_0^{\text{TE}}}{Z_1^{\text{TE}}+Z_0^{\text{TE}}},\qquad
r^{01}_{b}=\frac{Z_0^{\text{TM}}-Z_1^{\text{TM}}}{Z_0^{\text{TM}}+Z_1^{\text{TM}}},\qquad
t^{01}=1+r^{01},\qquad r^{10}=-r^{01}$$

i.e. the TE amplitude behaves as a transmission-line **voltage** and the TM amplitude as a
transmission-line **current** — the Michalski–Mosig analogy [R1]. Projecting the primary field onto
the two polarisations gives the incident amplitudes at $z=0$:

$$a_0^{d}=-\frac{\omega\mu_0}{2k_{z0}}\,\tilde I(k_x)\,e^{-jk_{z0}h}\,\frac{k_y}{k_\rho},
\qquad
b_0^{d}=-\tfrac12\,\tilde I(k_x)\,e^{-jk_{z0}h}\,\frac{k_x}{k_\rho}$$

Inside the slab, with $\tau^2=e^{-2jk_{z1}d}$, the standard multiple-reflection sum gives

$$D=\frac{t^{01}X_0^{d}}{1+r^{01}r^{12}\tau^{2}},\qquad U = D\,r^{12}\tau^{2},\qquad
\tilde R^{01}=\frac{r^{01}+r^{12}\tau^{2}}{1+r^{01}r^{12}\tau^{2}}$$

(the last is the generalised reflection coefficient, [R2] §2.1). The field at depth $z$ is then

$$\tilde E_x = (a^d+a^u)\frac{k_y}{k_\rho}+Z_1^{\text{TM}}(b^d-b^u)\frac{k_x}{k_\rho},\qquad
\tilde E_y = -(a^d+a^u)\frac{k_x}{k_\rho}+Z_1^{\text{TM}}(b^d-b^u)\frac{k_y}{k_\rho}$$

$$\tilde E_z = (b^d+b^u)\,\frac{k_\rho}{\omega\varepsilon_1}$$

with $a^d(z)=D_a e^{+jk_{z1}z}$ and $a^u(z)=D_a\,r^{12}_a\,e^{-jk_{z1}(2d+z)}$. **That last form
matters numerically**: written naively as $U e^{-jk_{z1}z}$ it is a product of an underflowing and an
overflowing exponential in the evanescent spectrum, and produces NaN. Combining the exponents
analytically keeps every term decaying.

### 5.4 Antenna current

A thin centre-fed dipole carries, to good approximation, a sinusoidal current — Balanis [R7],
**Eq. (4-56)** (verified in the source text):

$$I_e(z')=\hat{\mathbf a}_zI_0\sin\!\Big[k\Big(\frac{\ell}{2}-|z'|\Big)\Big]$$

A monopole is half of this by image theory, so with the feed at $x=0$ and the open tip at $x=L$:

$$I(x)=I_m\sin\!\big[k_e(L-x)\big],\qquad 0\le x\le L$$

which correctly vanishes at the tip. $k_e$, the wavenumber of the current wave *on the antenna*, is
discussed in §6.3. Uniform and triangular distributions are also implemented for the sensitivity study.

---

## 6. Numerical implementation

`matlab/layered_field.m`. Approximately 180 lines; the whole calculation is 15 s.

### 6.1 The key idea

Because the layer stack is horizontally invariant, the field at each depth is a **2-D inverse Fourier
transform** of the spectral field. So instead of evaluating Sommerfeld integrals point by point
(with all their branch-cut and pole difficulties), I build the spectral field on a $k_x,k_y$ grid and
take one `fft2` per depth. This is exact — not an approximation — and it makes the finite 12 mm
antenna free: the current distribution enters only as the factor $\tilde I(k_x)$.

Discretisation: $f(x,y)=\texttt{fft2}(\tilde f)/(N_x\Delta x\,N_y\Delta y)$, with the $k$-grid in FFT
order. Baseline grid $N=1024$, $\Delta x=0.25$ mm (window ±128 mm), $\Delta z=0.25$ mm, 120 depths.

### 6.2 Branch and singularity handling

* $\mathrm{Im}\,k_{zn}\le0$ is **forced**, not left to the principal branch — for lossless air the
  principal square root gives the wrong (growing) sign in the evanescent region.
* $k_\rho=0$ is a removable, direction-dependent singularity of measure zero. It is guarded; the
  single affected sample out of $10^6$ contributes nothing.
* A standoff $h>0$ multiplies the evanescent spectrum by $e^{-|k_{z0}|h}$, which is what makes the
  spectral integral converge briskly.

### 6.3 The current-wave number $k_e$ — the one genuine unknown

$k_e$ is not determined by the problem statement, and it is the largest structural assumption.
Three defensible hypotheses:

| | hypothesis | $\varepsilon_\text{eff}$ | $k_e$ (rad/m) | $L/\lambda_e$ | implied $R_\text{in}$ |
|---|---|---|---|---|---|
| **H1** | the 12 mm element is a **resonant quarter wave** | 6.5 | 130.9 | 0.250 | 18.6 Ω |
| H2 | bare wire on the interface, $\varepsilon_\text{eff}=(1+\varepsilon')/2$ | 37.6 | 315 | 0.602 | 41.9 Ω |
| H3 | the current wave sees air only | 1 | 51.3 | 0.098 | 12.9 Ω |

**H1 is the baseline**, for two reasons. First, practical 2.45 GHz ablation monopoles are 10–15 mm
because the conductor is *insulated* (coax dielectric, catheter wall), and King & Smith's insulated-
antenna theory [R3] shows the insulation — not the outer medium — sets the guided wavelength; a bare
quarter-wave in this gel would be only 3.5 mm. Second, H1 is the hypothesis consistent with your
generator spec: $P_\text{ref}/P_\text{fwd}\le 50/200$ gives $|\Gamma|\le 0.50$, and H1's
$R_\text{in}=18.6\ \Omega$ against 50 Ω gives $|\Gamma|=0.458$.

**How much weight to put on that second argument: not much.** It ignores $X_\text{in}$ entirely
(a real antenna has substantial reactance, which only increases $|\Gamma|$), it ignores conductor
loss, and any matching network in your chain invalidates it. Treat it as a plausibility check that
H1 is not absurd — not as a measurement. **The good news (§9A) is that it barely matters.**

### 6.4 Absolute scale from delivered power

You do not know the antenna's field strength, and you do not need to. Compute the field for a
reference current $I_m=1$ A, then close Poynting's theorem over the whole structure:

$$P_\text{total}= \underbrace{\int_\text{gel}\tfrac12\sigma_\text{eff}|\tilde{\mathbf E}|^{2}dV}_{P_\text{gel}}
\;+\;\underbrace{P_\uparrow}_{\text{radiated into air}}\;+\;\underbrace{P_\downarrow}_{\text{out of the floor}}$$

The two flux terms are evaluated in the spectral domain by Parseval's theorem — for a single
polarisation pair, $(\mathbf E\times\mathbf H^{*})_z = E_sH_p^{*}-E_pH_s^{*}$, and the normalisation
constant is identical to the FFT one. Evanescent components in lossless air contribute exactly zero
real power, as they must. Then

$$I_m = \sqrt{P_\text{del}/P_\text{total}(I_m{=}1)},\qquad |\mathbf E|\propto\sqrt{P_\text{del}}$$

For the baseline: $I_m = 4.02$ A, $R_\text{in}=2P_\text{total}/|I_\text{feed}|^2 = 18.6\ \Omega$.

**Power budget at 150 W delivered:** gel **149.58 W (99.7 %)**, radiated up into air 0.39 W (0.3 %),
out through the floor 0.03 W (<0.1 %). The gel is an extremely efficient absorber — which is the
single most useful sanity check on the whole model, and it means the answer is essentially
independent of what happens above and below.

---

## 7. Validation

`matlab/verify_layered_field.m`. Four independent checks.

**7.1 Homogeneous limit vs an analytic solution.** Setting all three layers to the same lossy medium
removes every interface, so the solver must reproduce the exact dyadic-Green's-function field of an
$x$-directed Hertzian dipole, $\mathbf E=-j\omega\mu\,I\,d\ell\,[\bar{\mathbf I}+\nabla\nabla/k^2]G\cdot\hat{\mathbf x}$
(coded independently in `hertzian_dipole_field.m`). Result:

| cut | max relative error | rms |
|---|---|---|
| vertical line $x=y=0$ | $1.4\times10^{-4}$ | $5.2\times10^{-5}$ |
| line along $x$, $z=-5$ mm (tests TM) | $1.3\times10^{-3}$ | $7.4\times10^{-4}$ |
| line along $y$, $x=0$ (tests TE) | $2.4\times10^{-3}$ | $1.4\times10^{-3}$ |
| diagonal $x=y$ (tests TE/TM mixing) | $2.4\times10^{-3}$ | $1.4\times10^{-3}$ |

This exercises the Weyl expansion, the polarisation projection, the FFT normalisation and the branch
choice simultaneously. It is the strongest check available.

**7.2 Grid and window convergence.**

| $N$ | $\Delta x$ | window | $P_\text{gel}$ | deviation in $|E|$ |
|---|---|---|---|---|
| 1024 | 0.25 mm | ±128 mm | 7.45 W | (reference) |
| 2048 | 0.25 mm | ±256 mm | 7.45 W | 0.14 % |
| 2048 | 0.125 mm | ±128 mm | 7.45 W | 0.00 % |
| 1536 | 0.20 mm | ±154 mm | 7.45 W | 0.25 % |

Converged to better than 0.3 % — discretisation is not a limiting error.

**7.3 Attenuation.** The deep-field decay (94.2 Np/m over 12→22 mm) correctly exceeds the plane-wave
$\alpha=58.98$ Np/m, the excess being geometric spreading.

**7.4 An unforced physical check.** The depth profiles (`fig5`) show a standing-wave ripple with
minima 7.0 mm apart. The predicted half-wavelength in the gel is
$\lambda_\text{gel}/2 = 14.07/2 = 7.04$ mm. Nothing in the model was tuned to produce that; it falls
out of the reflection at the gel/air floor, and it matches to 0.5 %.

---

## 8. Results

Figures in `figures/`:

| file | content |
|---|---|
| `fig1_Erms_dB.png` | **the main deliverable** — $E_\text{rms}$ on the $y=0$ plane, dB re max |
| `fig2_Erms_linear.png` | same, linear V/m with labelled contours |
| `fig3_SAR.png` | SAR map |
| `fig4_components.png` | $|E_x|,|E_y|,|E_z|$ separately |
| `fig5_depth_profiles.png` | depth profiles at 8 values of $x$ |
| `fig6_time_snapshots.png` | $|\mathbf E(\mathbf r,t)|$ at four phases |
| `fig7_axial_ratio.png` | polarisation ellipticity |
| `fig8_E_vs_time.png` | $|\mathbf E|$ vs $t$ at three points — shows the $2f$ ripple |
| `fig9_closed_form_fit.png` | numerical vs the three analytic approximations |

**Shape of the deposition.** The energy is confined to a lobe roughly $-5 \le x \le 17$ mm wide
(i.e. the antenna footprint plus about one $\lambda_\text{gel}/3$ at each end) reaching ~20 mm deep,
with a weak secondary bulge at 20–25 mm from the floor reflection. The peak is at the **feed** end,
not the centre, because the quarter-wave current distribution peaks there.

**Penetration on the vertical through the antenna centre ($x=6$ mm):**

| metric | value |
|---|---|
| $E_\text{rms}$ at 0.125 mm depth | 17 374 V/m |
| depth at which $E$ falls to $1/e$ | **2.66 mm** |
| depth at which $E$ falls to $1/10$ | **12.09 mm** |
| $E$ at the floor (30 mm) | 524 V/m = 3.0 % of surface |
| share of in-plane power deposited above 10 mm | **96.4 %** |

Note the near-surface $1/e$ distance (2.7 mm) is **far shorter than the plane-wave skin depth**
(16.9 mm). Near the antenna the decay is dominated by reactive near-field spreading, not by material
attenuation; the two only become comparable beyond ~15 mm. Quoting the skin depth as "the
penetration" would overstate the heated depth by roughly 6×.

**Polarisation.** The field is strongly elliptical in the near zone (axial ratio 0.5–0.9 within a few
mm) and becomes nearly linear at depth (0.02 at 20 mm). Practically: close to the antenna
$|\mathbf E(t)|$ ripples only modestly and never reaches zero; deep in the gel it swings almost fully
between 0 and $E_\text{peak}$ at 4.9 GHz.

---

## 9. Sensitivity analysis

`matlab/run_sensitivity.m` → `results/sensitivity_log.txt`. $E$ values are $E_\text{rms}$ in V/m at
$x=6$ mm; $z_{1/e}$ is the depth at which $E$ falls to $1/e$ of its near-surface value.

| case | $R_\text{in}$ (Ω) | % in gel | E@2 mm | E@10 mm | E@25 mm | $z_{1/e}$ (mm) |
|---|---|---|---|---|---|---|
| **A baseline: sinusoid, $L=\lambda_e/4$** | 18.6 | 99.7 | **7639** | **2337** | **431** | 2.66 |
| A sinusoid, bare wire on interface | 41.9 | 99.8 | 8294 | 2227 | 374 | 2.97 |
| A sinusoid, wire sees air only | 12.9 | 99.7 | 7590 | 2258 | 419 | 2.36 |
| A uniform current | 38.8 | 99.7 | 7016 | 2451 | 464 | 4.01 |
| A triangular current | 12.2 | 99.7 | 7578 | 2244 | 417 | 2.32 |
| B standoff $h$ = 0.10 mm | 33.6 | 99.9 | 6609 | 1877 | 337 | 0.84 |
| B standoff $h$ = 0.25 mm | 24.3 | 99.8 | 7328 | 2144 | 389 | 1.51 |
| B standoff $h$ = 1.00 mm | 13.4 | 99.6 | 7595 | 2503 | 476 | 4.40 |
| B standoff $h$ = 2.00 mm | 8.4 | 99.2 | 7166 | 2640 | 533 | 6.57 |
| B standoff $h$ = 5.00 mm | 3.3 | 95.0 | 5785 | 2605 | 608 | 9.90 |
| C 1.1 % NaCl, 25 °C | 18.8 | 99.7 | 7498 | 2162 | 359 | 2.59 |
| C 0.9 % NaCl, 20 °C | 18.5 | 99.7 | 7623 | 2333 | 431 | 2.69 |
| C 0.9 % NaCl, 37 °C | 18.8 | 99.7 | 7599 | 2251 | 390 | 2.55 |
| C muscle-equivalent target | 18.8 | 99.6 | 8127 | 2580 | 466 | 2.09 |
| C pure water, no salt | 17.7 | 99.7 | 8358 | 3395 | 1008 | 3.16 |
| D gel depth 20 mm | 18.2 | 99.6 | 7604 | 2187 | n/a¹ | 2.59 |
| D gel depth 40 mm | 18.5 | 99.7 | 7635 | 2327 | 502 | 2.65 |
| E metal (PEC) dish under gel | 18.4 | 99.7 | 7629 | 2310 | 580 | 2.64 |

¹ 25 mm lies outside a 20 mm slab; the tabulated value is the nearest in-range depth and is not comparable.

**A. Current distribution — the reassuring result.** Across *five* quite different assumptions,
$E$@2 mm spans 7016–8294 V/m (**±8 %**) and $E$@10 mm spans 2227–2451 (**±5 %**). Meanwhile
$R_\text{in}$ swings by a factor of 3.4. **Power normalisation is doing the work**: once the total
delivered power is fixed, the *shape* of the current has surprisingly little effect on the field in
the gel, because the gel absorbs 99.7 % of it regardless. This substantially de-risks the one
assumption I was most worried about (§6.3), and it means a full method-of-moments solution for the
current would change the answer by less than the material uncertainty.

**B. Standoff — the dominant sensitivity.** Lifting the wire from 0.1 mm to 5 mm changes $z_{1/e}$
from 0.84 mm to 9.90 mm (12×) and E@10 mm by +39 %. Physically it trades peak surface intensity for
depth uniformity. **If you want the model to be quantitatively right, measure this gap.** It is also
your main experimental control over deposition shape.

**C. Material.** Temperature (20–37 °C) is negligible (<5 %). Salt concentration matters at depth:
1.1 % vs 0.9 % changes E@25 mm by −17 %, and the salt-free case is +134 % — so **the NaCl content is
the dominant material variable for penetration depth**, exactly as expected from
$\sigma_i \propto N$. Moving to the muscle-equivalent target changes things by only 6–10 %, so the
glycerol question of §4.4 is *not* a first-order concern for the field map (it matters for whether
the phantom is a valid tissue surrogate — a different question).

**D/E. Boundaries.** The floor is irrelevant above ~20 mm depth. Below that, a metal dish raises
E@25 mm by 35 % over an air-backed dish. Slab thickness beyond 30 mm changes nothing measurable.

---

## 10. The final function — the deliverable

### 10.1 Derivation of the time dependence

With $\tilde E_i = |\tilde E_i|e^{j\phi_i}$ and $E_i(t)=|\tilde E_i|\cos(\omega t+\phi_i)$,

$$|\mathbf E(t)|^{2}=\sum_i|\tilde E_i|^{2}\cos^{2}(\omega t+\phi_i)
=\sum_i|\tilde E_i|^{2}\frac{1+\cos(2\omega t+2\phi_i)}{2}
=\tfrac12\Big[\big|\tilde{\mathbf E}\big|^{2}+\mathrm{Re}\big\{Q\,e^{j2\omega t}\big\}\Big]$$

since $\sum_i|\tilde E_i|^2e^{j2\phi_i}=\sum_i\tilde E_i^{2}=Q$. Hence

$$|\mathbf E(x,z,t)|=\sqrt{\tfrac12\Big[\big|\tilde{\mathbf E}\big|^{2}+|Q|\cos(2\omega t+\arg Q)\Big]}$$

Checks: for linear polarisation $|Q|=\big|\tilde{\mathbf E}\big|^2$ and this reduces to
$\big|\tilde{\mathbf E}\big|\,|\cos(\omega t+\phi)|$; for circular polarisation $Q=0$ and $|\mathbf E|$ is
constant. The period of $|\mathbf E|$ is $\pi/\omega = 1/(2f) =$ **204 ps** (the drive period is 408 ps).

### 10.2 Exact evaluation — `E_gel.m`

```matlab
addpath('matlab');
run_model                                   % once: builds results/field_model.mat

Emag        = E_gel(6e-3, -5e-3, 0);        % |E| at x=6mm, 5mm deep, t=0   [V/m]
Erms        = E_gel(6e-3, -5e-3);           % omit t -> time-RMS value
[Emag, S]   = E_gel(6e-3, -5e-3, 0);        % S: components, peak, min,
                                            %    axial ratio, SAR, phases
t  = linspace(0, 2/2.45e9, 500);
Et = E_gel(6e-3, -5e-3, t);                 % the waveform
```

Arguments are in **SI units (metres, seconds)**; $z$ is **negative** inside the gel. Inputs are
broadcast, so `E_gel(X, Z, t)` over meshgrids returns a whole map. Verified output at
$(6\,\text{mm},\,-5\,\text{mm})$: $E_\text{rms}=4320$, $E_\text{peak}=5871$, $E_\text{min}=1690$ V/m,
axial ratio 0.288, SAR 50 800 W/kg.

### 10.3 Compact analytic approximation

If you need a formula rather than a lookup, `fit_closed_form.m` fits three nested ansatz. The
physically correct one — a **line** source (the antenna is 12 mm, not a point) plus its image in the
gel floor — is:

$$\boxed{\;E_\text{rms}(x,z)\;\simeq\;K\left|\int_0^{L}\!\sin\!\big[k_e(L-s)\big]\Big[g(R_d,\psi_d)e^{-jk_1R_d}+\Gamma\,g(R_r,\psi_r)e^{-jk_1R_r}\Big]ds\right|\;}$$

$$g(R,\psi)=(\sin\psi)^{q}R^{-p},\qquad
R_d=\sqrt{(x-s)^2+(z-h)^2},\quad \sin\psi_d=\frac{|z-h|}{R_d}$$

$$R_r=\sqrt{(x-s)^2+(z+2d+h)^2},\quad \sin\psi_r=\frac{|z+2d+h|}{R_r}$$

| parameter | value | note |
|---|---|---|
| $K$ | 3924 | at $P_\text{del}=150$ W; scales as $\sqrt{P_\text{del}}$ |
| $q$ | 0.784 | pattern exponent (1 = ideal dipole) |
| $p$ | 1.144 | spreading exponent (1 = spherical) |
| $\Gamma$ | $0.170-0.253j$, $|\Gamma|=0.305$ | effective floor reflection |
| $k_1$ | $446.49-j\,58.98$ rad/m | **exact, not fitted** |
| $k_e$, $L$, $h$, $d$ | 130.9 rad/m, 12 mm, 0.5 mm, 30 mm | geometry |

**Accuracy** over $4\le R_d\le 28$ mm, depth ≥ 1 mm, $|x|\le45$ mm (18 142 points):
**1.41 dB rms, 70 % within 1 dB, 90 % within 2 dB**, worst case 11.2 dB (in the deep standing-wave
nulls and at extreme end-fire). For comparison the naive point-source forms give 3.36 and 3.33 dB
rms — the line-source structure is what buys the accuracy, because the current vanishes at the tip
and a point source badly over-predicts the field beyond $x=12$ mm.

Use this for quick estimates; use `E_gel.m` when you need the answer right.

---

## 11. Limitations — read before quoting numbers

1. **The standoff $h$ dominates.** ±1 mm of uncertainty here moves the near-surface field more than
   everything else combined (§9B). Measure it.
2. **The current distribution is prescribed, not solved.** A method-of-moments solution in the
   layered medium would remove this. §9A shows it would change $|E|$ by <8 %, so it is a refinement,
   not a correction — but it is the honest next step, and it would also give $X_\text{in}$ and hence
   a real prediction of $|\Gamma|$ to compare against your generator.
3. **Coax sheath currents are not modelled.** A real coax-fed monopole carries current on the outer
   conductor, which is the known cause of "backward heating" along the catheter in the ablation
   literature [R5,R6]. The model radiates only from the 12 mm element, so it will *under-predict*
   heating behind the feed ($x<0$).
4. **The Phytagel permittivity is inferred, not measured.** It is the saline value from a validated
   model (§4.2). No published measurement of this exact gel exists. **If you can get 20 minutes on an
   open-ended coaxial dielectric probe, one measurement removes the largest material uncertainty in
   the model.** And see §4.4 on the unverified glycerol content.
5. **Constant properties.** Both $\varepsilon'$ and $\sigma$ of aqueous media fall substantially as
   temperature rises; at these SAR levels (>10⁵ W/kg near the antenna) the gel heats in
   milliseconds, so a thermal run would need coupled EM–thermal iteration. The model is valid at
   $t=0^+$ / low duty cycle.
6. **Ideal geometry.** Perfectly flat surface, infinite lateral extent, no container walls, no
   irrigation, no air gap from surface curvature.
7. **The $R_\text{in}$/$|\Gamma|$ argument in §6.3 is a plausibility check, not evidence** — it
   ignores reactance, conductor loss and any matching network.

---

## 12. Running the code

```
matlab/
  saline_stogryn.m        Stogryn 1971 saline permittivity          [R8,R10]
  tissue_cole_cole.m      Gabriel 4-pole Cole-Cole tissue model     [R11,R12]
  model_config.m          every constant + geometry, each sourced
  layered_field.m         the solver (spectral 3-layer, FFT)        [R1,R2]
  hertzian_dipole_field.m analytic reference for validation         [R7,R16]
  verify_layered_field.m  the validation suite
  run_model.m             baseline run + 8 figures + results/*.mat
  run_sensitivity.m       18-case parameter study
  fit_closed_form.m       compact analytic reduction
  E_gel.m                 *** THE DELIVERABLE FUNCTION ***
```

```matlab
addpath('matlab');
verify_layered_field     % ~2 min  - validation (run this first)
run_model                % ~30 s   - baseline + figures
run_sensitivity          % ~5 min  - parameter study
fit_closed_form          % ~3 min  - analytic reduction
```

Tested on MATLAB R2025b. No toolboxes beyond base MATLAB are required.

---

## References

**[R1]** K. A. Michalski and J. R. Mosig, "Multilayered media Green's functions in integral equation
formulations," *IEEE Trans. Antennas Propag.*, vol. 45, no. 3, pp. 508–519, 1997.
doi:[10.1109/8.558666](https://doi.org/10.1109/8.558666) ·
[IEEE Xplore](https://ieeexplore.ieee.org/document/558666/)

**[R2]** W. C. Chew, *Waves and Fields in Inhomogeneous Media*. New York: IEEE Press, 1995. §2.1–2.2
(Weyl identity; generalised reflection coefficients for layered media).

**[R3]** R. W. P. King and G. S. Smith, *Antennas in Matter: Fundamentals, Theory, and Applications*.
Cambridge, MA: MIT Press, 1981. ISBN 0-262-11074-1.
[Internet Archive](https://archive.org/details/antennasinmatter00rono)

**[R4]** M. Gao, Z. Šipuš and A. K. Skrivervik, "Approximate Expressions to the Reactive Near-Field
Losses of Body-Implanted Antennas," arXiv:2205.04971, 2022.
doi:[10.48550/arXiv.2205.04971](https://doi.org/10.48550/arXiv.2205.04971) —
planar muscle phantom, Green's-function E-field map at 2.45 GHz.

**[R5]** P. Bernardi, M. Cavagnaro and J. C. Lin, "Distribution of SAR and temperature elevation
induced in a phantom by a microwave cardiac ablation catheter," *IEEE Trans. Microw. Theory Tech.*,
2004. [IEEE Xplore](https://ieeexplore.ieee.org/document/1318794/) —
2450 MHz cap-choke catheter antenna in a homogeneous dielectric phantom.

**[R6]** "In vitro and in vivo evaluation of the thermal patterns and lesions of catheter ablation
with a microwave monopole antenna," *J. Interv. Card. Electrophysiol.*, 2004.
doi:[10.1023/B:JICE.0000019264.14661.8e](https://doi.org/10.1023/B:JICE.0000019264.14661.8e) ·
[PubMed 15014211](https://pubmed.ncbi.nlm.nih.gov/15014211/) —
2.45 GHz monopole in contact with a myocardium-equivalent phantom surface.

**[R7]** C. A. Balanis, *Antenna Theory: Analysis and Design*, 3rd ed. Hoboken, NJ: Wiley, 2005.
Eq. (4-56) (sinusoidal current distribution, p. 170); §4.2, §3.6.
[full text (verified against, Eq. 4-56)](https://jontalle.web.engr.illinois.edu/Public/BOOKS/Balanis-AntennaTheoryAnalysis_Design3rdEd.05.pdf)

**[R8]** A. Stogryn, "Equations for calculating the dielectric constant of saline water," *IEEE Trans.
Microw. Theory Tech.*, vol. MTT-19, no. 8, pp. 733–736, 1971.
doi:[10.1109/TMTT.1971.1127617](https://doi.org/10.1109/TMTT.1971.1127617)

**[R9]** L. A. Klein and C. T. Swift, "An improved model for the dielectric constant of sea water at
microwave frequencies," *IEEE Trans. Antennas Propag.*, vol. AP-25, no. 1, pp. 104–111, 1977.
[PDF](https://aquarius.oceansciences.org/docs/ed_klien_swift_1978.pdf) — source of $\varepsilon_\infty=4.9$.

**[R10]** D. H. Gadani, V. A. Rana, S. P. Bhatnagar, A. N. Prajapati and A. D. Vyas, "Effect of
salinity on the dielectric properties of water," *Indian J. Pure Appl. Phys.*, vol. 50, pp. 405–410,
2012. [open-access PDF](https://nopr.niscpr.res.in/bitstream/123456789/14183/1/IJPAP%2050(6)%20405-410.pdf)
— the open reproduction of the Stogryn equations actually transcribed into the code.

**[R11]** C. Gabriel, "Compilation of the Dielectric Properties of Body Tissues at RF and Microwave
Frequencies," Report AL/OE-TR-1996-0037, Brooks AFB, TX, 1996.
[report text](https://niremf.ifac.cnr.it/docs/DIELECTRIC/Report.html)

**[R12]** S. Gabriel, R. W. Lau and C. Gabriel, "The dielectric properties of biological tissues: III.
Parametric models for the dielectric spectrum of tissues," *Phys. Med. Biol.*, vol. 41, no. 11,
pp. 2271–2293, 1996. doi:[10.1088/0031-9155/41/11/003](https://doi.org/10.1088/0031-9155/41/11/003).
Muscle parameter column transcribed from the open reproduction in *An Analytical Framework for
Frequency-Dependent Electromagnetic Power Absorption in Biological Tissues*,
[arXiv:2604.06225](https://arxiv.org/pdf/2604.06225), Table 1; verified by reproducing the IT'IS
value at 2.45 GHz (§4.5).

**[R13]** P. C. Qian et al., "Irrigated Microwave Catheter Ablation Can Create Deep Ventricular
Lesions Through Epicardial Fat With Relative Sparing of Adjacent Coronary Arteries,"
*Circ. Arrhythm. Electrophysiol.*, 2020.
doi:[10.1161/CIRCEP.119.008251](https://doi.org/10.1161/CIRCEP.119.008251) —
**source of the Phytagel phantom recipe; recipe text NOT verified against the PDF (§4.4).**

**[R14]** P. C. Qian et al., "Transcatheter non-contact microwave ablation may enable circumferential
renal artery denervation while sparing the vessel intima and media," *EuroIntervention*.
[link](https://eurointervention.pcronline.com/article/transcatheter-non-contact-microwave-ablation-may-enable-circumferential-renal-artery-denervation-while-sparing-the-vessel-intima-and-media)
— **same caveat.**

**[R15]** W. Zhang, D. Luan, J. Tang, S. S. Sablani, B. Rasco, H. Lin and F. Liu, "Dielectric
properties and other physical properties of low-acyl gellan gel as relevant to microwave assisted
pasteurization process," *J. Food Eng.*, vol. 149, pp. 195–203, 2015.
doi:[10.1016/j.jfoodeng.2014.10.017](https://doi.org/10.1016/j.jfoodeng.2014.10.017) —
1 % gellan with 0–300 mM NaCl; paywalled, numeric tables not obtained.

**[R16]** C. A. Balanis, *Advanced Engineering Electromagnetics*, 2nd ed. Wiley, 2012. Ch. 6
(dyadic Green's functions), used for the analytic Hertzian-dipole reference in §7.1.

**[R17]** IT'IS Foundation, *Tissue Properties Database — Dielectric Properties*.
[itis.swiss](https://itis.swiss/virtual-population/tissue-properties/database/dielectric-properties)
— used only as the cross-check target for §4.5.
