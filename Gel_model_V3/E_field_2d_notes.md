# E-field in the gel plane — 27 mm monopole, homogeneous lossy dielectric

Date: 2026-09-14 · Code: `E_field_2d.m` · Checks: `verify_E_field_2d.m` · Constants: `Experimental_setup_table.xlsx`

Same chain as the earlier iterations (`../E_field_model_writeup.md`: sinusoidal wire current → superposed
current elements → scale from delivered power → time dependence), with the three-layer Sommerfeld
solve replaced by the closed-form field of a current element in **one** unbounded lossy dielectric.
Five equations, nothing else.

## 0. What you get

```matlab
Et = E_field_2d(5e-3, -5e-3, deg2rad(0:45:315), 100);   % x, z [m], phase [rad], P_set [W]
Et.x, Et.z                                               % E_x, E_z [V/m] at each phase
[Et, S] = E_field_2d(x, z, phase, P_set);                % S: phasor, E_rms, E_peak, E_min, q [W/m^3]
```

Antenna frame: feed at $x=0$, tip at $x=27$ mm, wire axis $h=0.5$ mm above the surface, gel is
$z\le0$, plane $y=0$. $x,z$ may be arrays (a whole map). $|\mathbf E|\propto\sqrt{P_\text{del}}$
exactly, so one call covers every power.

## 1. The five equations

**(1) Helmholtz equation in the lossy dielectric.** Time convention $e^{+j\omega t}$.

$$\nabla^2\tilde{\mathbf E}+k_c^2\tilde{\mathbf E}=0,\qquad
k_c=\omega\sqrt{\mu_0\varepsilon_0\,(\varepsilon'-j\varepsilon'')},\qquad \mathrm{Im}\,k_c\le0$$

Source: Ellingson, *Electromagnetics Vol. 2*, [§3.3 Wave Equations for Lossy Regions](https://phys.libretexts.org/Bookshelves/Electricity_and_Magnetism/Electromagnetics_II_(Ellingson)/03%3A_Wave_Propagation_in_General_Media/3.03%3A_Wave_Equations_for_Lossy_Regions)
($\nabla^2\tilde{\mathbf E}-\gamma^2\tilde{\mathbf E}=0$, $\gamma^2=-\omega^2\mu\varepsilon_c$, i.e. $k_c=-j\gamma$) and
[§3.4 Complex Permittivity](https://phys.libretexts.org/Bookshelves/Electricity_and_Magnetism/Electromagnetics_II_(Ellingson)/03%3A_Wave_Propagation_in_General_Media/3.04%3A_Complex_Permittivity)
($\varepsilon_c=\varepsilon'-j\varepsilon''$; a *measured* $\varepsilon''$ already contains both conduction and relaxation loss).
Why it applies: the gel is a lossy dielectric, $\tan\delta=13.5/43.98=0.31\ll1$ is not a conductor, and
the drive is single-frequency, so the phasor form is exact. With the table's gel values:

| $\varepsilon_r$ | $k_c$ (1/m) | $\lambda_\text{gel}$ | field $1/e$ depth $1/\alpha$ | power $1/e$ depth | $\sigma_\text{eff}=\omega\varepsilon_0\varepsilon''$ | $L/\lambda_\text{gel}$ |
|---|---|---|---|---|---|---|
| $43.98-j13.5$ | $344.4-j51.7$ | 18.24 mm | 19.35 mm | 9.7 mm | 1.840 S/m | 1.48 |

**(2) Current on the monopole.** Balanis, *Antenna Theory* 3rd ed. [R7], **Eq. (4-56)** (thin linear
antenna, sinusoidal current; the monopole is one arm of it, §4.7 image theory — as in
`../E_field_model_writeup.md` §5.4):

$$I(x')=I_m\sin\!\big[k_e(L-x')\big],\ 0\le x'\le L;\qquad
I(x')=I(0)\cos\!\Big(\frac{\pi x'}{2L_s}\Big),\ -L_s\le x'<0$$

$k_e=k_c$: a bare wire immersed in one homogeneous medium carries the wavenumber of that medium
(King & Smith, *Antennas in Matter* [R3], bare dipole in a dissipative medium). With complex $k_c$ the
current is largest at the feed and decays $\sim e^{-\alpha L}=1/4$ toward the tip, which is what a
1.5 $\lambda$ wire in a lossy medium does. The second line is the return current on the bare coax
braid, `../v2_water_above/README.md` §1, kept because truncating $I$ at the feed puts a fictitious
point charge there; $L_s=L$ as in v2.1, `'L_sheath',0` removes it.

**(3) Field of the wire.** The solution of (1) driven by the line current is the superposition of
$x$-directed current elements $I(x')\,ds$ in the unbounded medium, Balanis [R7] §4.2
**Eqs. (4-8a)–(4-10c)** — exact for any $r$, and valid for complex $k$ because nothing in the derivation
assumes $k$ real:

$$\tilde{\mathbf E}(x,z)=-j\omega\mu_0\sum_{s} I(x'_s)\,ds\,\Big[\mathbf I+\frac{\nabla\nabla}{k_c^2}\Big]\frac{e^{-jk_cR_s}}{4\pi R_s}\cdot\hat{\mathbf x},
\qquad R_s=\sqrt{(x-x'_s)^2+(z-h)^2}$$

The kernel is `../matlab/hertzian_dipole_field.m` **verbatim** — the analytic reference that validated
the layered model to $10^{-4}$ (writeup §7). On $y=0$ only $E_x$ and $E_z$ exist. Segment length
$ds=0.1$ mm (halving it changes nothing to 4 digits).

**(4) Scale from delivered power.** `../E_field_model_writeup.md` §6.4: compute for $I_m=1$ A, then

$$P_\text{del}=\eta P_\text{set}=\int_{\rho\ge a}\tfrac12\,\sigma_\text{eff}\,|\tilde{\mathbf E}|^2\,dV,
\qquad I_m=\sqrt{P_\text{del}/P(I_m{=}1)}$$

With no reflection and an unbounded lossy medium, everything the antenna puts out is absorbed, so
the dissipation integral **is** the delivered power — no flux terms this time. The integral is
axisymmetric about the wire ($dV=2\pi\rho\,d\rho\,dx$) and starts at the conductor radius
$a=0.45$ mm, because $|\mathbf E|\sim1/\rho$ next to a bare wire and the log-divergence is cut off by
the conductor itself. $\eta=0.2$ from your statement ("runs at 20 % of the set power"); the table's
"power dropoff 0.7" is one edit if it meant 30 % delivered.

**(5) The oscillation.** `../E_field_model_writeup.md` §10.1, unchanged:

$$\mathbf E(x,z,t)=\mathrm{Re}\{\tilde{\mathbf E}\,e^{j\omega t}\}
=|\tilde E_x|\cos(\omega t+\phi_x)\,\hat{\mathbf x}+|\tilde E_z|\cos(\omega t+\phi_z)\,\hat{\mathbf z},
\qquad
|\mathbf E(t)|=\sqrt{\tfrac12\big[|\tilde{\mathbf E}|^2+|Q|\cos(2\omega t+\arg Q)\big]},\ Q=\tilde E_x^2+\tilde E_z^2$$

Period 408 ps; the magnitude ripples at $2f$ between $E_\text{peak}$ and $E_\text{min}$, the
semi-axes of the polarisation ellipse. For heating only the cycle mean matters:
$q=\tfrac12\sigma_\text{eff}|\tilde{\mathbf E}|^2=\sigma_\text{eff}E_\text{rms}^2$ (`S.q`).

## 2. Assumptions, and why each is acceptable here

| assumption | justification / cost |
|---|---|
| One homogeneous medium with the **gel's** $\varepsilon_r$, the water above ignored | The point of this iteration. Water $73-j8.8$ vs gel $43.98-j13.5$: normal-incidence $\|\Gamma\|=0.125$, i.e. 1.6 % of power reflected and the transmitted amplitude is $1+\Gamma\approx1.12$ — a ~10 % error on $\|E\|$ near the surface. The antenna's own surroundings are water, so the true current wavelength is $8.56/6.71=1.28\times$ longer than assumed; `'ke'` is exposed for that. |
| Bare thin wire, $a=0.45$ mm, axis $h=0.5$ mm above the surface | Same geometry as `../matlab/model_config.m`. Insulation would lower the near-field loss and lengthen the current wave; not in the table, so not modelled. |
| $k_e=k_c$ | King & Smith [R3]. The earlier quarter-wave hypothesis (H1, writeup §6.3) was argued from *insulation* of a 12 mm element; a 27 mm wire in one dielectric has no such argument. Sensitivity at (5, −5) mm: $E_\text{rms}$ 617 V/m baseline, 808 with $k_e=\pi/2L$, 745 without the braid arm. |
| Monopole equations only | (2) and (3) are the thin-wire sinusoidal-current model of a linear monopole/dipole; no slot, choke, helical or aperture formulas are used. |
| Microwave, not RF ablation | This is a radiating problem — (1) with $k_c\ne0$, heating by the **total** dielectric loss $\varepsilon''$ (ionic + dipolar). RF ablation would be the quasi-static Laplace problem for a conduction current; none of that is used. |
| $\varepsilon''=13.5$ is the total loss factor | Ellingson §3.4. The table's separate $\sigma=1.7$ S/m corresponds to $\varepsilon''=12.5$, 8 % lower; `'epsr'` is one edit. |
| No reflection, $\eta=0.2$ | User statement. Every $\|E\|$ scales as $\sqrt{P_\text{del}}$, so 10–30 W delivered is 0.71–1.22 × the 20 W numbers below. |

## 3. Verification (`verify_E_field_2d.m`)

| check | result |
|---|---|
| field = Σ `hertzian_dipole_field.m` over the same current | identical (rel err 0) |
| energy: Poynting flux out of a closed cylinder round the wire = dissipation outside it | 117.5 W vs 118.1 W, **−0.5 %** (quadrature) |
| scale integral $P(I_m{=}1)=156.6$ W: $ds/2$, grid ×2, extent 5/α | 0.00 %, −0.01 %, +0.15 % |
| sampled $\|\mathbf E(t)\|$ vs the closed form of (5) | rel err $5\times10^{-16}$ |

## 4. Result at $P_\text{set}=100$ W ($P_\text{del}=20$ W), $I_m=0.357$ A, $|I_\text{feed}|=0.678$ A

| $x$ (mm) | $z$ (mm) | $\|\tilde E_x\|$ | $\|\tilde E_z\|$ | $\phi_x$ | $\phi_z$ | $E_\text{rms}$ (V/m) | $q$ (W/m³) |
|---|---|---|---|---|---|---|---|
| 5 | −2 | 808 | 1658 | 170.8° | 156.6° | 1304 | 3.13e6 |
| 5 | −5 | 609 | 624 | 137.1° | 117.7° | 617 | 7.00e5 |
| 5 | −10 | 363 | 252 | 58.9° | 34.9° | 313 | 1.80e5 |
| 13.5 | −5 | 255 | 401 | −10.8° | −18.5° | 336 | 2.08e5 |
| 25 | −5 | 257 | 251 | 154.8° | 116.2° | 254 | 1.19e5 |
| −5 (under the braid) | −5 | 1200 | 75 | −155.0° | 167.7° | 850 | 1.33e6 |

The waveform at (5, −5) mm, $E_\text{peak}=860$, $E_\text{min}=147$ V/m:

| phase $\omega t$ | 0° | 45° | 90° | 135° | 180° | 225° | 270° | 315° |
|---|---|---|---|---|---|---|---|---|
| $E_x$ (V/m) | −446 | −609 | −415 | 22 | 446 | 609 | 415 | −22 |
| $E_z$ (V/m) | −290 | −596 | −553 | −186 | 290 | 596 | 553 | 186 |

## References

[R3] R. W. P. King, G. S. Smith, *Antennas in Matter*, MIT Press 1981 — [archive.org](https://archive.org/details/antennasinmatter00rono).
[R7] C. A. Balanis, *Antenna Theory*, 3rd ed., Wiley 2005: Eq. (4-56); §4.2 Eqs. (4-8a)–(4-10c); §4.7 — [full text](https://jontalle.web.engr.illinois.edu/Public/BOOKS/Balanis-AntennaTheoryAnalysis_Design3rdEd.05.pdf).
S. W. Ellingson, *Electromagnetics Vol. 2*, LibreTexts, §3.3–3.4 (links above).
Internal: `../E_field_model_writeup.md` §5.4, §6.4, §7, §10.1; `../v2_water_above/README.md` §1; `../matlab/hertzian_dipole_field.m`; `Experimental_setup_table.xlsx`.
