# References — Dielectric properties of 0.9–1.1 % saline Phytagel

**Status of the search (as of 2026-09-05):** no published source states, in exact words, both the
real and imaginary permittivity of a 0.9–1.1 % saline Phytagel (gellan gum) gel. The value has to be
built from three separate groups of sources: (A) papers that specify the exact material, (B) the
closest measured gellan-gum-plus-salt dielectric data, and (C) primary data for the saline itself,
which dominates the result.

---

## A. Sources that specify the exact material (0.9 % saline + Phytagel)

These establish that the material is a recognised, published phantom formulation. **Neither paper
publishes a measured ε′ or ε″ for it** — they describe it qualitatively as "muscle equivalent".

**A1.** Qian PC, et al. *Irrigated Microwave Catheter Ablation Can Create Deep Ventricular Lesions
Through Epicardial Fat With Relative Sparing of Adjacent Coronary Arteries.*
Circulation: Arrhythmia and Electrophysiology. 2020. doi:10.1161/CIRCEP.119.008251
<https://www.ahajournals.org/doi/10.1161/CIRCEP.119.008251>

> "a 30:70 mixture of glycerol and 0.9% saline by volume hardened with 1% weight per volume gellan
> (Phytagel)"

**A2.** Qian PC, et al. *Transcatheter non-contact microwave ablation may enable circumferential
renal artery denervation while sparing the vessel intima and media.* EuroIntervention.
<https://eurointervention.pcronline.com/article/transcatheter-non-contact-microwave-ablation-may-enable-circumferential-renal-artery-denervation-while-sparing-the-vessel-intima-and-media>

> "a microwave muscle equivalent phantom layer consisting of 70% 0.9% saline and 30% glycerol mixed
> volume per volume, solidified with 1% weight per volume gellan (Phytagel; Sigma-Aldrich
> Corporation)"

> ⚠️ **VERIFY BEFORE QUOTING.** Both full texts are paywalled (HTTP 403 from EuroIntervention, AHA
> Journals and PubMed). The two quotations above were taken from the search index of those pages,
> not from a PDF that was read directly. Check them against the published PDF via the library before
> putting them in the report.

---

## B. Closest measured dielectric data for a gellan (Phytagel) gel with salt

**B1.** Zhang W, Luan D, Tang J, Sablani SS, Rasco B, Lin H, Liu F. *Dielectric properties and other
physical properties of low-acyl gellan gel as relevant to microwave assisted pasteurization
process.* Journal of Food Engineering. 2015;149:195–203.
<https://www.sciencedirect.com/science/article/abs/pii/S0260877414004348>

> "various amounts of sucrose (0, 0.1, 0.3 and 0.5 g/mL solution) and salt (0, 100, 200, and
> 300 mM) added to 1% gellan gel (with 6 mM Ca²⁺ addition)"

> "Results showed that sucrose and salt addition were effective for the adjustment of the dielectric
> constants and loss factor of gellan gels, respectively."

**Why this is the best match available:** 1 % gellan is the same chemistry and concentration as the
Phytagel in A1/A2, and 200 mM NaCl ≈ 1.17 % w/v — i.e. the 100–200 mM points bracket the
0.9–1.1 % range of interest. Reports ε′ and ε″ at 915 MHz over 22–100 °C, plus regression equations
for ε′ and ε″ vs. salt content, sucrose content and temperature. Paywalled; ScienceDirect and the
WSU-hosted PDF both refused direct access, so the numeric table has to be obtained through the
university library.

**B2.** *Effect of temperature on the dielectric properties of low acyl gellan gel.* Journal of Food
Engineering. 2012. <https://www.sciencedirect.com/science/article/abs/pii/S0260877412002087>

Companion study: 1 % gellan gel with 0.17 % and 0.3 % CaCl₂, measured 0.2–20 GHz over 10–100 °C.
Useful for the frequency dependence across the whole microwave band rather than just 915 MHz.

---

## C. Primary data for the saline, which sets the values

Phytagel at 1 % w/v is a >98 % water gel, so the complex permittivity is essentially that of the
saline solution. These are the standard citable sources.

**C1.** Peyman A, Gabriel C, Grant EH. *Complex permittivity of sodium chloride solutions at
microwave frequencies.* Bioelectromagnetics. 2007;28(4):264–274.
<https://www.researchgate.net/publication/6598547_Complex_permittivity_of_sodium_chloride_solutions_at_microwave_frequencies>

> "measured at 20 °C at concentrations between 0.001 and 5 mol/L and over a frequency range
> 0.13–20 GHz, and results were combined with literature values to derive empirical equations to
> predict the dielectric behavior of sodium chloride solutions between 0 and 5 mol/L and 5 °C–35 °C"

0.9 % NaCl = 0.154 mol/L and 1.1 % NaCl = 0.188 mol/L, both well inside this range. **This is the
paper to cite for the actual numbers.**

**C2.** Stogryn A. *Equations for calculating the dielectric constant of saline water.* IEEE Trans
Microwave Theory Tech. 1971;19(8):733–736. (Also Klein & Swift 1977.) Standard closed-form model for
the complex relative permittivity of saline water; convenient for putting ε*(f) directly into a
simulation.

**C3.** Gabriel S, Lau RW, Gabriel C. *The dielectric properties of biological tissues: II.
Measurements in the frequency range 10 Hz to 20 GHz.* Physics in Medicine and Biology.
1996;41(11):2251–2269. doi:10.1088/0031-9155/41/11/002
Reference tissue values to compare the phantom against.

---

## D. Supporting / context references

**D1.** Costanzo S, Cioffi V, Qureshi AM, Borgia A. *Gel-Like Human Mimicking Phantoms: Realization
Procedure, Dielectric Characterization and Experimental Validations on Microwave Wearable Body
Sensors.* Biosensors (Basel). 2021;11(4):111. doi:10.3390/bios11040111
<https://pmc.ncbi.nlm.nih.gov/articles/PMC8068187/>

Open access; verified directly. Saline-loaded gelatine (not Phytagel) phantoms, measured 500 MHz –
5 GHz. Measured values at 2.4 GHz (Table 3):

| Tissue | Measured ε′ | Measured σ (S/m) |
|---|---|---|
| Fat | 5.43 | 0.06 |
| Skin | 42.58 | 1.07 |
| Blood (1 % salt) | 57.01 | 1.59 |
| Muscle | 53.61 | 1.40 |

**D2.** Chen WJ, Wang Q, Kim CY. *Gel Phantom Models for Radiofrequency and Microwave Ablation of
the Liver.* Digestive Disease Interventions. 2020;4(3):303–310. doi:10.1055/s-0040-1716737
<https://pmc.ncbi.nlm.nih.gov/articles/PMC8297667/>

Verified directly. Confirms the gap in the literature:

> "Various other materials have been used to make gel phantoms, including TX-150, hydroxyethyl
> cellulose, gellan gum, carrageenan, and alginates, but their properties and use in thermal
> ablation studies are limited."

**D3.** Tejo-Otero A, et al. *Soft-Tissue-Mimicking Using Hydrogels for the Development of
Phantoms.* Gels. 2022;8(1):40. doi:10.3390/gels8010040
<https://pmc.ncbi.nlm.nih.gov/articles/PMC8774477/>

Verified directly. Material description only — mechanical properties, no dielectric data:

> "Phytagel is produced from a bacterial substrate that is composed of rhamnose, glucuronic acid,
> and glucose."

**D4.** McGarry CK, et al. *Tissue mimicking materials for imaging and therapy phantoms: a review.*
Physics in Medicine & Biology. 2020;65(23):23TR01. doi:10.1088/1361-6560/abbd17
<https://iopscience.iop.org/article/10.1088/1361-6560/abbd17>

General review of tissue-mimicking materials. Checked — contains no gellan/Phytagel dielectric data.

---

## E. Working values for the model

**Derived, not quoted.** Single-Debye water model (ε_s = 78.3, ε_∞ = 5.2, τ = 8.27 ps at 25 °C) with
the saline DC conductivity added as a loss term:

    ε′(f)  = ε_∞ + (ε_s − ε_∞)/(1 + (ωτ)²)
    ε″(f)  = (ε_s − ε_∞)·ωτ/(1 + (ωτ)²)  +  σ/(ωε₀)

Taking σ ≈ 1.5 S/m for 0.9 % NaCl and σ ≈ 1.8 S/m for 1.1 % NaCl at 25 °C:

| Composition | 915 MHz | 2.45 GHz |
|---|---|---|
| 0.9 % saline (σ ≈ 1.5 S/m) | ε* ≈ 78 − j33 | ε* ≈ 77 − j20 |
| 1.1 % saline (σ ≈ 1.8 S/m) | ε* ≈ 78 − j39 | ε* ≈ 77 − j23 |

Label these in the report as **estimated from saline data (C1/C2)**, not as measured values for the
gel. The 1 % Phytagel network is expected to lower ε′ by a few percent and slightly lower σ relative
to free saline, but no published measurement quantifies this for the 0.9–1.1 % saline case.

Note for the waveform derivation: the σ/(ωε₀) conduction term dominates ε″ below ~2 GHz, so at
915 MHz this material behaves as a conduction-dominated lossy medium; by 2.45 GHz the dipolar and
conduction loss mechanisms are comparable.
