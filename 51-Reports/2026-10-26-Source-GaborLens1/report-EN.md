# LhARA Stage 1 beam-line simulation: the laser-driven proton source and the first Gabor lens

| Item | Value |
|---|---|
| Date | 2026-10-26 |
| Code | `LhARAlinearOptics`, branch `claude/busy-euler-mjur5i` |
| Beam line | `11-Parameters/LhARABeamLine-Params-LsrDrvn-Gabor.csv` (TNSA laser-driven source + Gabor lenses; the Stage 1 main configuration, unmodified) |
| Events | 50 000 protons, fixed random seed 20261026 (fully reproducible) |
| Run time | Simulation 303 s; official envelope script `plotBeam.py` 44 s; this report's analysis about 4 min |
| Analysis script | `sourceGL1Analysis.py` in this directory (see Section 8 to reproduce) |

> The Chinese version is `report-CN.md`. Figure labels are in English so that both reports share the same figures.

---

## Summary

1. **Laser-driven proton source.** The code's default laser has 10 J in 40 fs (250 TW), a 1.5 µm focal-spot radius and an intensity of 2.1×10²¹ W/cm². The TNSA model then gives a hot-electron temperature Te = 10.5 MeV and a proton cut-off energy Kmax = 33.9 MeV. The generated spectrum falls roughly exponentially: median 5.4 MeV, mean 8.4 MeV, and only 3.55% of the protons lie between 14 and 16 MeV (around the 15 MeV reference energy). The source size is σ = 1.5 µm. The parameter `rpmax` restricts the emission directions to r′ < 5 mrad, and the RMS geometric emittance is 2.8×10⁻³ mm·mrad.
2. **Two implicit assumptions of the source model.**
   - `rpmax` keeps only about **0.08%** of the protons emitted into the TNSA cone, so every "transmission" in this report refers to that small-angle subset.
   - Energies are sampled from the angle-integrated spectrum, whereas higher-energy protons are more strongly concentrated inside a 5 mrad cone; the simulated spectrum is therefore softer than the true spectrum inside the cone.
3. **Source → entrance of the first Gabor lens (GL1)** (s = 0.25 m). The two nozzle apertures (radii 2.0 mm and 2.87 mm) are fully transparent to this 5 mrad beam, so the transmission is **100%**. The physical RMS beam size at the GL1 entrance is 0.63 mm.
4. **GL1**: length 0.857 m, electron density n_e = 5.1×10¹⁵ m⁻³, equivalent field B₀ = 1.40 T.
   - For 15 MeV protons: k = 1.55 m⁻², phase advance 1.07 rad, effective focal length 0.92 m.
   - GL1 reduces the RMS divergence of the whole beam from 1.88 to 1.03 mrad, but 15 MeV protons still diverge slightly after the lens (x′/x = +0.26 m⁻¹, i.e. a virtual source 3.9 m upstream of the exit). Only protons of about **11.3 MeV** leave GL1 parallel. This is consistent with the design, in which the capture is done by two Gabor lenses together.
5. **Strong chromatic focusing.** The focusing strength scales as k ∝ γ/p². Protons below 11.3 MeV are therefore over-focused (converging), and those below about 2.7 MeV cross the axis inside the lens.
   - The projected RMS emittance of the whole beam grows from 0.0028 to 1.35 mm·mrad (**about 480×**).
   - In ±0.25 MeV energy bands it grows only 4–23×, and for a single energy a linear lens conserves emittance exactly. The emittance growth is therefore almost entirely chromatic, caused by the broad spectrum.
6. **Electrons at the source** (computed with the model used by the code):
   - The hot-electron temperature is 10.5 MeV (the small-spot correction brings it below the 15.6 MeV of pure ponderomotive scaling), and there are about 1.8×10¹² hot electrons.
   - The rear-surface sheath has a radius of 3.8 µm, a density of 3.2×10²⁷ m⁻³ (about 1.9 times the critical density), a field of about 2.5×10¹³ V/m and a Debye length of 0.42 µm.
   - An independent estimate with Mora's isothermal-expansion model gives Kmax = 24–35 MeV, consistent with the code's 33.9 MeV.
7. **Can electrons be used to diagnose the proton source distribution? Partly.**
   - In the model, the shape of the proton spectrum is set by Te, and the cut-off energy by the laser-to-electron conversion efficiency together with the size of the electron sheath (electron divergence, target thickness). Measuring the electron temperature and sheath size therefore constrains the proton spectrum.
   - The constraint is model dependent and demands high precision: Kmax changes by about −1.3 MeV per degree of electron divergence, so predicting Kmax to ±1 MeV requires the divergence to ±0.8°.
   - Electrons cannot give the proton angular distribution or the details of the proton phase space.
   - The recommended use is as a shot-to-shot monitor, calibrated against direct proton measurements.
8. **Two code-level findings that affect the interpretation.**
   - For protons far from 15 MeV, the transverse positions given by the code are scaled by p/p₀ (Section 4.3).
   - 91% of the "losses" along the whole beam line come from the code's longitudinal cut (|z| > 5 m), not from the collimator (Section 5).

---

## 1. Simulation set-up

- **Beam line**: the main configuration file, with no parameter changed. The source is `SourceMode = 0` (parameterised TNSA). The file sets only `rpmax = 0.005`; all other laser parameters take the code defaults (Section 2.2).
- **Reference particle**: a 15 MeV proton, p₀ = 168.44 MeV/c, Bρ = 0.562 T·m. It was checked that the reference momentum is constant along the whole line (the RF cavities do not accelerate the reference particle at their reference phase).
- **Procedure**:
  1. `sourceGL1Analysis.py simulate` generates and tracks 50 000 protons with a fixed seed and writes the binary data file;
  2. the official scripts `plotBeam.py` and `visualiseBeam.py` process the same data file;
  3. `sourceGL1Analysis.py analyse` reads the data file and produces the tables and figures of this report.
- **Envelope inside elements**: the code records particles only at element exits, so the envelope inside each element is reconstructed in 20 steps using the element's own transfer matrix (for the Gabor lens, the code evaluates the matrix for each particle's energy). At each element exit the reconstruction is compared with the recorded particles: it agrees exactly in drifts, and to within 0.26% inside the lenses (from the interpolation in energy).
- **Statistical precision**: with 50 000 events, fractions carry an error of about ±0.08% (e.g. 3.55% ± 0.08%), and RMS sizes about 0.3%.

---

## 2. The laser-driven proton source

### 2.1 Model (TN Section 5)

The code uses the TNSA model of Fuchs et al. (Nature Physics 2, 48 (2006)):

- **Spectrum**: dN/dK ∝ exp(−√(2K/Te)) / √(2K·Te), for Kmin ≤ K ≤ Kmax.
- **Hot-electron temperature Te**: from ponderomotive scaling, with the small-focal-spot correction of Dover et al. (HEDP 37, 100847 (2020)).
- **Cut-off energy Kmax**: from t_laser/t₀ = X(1 + 1/(2(1−X²))) + ¼ ln((1+X)/(1−X)), with Kmax = X²·K∞. Here t₀ is the proton "ballistic time" across the sheath and K∞ is the limiting energy for infinitely long acceleration.
- **Emission direction**: the maximum emission angle is θ_max(K) = 20° − 15°·K/Kmax; inside this cone r′ is sampled from g(r′) ∝ r′ₘ² − r′². The condition r′ < `rpmax` is then applied.
- **Transverse position**: x and y are independent Gaussians with σ equal to the focal-spot radius r0.

### 2.2 Parameters (values used by the code)

| Parameter | Value | Origin |
|---|---|---|
| Laser wavelength λ | 0.8 µm | default |
| Laser power P | 2.5×10¹⁴ W (10 J / 40 fs) | default |
| Strehl ratio | 0.6 | default |
| Focal-spot radius r0 | 1.5 µm | default |
| Pulse duration τ | 40 fs | default |
| Target thickness d | 5 µm | default |
| Electron divergence half-angle θe | 25° | default |
| Intensity I = P·Strehl/(π r0²) | 2.1×10²¹ W/cm² | derived |
| Normalised vector potential a0 | 31.5 | derived |
| Hot-electron temperature Te | **10.50 MeV** | derived |
| Limiting energy K∞ | 94.8 MeV | derived |
| Ballistic time t₀ | 28.4 fs | derived |
| t_laser / t₀ | 1.41 | derived |
| X = √(Kmax/K∞) | 0.598 | derived |
| Kmin / Kmax | 1 MeV / **33.87 MeV** | default / derived |
| Angular parameters SigmaThetaS0, SlopeThetaS | 20°, 15° | default |
| rpmax | 0.005 (5 mrad) | beam-line file |

The full list is in `results/source_parameters.csv`.

### 2.3 Generated distributions

![Figure 1](figures/fig1_source_spectrum.png)

**Figure 1.** Left: the generated spectrum (histogram) agrees with the model curve. Right: the model's angle-integrated spectrum (the one sampled by the code) compared with the spectrum inside the cone r′ < 5 mrad.

![Figure 2](figures/fig2_source_angles.png)

**Figure 2.** Left and centre: emission directions fill a disc of radius 5 mrad uniformly; the distribution of r′ is ∝ r′, and σ(xₚ) = 2.51 mrad. Right: the model's emission cone (5°–20°) is far wider than the 5 mrad (0.29°) cut.

![Figure 3](figures/fig3_source_transverse.png)

**Figure 3.** Left: the transverse position at the source is Gaussian with σ = 1.50 µm. Centre and right: the trace-space slope x′ = pₓ/p₀ = xₚ·p/p₀ is proportional to the momentum, so higher-energy protons cover a wider range of x′.

Beam parameters at the source:

| Quantity | Value |
|---|---|
| Kinetic energy: mean / median / range | 8.43 MeV / 5.42 MeV / 1.00–33.87 MeV |
| Fraction with 14 < K < 16 MeV | 3.55% ± 0.08% |
| σx = σy | 1.50 µm |
| σx′ = σy′ (trace space, relative to p₀) | 1.88 mrad |
| εx = εy (RMS, geometric) | 2.82×10⁻³ mm·mrad |
| αx, βx | 0.002, 0.80 mm (x and x′ uncorrelated) |

### 2.4 Interpretation and caveats

- **Spectrum.** The TNSA spectrum is dominated by low-energy protons; only 3.55% lie within ±1 MeV of the 15 MeV reference energy. The beam line must therefore pick a narrow band out of a very broad spectrum, which is the origin of the chromatic effects and of the energy selection downstream.
- **What the `rpmax` cut means.** The model's emission half-angle is about 19.6° at 1 MeV, 13.4° at 15 MeV and 5° at Kmax. The cut r′ < 5 mrad keeps 0.04%–0.66% of these protons (about **0.08%** weighted by the spectrum). The "transmissions" in this report are therefore defined relative to this small-angle subset. To estimate the number of protons per laser shot delivered to the end station, multiply by this fraction and by the total TNSA yield.
- **Softer spectrum.** The code first samples the energy from the angle-integrated spectrum and then forces the direction into the 5 mrad cone, without making the acceptance depend on energy. Higher-energy protons are emitted into narrower cones, so they should make up a larger fraction inside 5 mrad (Figure 1, right). The code therefore underestimates the fraction of high-energy protons in the r′ < 5 mrad subset.
- **Source size.** The code takes σ = r0 = 1.5 µm (the laser spot), whereas the electron sheath on the rear surface has a radius of 3.8 µm in the model (Section 3). Experimentally, TNSA protons are usually emitted from a region larger than the laser spot, so the source size may be underestimated. At 25 cm, however, the source size is negligible compared with the beam size (0.0015 mm versus 0.6 mm).

---

## 3. Electrons at the laser source

> The code **does not track electrons**. Electrons affect the protons only through a few source-model parameters: the hot-electron temperature Te, the laser-to-electron conversion efficiency η, the electron divergence θe and the target thickness d. This section computes the electron parameters with the formulae implemented in the code (`Source.calculateTe`, `Source.calculatet0`) and the Fuchs model cited in TN Section 5, and checks them against Mora's model. The numbers are in `results/electron_source_model.csv`.

### 3.1 Physical picture

1. **Laser heating of electrons.** The intensity of 2.1×10²¹ W/cm² corresponds to a0 = 31.5 ≫ 1, deep in the relativistic regime. Electrons gain a transverse momentum of order a0·mₑc in the laser field (ponderomotive / J×B heating), forming a population of "hot" electrons with MeV temperatures.
2. **Small-spot effect.** The transverse excursion of an electron in the laser field, y0 = a0·λ/2π ≈ 4.0 µm, is larger than the effective spot radius r0/(2 ln 2) ≈ 1.1 µm. Electrons leave the focus before acquiring the full ponderomotive momentum and reach only 0.68 of it. Te therefore drops from the 15.6 MeV of pure ponderomotive scaling to **10.5 MeV** (Dover et al. 2020).
3. **Energy conversion.** The conversion-efficiency scaling η = 1.2×10⁻¹⁵·I^0.75 gives 11.9 at this intensity, which is larger than 1 and unphysical; the code caps it at 0.5. The TN notes that this scaling applies only for I < 3.1×10¹⁹ W/cm², while the intensity here is about 70 times higher. The number of hot electrons and K∞ therefore carry a large model uncertainty.
4. **Number of hot electrons and transport through the target.** N_e = η·E_spot/Te ≈ 1.8×10¹² (E_spot = Strehl × 10 J = 6 J). The electrons cross the 5 µm target with a half-angle of 25° and emerge from a region of radius B = r0 + d·tan θe = 3.8 µm on the rear surface, i.e. an area of 4.6×10⁻¹¹ m².
5. **Rear-surface sheath.** The hot-electron density is n_e0 = N_e/(c·τ·S) ≈ 3.2×10²⁷ m⁻³, about 1.9 times the critical density (1.74×10²⁷ m⁻³). Electrons crossing the rear surface set up a Debye sheath with a Debye length λ_D ≈ 0.42 µm and a field of order √(n_e0·Te/ε₀) ≈ 2.5×10¹³ V/m (25 TV/m). This field ionises the hydrogen-containing contamination layer on the rear surface and accelerates protons along the target normal: this is TNSA.
6. **Electron refluxing.** The sheath potential pulls most of the hot electrons back into the target, where they recirculate (reflux); only a small fraction of the most energetic electrons escapes into vacuum. This matters for electron diagnostics (Section 3.6).
7. **Time scales.** The laser pulse lasts 40 fs, the proton ballistic time is t₀ = 28.4 fs, and the product of the ion plasma frequency and the pulse length is ω_pi·τ ≈ 3.0. The protons gain most of their energy within about one laser-pulse duration.

### 3.2 Electron parameters from the model

| Quantity | Value |
|---|---|
| a0 | 31.5 |
| Transverse excursion y0 / effective spot radius | 4.0 µm / 1.08 µm |
| Small-spot correction p_t/p_t0 | 0.683 |
| Te (pure ponderomotive / corrected) | 15.6 MeV / **10.5 MeV** |
| Conversion efficiency η (raw scaling / used by code) | 11.9 / 0.5 |
| Number of hot electrons N_e | 1.8×10¹² |
| Sheath radius B, area | 3.83 µm, 4.6×10⁻¹¹ m² |
| Hot-electron density n_e0 | 3.2×10²⁷ m⁻³ (1.85 n_c) |
| Sheath-field scale √(n_e0 Te/ε₀) | 2.5×10¹³ V/m |
| Debye length λ_D | 0.42 µm |
| Ion-acoustic speed c_s | 0.106 c |
| ω_pi·τ_laser | 2.99 |
| Ballistic time t₀ | 28.4 fs |

### 3.3 How the electrons determine the proton source

In the model used by the code, the energy distribution of the source protons is fully determined by the laser and electron parameters:

| Electron quantity | Proton quantity affected | Relation |
|---|---|---|
| Hot-electron temperature Te | Shape (slope) of the spectrum | dN/dK ∝ exp(−√(2K/Te))/√K |
| Conversion efficiency η | Limiting energy K∞ | K∞ = 2mₑc²·√(η·P·Strehl/P_R) |
| Electron divergence θe, target thickness d | Sheath radius B → ballistic time t₀ → cut-off Kmax | t₀ = B/v∞, Kmax = X²·K∞ |
| Electron number N_e, density n_e0 | Total proton yield (normalisation of the spectrum) | dN/dK ∝ n_e0·c_s·t_laser·S |

The proton **angular distribution** (SigmaThetaS0, SlopeThetaS) and **transverse size** (σ = r0) are independent empirical inputs in the code, with no link to the electron parameters.

### 3.4 Cross-check with Mora's isothermal-expansion model

Mora's model of isothermal plasma expansion (PRL 90, 185002 (2003), in `92-Literature/`) gives an independent estimate of the cut-off energy:

Kmax = 2·Z·Te·[ln(τ_N + √(τ_N² + 1))]², where τ_N = ω_pi·t/√(2e) and e is Euler's number.

| Acceleration time t | Kmax (Mora) |
|---|---|
| t = τ_laser | 23.9 MeV |
| t = 1.3 τ_laser | 34.6 MeV |
| Code model (Fuchs) | 33.9 MeV |

The two independent models agree on the scale of the cut-off energy, so the code's source model is self-consistent for these parameters.

### 3.5 Sensitivity

Kmax as a function of electron divergence and target thickness is computed with the code's own functions. The effect of changing Te is obtained by reweighting: instead of re-simulating, each of the 50 000 simulated protons is weighted by the ratio of the new to the old spectrum.

![Figure 8](figures/fig8_electrons.png)

**Figure 8.** Left: cut-off energy versus electron divergence and target thickness. Centre: proton spectrum for different Te. Right: fraction of protons with 14 < K < 16 MeV (blue) and beam size at the GL1 exit (red) versus Te.

**Sensitivity of the cut-off energy to the electron divergence** (target thickness 5 µm; the 14–16 MeV fraction is the analytic model value):

| θe | 10° | 15° | 20° | **25°** | 30° | 35° | 40° | 45° |
|---|---|---|---|---|---|---|---|---|
| Kmax (MeV) | 57.0 | 48.6 | 40.8 | **33.9** | 27.9 | 22.7 | 18.3 | 14.5 |
| Fraction 14–16 MeV | 3.42% | 3.48% | 3.56% | 3.67% | 3.81% | 4.00% | 4.24% | 1.31% |
| σx at GL1 exit (mm), reweighted | – | – | – | 1.44 | 1.29 | 1.15 | 1.00 | 0.84 |

- Near the default value, dKmax/dθe = **−1.29 MeV per degree**.
- The target thickness matters just as much: at θe = 25°, d = 2, 5 and 10 µm give Kmax = 56.0, 33.9 and 16.0 MeV.
- A larger θe means a larger sheath, slower acceleration and a lower Kmax. By θe ≈ 45°, Kmax falls below 15 MeV and protons at the reference energy almost disappear.

**Sensitivity of the spectrum and the beam to the electron temperature** (Kmax fixed, reweighted):

| Te/Te0 | 0.7 | 0.8 | 0.9 | 1.0 | 1.1 | 1.2 | 1.3 |
|---|---|---|---|---|---|---|---|
| Te (MeV) | 7.35 | 8.40 | 9.45 | 10.50 | 11.55 | 12.60 | 13.65 |
| Fraction 14–16 MeV | 3.18% | 3.33% | 3.45% | 3.55% | 3.64% | 3.72% | 3.79% |
| σx at GL1 exit (mm) | 1.33 | 1.37 | 1.41 | 1.44 | 1.46 | 1.48 | 1.51 |
| Fraction reaching the end of the line | 2.44% | 2.55% | 2.65% | 2.73% | 2.80% | 2.86% | 2.91% |

A ±10% change in Te changes the 14–16 MeV fraction by only about ±3% (relative) and the beam size at the GL1 exit by about ±2%. With Kmax fixed, the number of protons delivered near 15 MeV is therefore not very sensitive to Te.

### 3.6 Can electrons be used to diagnose the proton source distribution?

**Conclusion: partly. Electrons constrain the proton energy distribution and the size of the emitting region indirectly, but they cannot replace direct proton diagnostics.**

**What electrons can constrain**
1. **Spectral shape.** In the model the spectral slope is set directly by Te. Measuring the hot-electron temperature (for example the spectrum of escaping electrons with a magnetic spectrometer, or the bremsstrahlung hard-X-ray spectrum of electrons in the target) constrains the shape of the proton spectrum.
2. **Cut-off energy.** Kmax is set jointly by the sheath size (θe, d) and the conversion efficiency. Measuring the electron divergence and the size of the emitting region on the rear surface constrains Kmax. Common techniques include Kα imaging with a buried fluorescence layer (electron transport and divergence inside the target) and imaging of optical transition radiation (OTR/CTR) from the rear surface (size and time structure of the region where electrons exit the target).
3. **Size of the emitting region.** The size of the rear-surface electron sheath (3.8 µm in the model) is a better proxy for the true TNSA emitting region than the proton source size used by the code (1.5 µm). Sheath imaging could therefore be used to improve the transverse size in the source model.
4. **Shot-to-shot monitoring.** Electron signals arrive before the protons and can be measured without intercepting the proton beam. For a high-repetition-rate facility such as LhARA, they could monitor the state of the source on every laser shot.

**What electrons cannot give**
1. **Proton angular distribution.** The angular parameters in the code are empirical and not linked to the electrons in the model. In reality the angular distribution depends on the spatial shape of the sheath; inferring it from electrons would need an additional model and would carry large errors.
2. **Details of the proton phase space**, such as the virtual-source position, laminarity and energy–angle correlations. These require direct proton measurements (radiochromic-film stacks, Thomson-parabola spectrometers, time of flight).

**Precision required and difficulties**
- Kmax is very sensitive to the electron parameters. Predicting it to ±1 MeV requires θe to ±0.8°, and the target thickness and conversion efficiency to comparable precision. These quantities are themselves hard to measure accurately, so a Kmax prediction from electrons can only be semi-quantitative.
- The constraint from Te is robust but weak: a ±10% error on Te corresponds to only about ±3% in the fraction of protons near 15 MeV.
- **Escaping electrons are not the sheath electrons.** The target charges up and only a small, most energetic fraction of the electrons escapes. The relation between the measured escaping-electron temperature and the sheath Te is empirical and must be cross-calibrated against proton measurements.
- The model itself is uncertain: the conversion-efficiency scaling is used far outside its range of validity (about 70 times the intensity limit).

**How to do it with this code.** The source section of the beam-line file accepts `Te`, `Kmax`, `DivAngle` (degrees) and `Thickness` (metres) directly. Parameters that are given are used as they are; only missing ones are computed from the laser parameters. For example:

```
1,Source,Source,Parameterised TNSA,Te,9.5,MeV,hot-electron temperature measured with an electron spectrometer
1,Source,Source,Parameterised TNSA,DivAngle,30,degrees,electron divergence half-angle from K-alpha imaging
1,Source,Source,Parameterised TNSA,Thickness,5.E-6,m,target thickness
```

A suggested study: build an "electron-constrained source" from measured electron parameters, simulate it, and compare with the directly measured proton spectrum. The size of the discrepancy is a quantitative measure of how far the electron diagnostic can be trusted.

---

## 4. From the source to the first Gabor lens

### 4.1 Beam-line elements (capture section)

| Element | s (exit, m) | Length | Notes |
|---|---|---|---|
| Source | 0 | – | TNSA source |
| Interface: Drift 1 | 0.05 | 0.05 m | |
| Interface: Aperture 1 | 0.05 | – | circular, radius 2.0 mm |
| Interface: Drift 2 | 0.10 | 0.05 m | |
| Interface: Aperture 2 | 0.10 | – | circular, radius 2.87 mm |
| Capture: Drift 1 | 0.25 | 0.15 m | to GL1 |
| **Capture: Gabor lens 1** | **1.107** | **0.857 m** | **GL1** |
| Capture: Drift 2, Drift 3 | 1.257, 1.407 | 0.15 m each | |
| Capture: Gabor lens 2 | 2.264 | 0.857 m | GL2 |
| Capture: Drift 4 | 2.414 | 0.15 m | end of capture section |

### 4.2 Parameters of the first Gabor lens

| Quantity | Value |
|---|---|
| Length L | 0.857 m |
| Strength in the beam-line file | 2.4917 (the unit column says rad/m) |
| Conversion inside the code | Strength × Bρ(15 MeV) = 1.40, i.e. an equivalent solenoid field B₀ = 1.40 T |
| Electron density n_e | 5.1×10¹⁵ m⁻³ (5.1×10⁹ cm⁻³) |
| Focusing strength k(15 MeV) | 1.55 m⁻² (√k = 1.246 m⁻¹, exactly Strength/2) |
| Phase advance √k·L (15 MeV) | 1.07 rad (61°) |
| Effective focal length f = −1/M₂₁ (15 MeV) | 0.92 m |
| Transfer matrix (identical in x and y, no coupling) | M = [[0.482, 0.703 m], [−1.091 m⁻¹, 0.482]] |
| Energy dependence of the focusing | k = e²·m·γ·n_e/(2ε₀p²) ∝ γ/p² |

### 4.3 Envelope

![Figure 4](figures/fig4_envelope_capture.png)

**Figure 4.** RMS envelope (σx and σx′) from the source to the end of the capture section. Solid black: all protons as given by the code; dashed black: all protons with positions rescaled by p₀/p; coloured lines: energy bands. Orange shading marks the Gabor lenses and dotted lines the apertures.

| Location | s (m) | Transmission | σx, code (mm) | σx, physical (mm)* | σx′ (mrad) | εx (mm·mrad) | αx | βx (m) |
|---|---|---|---|---|---|---|---|---|
| Source | 0 | 100% | 0.0015 | 0.0015 | 1.88 | 0.0028 | 0.00 | 0.0008 |
| Aperture 1 | 0.05 | 100% | 0.094 | 0.125 | 1.88 | 0.0028 | −62.4 | 3.12 |
| Aperture 2 | 0.10 | 100% | 0.188 | 0.251 | 1.88 | 0.0028 | −124.9 | 12.5 |
| **GL1 entrance** | 0.25 | 100% | 0.469 | **0.627** | 1.88 | 0.0028 | −312.2 | 78.1 |
| **GL1 exit** | 1.107 | 100% | 1.436 | **1.406** | 1.03 | **1.35** | −0.43 | 1.53 |
| 0.15 m after GL1 | 1.257 | 100% | 1.504 | 1.451 | 1.03 | 1.35 | −0.55 | 1.67 |
| GL2 exit | 2.264 | 100% | 2.017 | 2.130 | 0.70 | 1.20 | −0.63 | 3.38 |
| End of capture section | 2.414 | 64.0% | 2.018 | 2.090 | 0.74 | 1.32 | −0.52 | 3.09 |

\* "Physical" multiplies each proton's position by p₀/p; it is meaningful only for s ≥ 0.05 m (see below). The y results equal the x results within statistics. The full table is `results/beam_parameters_capture.csv`.

**Why positions are scaled by p/p₀.** The code defines the trace-space slope as x′ = pₓ/p₀, relative to the reference momentum. Drifts use a fixed matrix at the reference momentum (x → x + L·x′), while the Gabor lens uses the standard thick-lens matrix evaluated with the particle's own momentum. Combined, these mean that for a proton of momentum p the position given by the code is p/p₀ times the physical position, while the slope x′ itself is correct. Therefore:
- **Ratios are unaffected**, e.g. x′/x, focus positions and which energy leaves parallel (Figure 6).
- **Results near 15 MeV are accurate.**
- **Beam sizes for protons far from 15 MeV are biased**: they are too small by a factor 0.57 at 5 MeV and too large by a factor 1.42 at 30 MeV. At the GL1 entrance, the code's RMS size of the whole beam (0.47 mm) is 25% smaller than the physical one (0.63 mm). The physical value, 0.63 mm ≈ 0.25 m × 2.5 mrad, is exactly what a point source gives at 25 cm.
- This is a validity limit of linear optics for a ±100% momentum spread and should be stated in any thesis or paper.

### 4.4 Trace space

![Figure 5](figures/fig5_trace_space_GL1.png)

**Figure 5.** (x, x′) at the GL1 entrance, at the GL1 exit and 0.15 m after the exit; the colour shows the kinetic energy.

- **Entrance**: all protons lie on the line x′ = x/0.25 m, the signature of a point source after a drift. The emittance here is tiny.
- **Exit**: the distribution opens into a "bow tie".
  - High-energy protons (yellow-green) lie on a line of positive slope: they are under-focused and still diverging.
  - Medium- and low-energy protons (blue-purple) are rotated to negative slopes: they are converging.
  - Protons of 1–3 MeV (dark purple) are strongly over-focused and form a vertical band around x ≈ 0; protons of about 2.7 MeV are imaged onto the axis exactly at the exit.
- **0.15 m after the exit**: the shear continues and the high- and low-energy parts separate further.

### 4.5 Chromatic focusing

![Figure 6](figures/fig6_chromatic_GL1.png)

**Figure 6.** Left: x′/x at the GL1 exit versus energy. Points are simulated protons and the red line is the analytic result from the Gabor-lens transfer matrix; they coincide exactly, which validates the analysis. Right: effective focal length and phase advance of GL1 versus energy.

- x′/x > 0 means diverging, = 0 parallel and < 0 converging.
- **Protons of about 11.3 MeV leave parallel.** At 15 MeV, x′/x = +0.26 m⁻¹: still slightly diverging, with a virtual source about 3.9 m upstream of the exit.
- The lower the energy, the stronger the convergence: at 10 MeV, x′/x ≈ −0.17 m⁻¹ and the protons focus about 5.9 m after the exit; at 5 MeV, x′/x ≈ −2.2 m⁻¹ and they focus about 0.45 m after the exit.
- Protons of about 2.7 MeV are imaged onto the axis exactly at the exit (x = 0, so x′/x diverges there: the vertical line in the figure). Lower-energy protons cross the axis inside the lens, and protons of about 1.1 MeV once again leave parallel.
- The phase advance inside the lens exceeds π/2 below about 7 MeV and reaches π at about 1.8 MeV. At a phase advance of π, M₂₁ = 0 and the effective focal length diverges (right panel).
- The effective focal length is 0.92 m at 15 MeV and about 1.6 m at 30 MeV: it varies strongly with energy.

### 4.6 The emittance growth is chromatic

| Energy band | Protons | εx at GL1 entrance (mm·mrad) | εx at GL1 exit (mm·mrad) | Growth |
|---|---|---|---|---|
| All (1–34 MeV) | 50 000 | 0.0028 | 1.35 | ×480 |
| 5 ± 0.25 MeV | 1548 | 0.0021 | 0.048 | ×23 |
| 10 ± 0.25 MeV | 751 | 0.0031 | 0.040 | ×13 |
| 15 ± 0.25 MeV | 433 | 0.0038 | 0.030 | ×8 |
| 20 ± 0.25 MeV | 277 | 0.0040 | 0.025 | ×6 |
| 30 ± 0.25 MeV | 165 | 0.0056 | 0.021 | ×4 |

- The narrower the band and the higher its energy (i.e. the smaller its relative energy spread), the smaller the growth.
- For a single energy, the Gabor-lens transfer matrix has unit determinant and the emittance is conserved exactly.
- The 480-fold growth of the projected emittance is therefore almost entirely chromatic: protons of different energies are rotated into different directions in phase space.
- For TNSA beams, the projected emittance of the whole beam is thus not a good figure of merit. More meaningful are **the emittance within the selected energy band** and **the transmission of that band**.

### 4.7 Summary of this section

- Nothing is lost between the source and the GL1 entrance: the nozzle apertures are fully transparent to the 5 mrad beam.
- GL1 is a strong lens (focal length 0.92 m at 15 MeV), but on its own it does not make the 15 MeV beam parallel. Capture is done by GL1 and GL2 together: at the GL2 exit the RMS divergence of the whole beam is reduced further, to 0.70 mrad.
- The strong chromatic focusing of GL1 separates protons of different energies in position and angle. This is the physical basis of the energy selection downstream (collimator + RF cavities), and it is also why the projected emittance grows so much.

---

## 5. Overview of the whole beam line (brief)

![Figure 7](figures/fig7_full_line.png)

**Figure 7.** Left: transmission along the whole beam line (log scale). Right: energy spectra at the source and at the end of the line.

- The beam line is about 21 m long. **1363 protons reach the end, 2.73% ± 0.07% of those generated**, with energy **14.93 ± 0.44 MeV** (RMS): the line selects a very narrow band.
- Beam at the end (code values): σx = 2.09 mm, σy = 1.28 mm; εx = 0.58 mm·mrad, εy = 0.54 mm·mrad; βx = 7.5 m, βy = 3.0 m.
- **Classification of losses** (re-applying the checks of `BeamLineElement.Transport` in order):

| Reason | Protons | Fraction of generated | Energy range | Location |
|---|---|---|---|---|
| Longitudinal cut \|z\| > 5 m | 45 710 | **91.4%** | 1.0–33.9 MeV | s = 2.26–14.5 m |
| RPLCswitch at the entrance of the arc (cause not yet understood) | 1 727 | 3.5% | 13.1–16.9 MeV | s = 14.55 m |
| Energy-selection collimator (aperture) | 910 | 1.8% | 10.4–19.6 MeV | s = 5.85 m |
| Vacuum chamber (radius 35 mm) | 290 | 0.6% | 14.1–15.9 MeV | s = 16.0–17.2 m |

- **The longitudinal cut dominates.** The code requires |z| ≤ 5 m at the entrance of every element. In the linear approximation, z grows by δ/(β₀²γ₀²) per metre, so protons further from 15 MeV are removed earlier. The transmission and the spectrum at the end of the line are therefore determined mainly by this cut, with the collimator contributing only a small part. These are not losses on a physical aperture; **before drawing conclusions from the full-line transmission, it is worth confirming the physical meaning of this cut with the supervisor.**
- Another group of protons close to 15 MeV is removed at the RPLCswitch at the entrance of the arc. The cause has not yet been identified and is recorded in the to-do list.
- All of these occur at s > 2.2 m and do not affect the source and GL1 (s ≤ 1.26 m), which are the subject of this report.

Output of the official tools (first page):

![Figure 9](figures/fig9_official_beam_progression.png)

**Figure 9.** Beam evolution along the whole line from `plotBeam.py` (full PDF: `99-Scratch/LhARA-LsrDrvn-Gabor-BeamProgression.pdf`; data: `results/plotBeam_full_line.csv`).

![Figure 10](figures/fig10_official_visualise.png)

**Figure 10.** Trajectories of 100 protons from `visualiseBeam.py`.

---

## 6. Conclusions

1. With the code's default laser parameters (10 J, 40 fs, 1.5 µm spot), the TNSA source produces an exponential-like spectrum with Te = 10.5 MeV and Kmax = 33.9 MeV, of which only 3.55% lies within 15 ± 1 MeV.
2. The main configuration uses `rpmax` to restrict the source to 5 mrad, which amounts to simulating only about 0.08% of the protons emitted into the TNSA cone. This subset reaches GL1 through the nozzle without any loss.
3. GL1 (n_e = 5.1×10¹⁵ m⁻³) has a focal length of 0.92 m at 15 MeV and roughly halves the beam divergence, but the 15 MeV beam still diverges slightly after the lens; protons of about 11.3 MeV leave it parallel.
4. The chromatic focusing of GL1 (k ∝ γ/p²) is the dominant physical effect in this section: low-energy protons are over-focused and the projected emittance grows about 480 times, while the growth within narrow energy bands is much smaller. Chromaticity is both the basis of the energy selection and the main cause of the loss of beam quality.
5. The electrons at the source determine the proton spectrum: Te sets the slope and the sheath size sets the cut-off. Constraining the proton spectrum through electron diagnostics is **partly feasible**: electrons can serve as a shot-to-shot monitor and constrain the spectral parameters, but the approach is model dependent, demands high precision (θe to ±0.8°), cannot give the proton angular distribution or phase space, and must be cross-calibrated against direct proton diagnostics.

---

## 7. Limitations and suggested next steps

**Limitations**
- Linear optics without space charge; near the target the beam density is high and space charge is not negligible.
- The Gabor lens is treated as having an ideal, uniform electron density.
- With a ±100% momentum spread, the positions of protons far from 15 MeV are scaled by p/p₀ (Section 4.3).
- In the source model the conversion efficiency is used far outside its range of validity, the angular distribution is empirical, and the source size is taken as the laser spot.
- The laser parameters are the code defaults and are not necessarily the current LhARA design values.

**Suggestions**
1. Confirm the LhARA laser parameters with the supervisor, put them into your own copy of the beam-line file and repeat this analysis (commands in Section 8).
2. Remove the `rpmax` cut, or open it up to the nozzle acceptance (about 29 mrad), and study how the capture efficiency and emittance depend on the acceptance.
3. Study how the electron density and position of GL1 (currently 0.25 m from the source) affect the parallelism of the 15 MeV beam, e.g. under what conditions a single lens makes the 15 MeV beam parallel.
4. Clarify the physical meaning of the |z| > 5 m cut and of the losses at the RPLCswitch before evaluating the full-line transmission.
5. Following Section 3.6, build an "electron-constrained source" from measured electron parameters and compare it with proton measurements.

---

## 8. Files and reproduction

**Reproduction** (from the top directory of the repository, with the environment set up):

```bash
source startup.bash
python 51-Reports/2026-10-26-Source-GaborLens1/sourceGL1Analysis.py simulate -n 50000 -s 20261026 -d 99-Scratch/LhARA-LsrDrvn-Gabor.dat
python 51-Reports/2026-10-26-Source-GaborLens1/sourceGL1Analysis.py analyse -d 99-Scratch/LhARA-LsrDrvn-Gabor.dat
python 03-Scripts/plotBeam.py -i 99-Scratch/LhARA-LsrDrvn-Gabor.dat -n 50000 -o 99-Scratch/LhARA-LsrDrvn-Gabor-beam.csv
python 03-Scripts/visualiseBeam.py -i 99-Scratch/LhARA-LsrDrvn-Gabor.dat -o 99-Scratch/LhARA-LsrDrvn-Gabor-vis.pdf -n 100
```

The data file (84 MB) is kept in `99-Scratch/` and is not committed to git; the commands above regenerate it exactly.

**Result files** (`results/`)

| File | Content |
|---|---|
| `source_parameters.csv` | Source parameters (inputs and derived values) |
| `electron_source_model.csv` | Electron and sheath parameters, Mora cross-check |
| `electron_sensitivity_Te.csv`, `electron_sensitivity_angle.csv` | Sensitivity to the electron temperature and divergence |
| `beam_parameters_capture.csv` | Beam parameters along the capture section (`*_corr` columns: positions rescaled by p₀/p, meaningful only for s ≥ 0.05 m) |
| `GL1_chromatic.csv` | GL1 transfer matrix, focal length and x′/x versus energy |
| `transmission_full_line.csv` | Number of protons at each element of the whole line |
| `loss_classification.csv` | Classification of losses |
| `plotBeam_full_line.csv` | Output of the official `plotBeam.py` |
| `summary.json` | All key numbers |
