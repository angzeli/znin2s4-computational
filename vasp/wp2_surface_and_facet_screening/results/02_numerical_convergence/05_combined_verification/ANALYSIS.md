# WP2 Stage 02 — Final Combined Verification

## Executive conclusion

**VERDICT: PASS WITH DOCUMENTED LIMITATION.** Stage 02 can close for the specified clean β(001) screening model: **four complete ZnIn₂S₄ layers, full ionic relaxation, 500 eV, Γ-centred 10×10×1 sampling and nominal 20 Å vacuum**. The existing accepted **J34.1** already contains this combined recipe. A duplicate slab calculation in `05_combined_verification` is not required.

The matched bulk static **J39.1 passes the fixed-geometry electronic-convergence checks**. Its E₀ reference is **−28.75482061 eV per formula unit**. With this reference, the 2L→4L paired surface-excess change is **+1.871971 meV/Å² (FAIL)**, whereas **4L→6L is +0.246638 meV/Å² (PASS)** against the predefined 1 meV/Å² criterion. F and energy-without-entropy comparisons give the same verdict. Both face-specific work functions and all four vacuum-aligned frontier extrema pass the predefined 0.05 eV screening target for 4L→6L. The relaxed-region comparison also passes for those measured observables; FULL remains the preferred treatment.

**No unresolved Stage 02 acceptance criterion requires an 8L slab.** An 8L calculation could add information about the interior profile or the longer thickness trend, but that information is not necessary to choose the present screening model. The documented limitation is that the **full interior electrostatic profile is still thickness-dependent**, despite agreement of the measured surface descriptors. This closure does not establish a bulk-like interior, rigorous bulk band edges, universal surface-energy accuracy or convergence of downstream adsorption/interface properties.

Updated on 22 September 2026 against the archived J39.1 rerun. Acceptance uses the predefined [ROADMAP §5.6.6](../../../ROADMAP.md#566-convergence-acceptance) thresholds. This report integrates the native endpoint evidence with the [thickness](../03_slab_thickness/THICKNESS_ANALYSIS.md) and [relaxed-region](../04_relaxed_region/ANALYSIS.md) analyses.

## Matched β bulk reference

### Provenance and geometry

Source: [bulk_reference](../../../calculation/02_numerical_convergence/05_combined_verification/bulk_reference), **J39.1**, attempt `8447cc08b9f04f7482fb3081d943cee8`. CMW records normal completion with exit 0; archived OUTCAR, OSZICAR and runtime metadata match this registered attempt. All three endpoint energy definitions reproduce the earlier J38.1 result to the printed 10⁻⁸ eV precision. The scientific checks below use the archived J39.1 outputs.

The POSCAR parses as **S₄In₂Zn: seven atoms, one ZnIn₂S₄ formula unit**. Its cell and coordinates match the accepted [WP1 β CONTCAR](../../../../wp1_polymorph_polytype_benchmark/calculation/01_geometry_optimisation/beta/CONTCAR). This is the accepted periodic β parent used to construct the slabs, with no Selective Dynamics constraints or added vacuum. The 2.728575 Å interlayer gap belongs to the bulk structure.

Cell vectors, in Å, are `(3.8782842923155867, 0, 0)`, `(−1.9391421461577933, 3.3586927202930394, 0)` and `(0, 0, 12.1484331332749314)`. Thus a≈b=3.8782842923 Å, c=12.1484331333 Å, α=β=90°, γ=120° and **V=158.2450674694 Å³**. Coordinates, lattice and volume are finite. There is exactly one fixed-geometry evaluation; final CONTCAR and XML agree to approximately 6.1×10⁻⁸ Å, and input/final geometry differences are only output-rounding effects. No ionic move or geometry optimization occurred.

### Executed model and static convergence

| Setting | Native-output value |
| --- | --- |
| VASP / PAW | VASP 6.6.1; PAW_PBE S, In_d, Zn, each dated 06Sep2000 |
| Functional / dispersion | PBE (`GGA=PE`); D3(BJ), `IVDW=12` |
| Basis / projections | ENCUT=500 eV; PREC=Accurate; LREAL=False; LASPH=True; ADDGRID=False |
| Spin / occupations | ISPIN=1; ISMEAR=0; SIGMA=0.05 eV; NELECT=62 |
| SCF algorithm / mixing | ALGO=Normal (IALGO=38); AMIN=0.01; AMIX=0.40; BMIX=1.00 |
| SCF budget / threshold | EDIFF=10⁻⁶ eV; NELM=120; NELMIN=2 |
| Initialization | ISTART=0; ICHARG=2; fresh atomic-density start; no inherited wavefunction |
| Sampling / symmetry | Γ-centred 10×10×4; 50 irreducible k-points; ISYM=2 |
| Static / electrostatics | IBRION=−1; NSW=0; LDIPOL=False; IDIPOL=0; no active slab DIPOL |
| Parallel execution | 8 MPI ranks × 1 thread; NCORE=4; KPAR=1 |

OUTCAR independently confirms these settings; XML is used where tags are available. The S/In_d/Zn POTCAR triplet is byte-identical to each corresponding repeat in the accepted slabs. The bulk and slabs share the relevant physical/electronic settings. Appropriate differences are three-dimensional 10×10×4 rather than 10×10×1 sampling, bulk ISYM=2 rather than slab ISYM=0, no bulk dipole correction, and static rather than relaxing ions. Cell-dependent basis dimensions and band counts also differ. J36.1 alone has a larger NELM ceiling of 300; that changes the completion budget, not the model or convergence threshold.

The archived rerun uses `SYSTEM = WP2-P02-beta-bulk-matched-static-500-k10x10x4`, consistent with the **10×10×4** mesh in KPOINTS and OUTCAR. Relative to J38.1, the source INCAR changes only the descriptive label, comments and whitespace; the scientific settings are unchanged.

The bulk SCF genuinely reaches EDIFF in **15 electronic iterations**, well below NELM=120. Final **dE=−3.5740×10⁻⁷ eV** and **d eps=−6.2225×10⁻⁷ eV** both meet the criterion, with an explicit OUTCAR termination marker. This is the two-condition electronic stopping test described by [VASP](https://vasp.at/wiki/EDIFF), not an inference from exit 0. Final rms is **0.00188**; the last available rms(c) is **0.00216 at iteration 14** (the terminal row omits it). Late charge residuals decrease 0.00846→0.00518→0.00216 over steps 12–14; late energy changes also decrease. Early oscillation settles into a damped approach rather than recurrent late amplification.

All electronic diagnostics are finite. There is no BRMIX, non-Hermitian, NaN/Inf, EDDDAV, ZHEGV/ZPOTRF or serious internal/fatal-error evidence; stderr is empty. Complete XML and a normal OUTCAR timing footer are present. VASP elapsed time is **44.020 s**, with **44.578 s** for the managed attempt. Seven atoms, 38 bands, bulk symmetry and one SCF evaluation explain why this run is much faster than a slab relaxation.

The final maximum force is **0.00779471 eV/Å**, with RMS force 0.00456375 eV/Å. Stress diagonal components are **−0.053980, −0.053980, +0.361216 kbar**. These small residuals are consistent with the accepted geometry. The calculation evaluates that geometry without reoptimizing it. Weighted EIGENVAL occupations reproduce 62 electrons; EIGENVAL/XML eigenvalues and occupations agree within their printing precision.

### Energy definitions and normalization

OUTCAR and XML agree for all three endpoint energies. Because this bulk cell contains exactly one formula unit, **e_bulk equals the cell energy**, without another divisor.

| Definition | E_bulk (eV/cell) | e_bulk (eV/formula unit) |
| --- | --- | --- |
| F / TOTEN | -28.75482063 | -28.75482063 |
| Energy without entropy | -28.75482059 | -28.75482059 |
| E₀ = energy(sigma→0), primary | -28.75482061 | -28.75482061 |

Use the final ion–electron total energy, including dispersion. The last DAV free energy is −26.63857234 eV; its difference from final F is −2.11624829 eV, consistent with the printed D3 term (−2.11625 eV). Substituting the DAV value would therefore omit dispersion and invalidate the surface subtraction.

E₀ is the primary comparison, consistently for bulk and slabs. Its notation does not establish convergence with respect to smearing width: these are all Gaussian σ=0.05 eV calculations, and F/E_without_entropy remain useful diagnostics ([VASP smearing guidance](https://vasp.at/wiki/Smearing_technique)). Printed digits permit reproducible subtraction; they do not imply eight-decimal physical accuracy.

## Accepted slab endpoints

| Slab / job | Accepted native directory | Formula units / atoms | NELM | Converged SCFs / evaluations | Final max force (eV/Å) |
| --- | --- | --- | --- | --- | --- |
| 2L / J33.1 | [2L output](../../../calculation/02_numerical_convergence/03_slab_thickness/01_relaxation/beta001_2L) | 2 / 14 | 120 | 33 / 33 | 0.01355518 |
| 4L / J34.1 | [4L output](../../../calculation/02_numerical_convergence/03_slab_thickness/01_relaxation/beta001_4L) | 4 / 28 | 120 | 14 / 14 | 0.01802549 |
| 6L / J36.1 | [6L output](../../../calculation/02_numerical_convergence/03_slab_thickness/01_relaxation/beta001_6L/retry_nelm300_r1) | 6 / 42 | 300 | 16 / 16 | 0.01596172 |

**J35.1 is excluded.** The files at the 6L parent level preserve the failed attempt. All 6L scientific endpoints here come exclusively from `beta001_6L/retry_nelm300_r1`, confirmed against the J36.1 CMW identity and registered native outputs. The retry archive contains outputs rather than another complete input set; its original POSCAR/KPOINTS/POTCAR were checked through the source location recorded in its own RUN_METADATA. Its starting POSCAR is byte-identical to the original 6L geometry, and executed OUTCAR confirms NELM=300. No failed J35 energy enters a table or fit.

Formula-unit counts come from actual compositions Zn₂In₄S₈, Zn₄In₈S₁₆ and Zn₆In₁₂S₂₄. The unrounded one-face area is common to all slabs:

**A = |a × b| = 13.025965219827203 Å².**

Every accepted force evaluation has a genuine EDIFF marker and both final electronic energy changes below 10⁻⁶ eV. Electronic counts match XML; no accepted cycle exhausts NELM. All three contain explicit ionic convergence and normal timing footers, with finite trajectories and no serious numerical warning identified. Their final CONTCARs agree with force-evaluated XML structures within 4.5×10⁻⁷ Å. OUTCAR/XML endpoint energies and all force blocks agree within printed precision. Thus NSW or walltime exhaustion did not supply these endpoints.

| Energy (eV/slab cell) | 2L / J33.1 | 4L / J34.1 | 6L / J36.1 |
| --- | --- | --- | --- |
| F / TOTEN | -56.96701880 | -114.45163599 | -171.95807850 |
| Energy without entropy | -56.96444866 | -114.45034546 | -171.95675999 |
| E₀, primary | -56.96573373 | -114.45099072 | -171.95741924 |

| SCF / runtime diagnostic | 2L | 4L | 6L |
| --- | --- | --- | --- |
| First SCF iterations | 31 | 75 | 215 |
| Subsequent iterations: range; median | 4–20; 7 | 4–19; 6 | 4–16; 7 |
| VASP elapsed time (s) | 4630.193 | 12568.158 | 57253.738 |

## Paired surface-excess energies

For each energy definition separately:

`Γ_pair(N) = [E_slab(N) − N e_bulk] / A`.

Here N counts ZnIn₂S₄ formula units and A is the **one-face area**. This is the excess associated with the **pair of inequivalent β(001) faces**. It is not a single-face surface energy, and Γ_pair/2 is not assigned to either face. All arithmetic uses unrounded native energies and the cell-derived area; an independent decimal calculation checks the normalization and differences.

### Primary E₀ result

| Slab | E₀,slab − N e₀,bulk (eV) | Γ_pair (eV/Å²) | Γ_pair (meV/Å²) |
| --- | --- | --- | --- |
| 2L | 0.54390749 | 0.041755638 | 41.755638 |
| 4L | 0.56829172 | 0.043627609 | 43.627609 |
| 6L | 0.57150442 | 0.043874247 | 43.874247 |

| Refinement | ΔΓ_pair (meV/Å²) | Predefined criterion | Verdict |
| --- | --- | --- | --- |
| 2L→4L | +1.871971 | \|ΔΓ_pair\| ≤ 1 meV/Å² | FAIL |
| 4L→6L | +0.246638 | \|ΔΓ_pair\| ≤ 1 meV/Å² | PASS |

The **2L model fails** the energetic criterion by 0.871971 meV/Å². The **4L production candidate passes** its refinement to 6L using the same threshold and independent bulk reference.

### Energy-definition robustness and numerical margin

| Definition | Γ_pair(2L) | Γ_pair(4L) | Γ_pair(6L) | ΔΓ 2L→4L | ΔΓ 4L→6L |
| --- | --- | --- | --- | --- | --- |
| F / TOTEN | 41.656987 | 43.578078 | 43.823645 | +1.921091 | +0.245567 |
| Without entropy | 41.854290 | 43.677139 | 43.924849 | +1.822850 | +0.247709 |
| E₀ | 41.755638 | 43.627609 | 43.874247 | +1.871971 | +0.246638 |

All entries in this diagnostic table are **meV/Å²**. Every definition gives FAIL for 2L→4L and PASS for 4L→6L. The full F-to-without-entropy spread of the 4L→6L change is only **0.002142 meV/Å²**. E₀ is not close to the acceptance boundary: its margin is **0.753362 meV/Å²**, equivalent to **9.813265 meV** in the 6L-minus-4L excess-energy numerator.

For a bulk-reference shift δe_bulk at fixed slab endpoints, `δ(ΔΓ_4→6) = −2 δe_bulk/A`. A 1 meV/formula-unit shift changes this refinement metric by 0.15354 meV/Å²; the nearest boundary would require the bulk reference to become approximately **4.906633 meV/formula unit more negative**, with other quantities unchanged. This is a sensitivity calculation, not an estimated error bar. The historical WP1 OUTCAR E₀ is only 0.189120 meV/formula unit away, but that run also differs in mesh and ADDGRID, so it is supporting consistency rather than a controlled bulk mesh-convergence proof.

Finite ionic tolerance (0.02 eV/Å), SCF tolerance and sampling still limit accuracy. A force threshold alone does not bound the remaining relaxation energy; no independent tighter restart of these production endpoints is claimed. Nevertheless, the observed refinement is well inside the fixed project threshold and robust to the tested energy definitions. The new bulk run validates a reference **at 500 eV and 10×10×4**; it does not independently establish the complete-basis or dense-k limit of absolute Γ_pair. The reported pass applies to thickness refinement at this specified numerical model.

## Relation to slab-only energy scaling

The independent bulk reference remains authoritative. The previous slab-only diagnostics recompute as follows using E₀:

| Diagnostic | Value (eV/formula unit) | Difference from independent e_bulk (meV/formula unit) |
| --- | --- | --- |
| Independent matched bulk | -28.75482061 | 0 |
| ε₂→₄ = (E₄−E₂)/2 | -28.742628495 | +12.192115 |
| ε₄→₆ = (E₆−E₄)/2 | -28.753214260 | +1.606350 |
| Three-point least-squares slope | -28.7479213775 | +6.899233 |

The successive increment moves toward the independent bulk value: its mismatch falls from **12.192115 to 1.606350 meV/formula unit**. The curvature is `E₆−2E₄+E₂ = −21.171530 meV/cell`, or **−1.625333 meV/Å²** after dividing by A. It is the difference between successive Γ_pair increments, not either increment itself.

The fitted slope was qualitatively useful in suggesting improving extensivity, but was **6.899233 meV/formula unit less negative** than the independent bulk result. Using the fitted slope as the reference would give adjacent changes of approximately ±0.812666 meV/Å², apparently passing both steps. The equal-and-opposite changes arise from fitting the same three equally spaced thicknesses; their agreement with the threshold is not an independent convergence test. The independent reference materially clarifies the conclusion: **2L fails and 4L passes**. Agreement of the later increment with e_bulk supports this outcome but is not a separate proof of the asymptotic limit.

## Vacuum-referenced electronic convergence

Native LOCPOT/CHGCAR profiles were recomputed for the three accepted endpoints; their headers match the final geometries and their grids agree. LOCPOT contains the LVHAR Hartree-plus-ionic potential. Each face is referenced separately as `Φ_side = V_vac,side − E_F`, following [VASP's work-function method](https://vasp.at/wiki/Computing_the_Workfunction). A 6 Å setback from the outermost atom and 1 Å exclusion from the cell boundary select depleted, flat vacuum windows; no mean over the two faces is used.

All six windows pass the established extraction checks: width ≥2 Å, |slope| ≤0.005 eV/Å, detrended residual ≤0.02 eV, range ≤0.03 eV and maximum |density| ≤10⁻⁵ e/Å³. Actual widths are ≥2.79287 Å; maximum |slope| is 0.00039349 eV/Å, maximum range 0.00118307 eV and maximum |density| 1.07344×10⁻⁶ e/Å³. Integrated densities reproduce the required 124/248/372 electrons. The prior 5-versus-6 Å window sensitivity is below 0.000254 eV. These checks support the reference extraction; they are not a proof of every periodic-image energy error.

| Quantity (eV) | 2L | 4L | 6L | 2L→4L | 4L→6L |
| --- | --- | --- | --- | --- | --- |
| Lower-face Φ | 3.983900 | 3.767710 | 3.767954 | -0.216190 | +0.000244 |
| Upper-face Φ | 6.833582 | 6.930883 | 6.948772 | +0.097301 | +0.017889 |

**2L→4L fails 0.05 eV on both faces; 4L→6L passes on both.** The final face asymmetry persists: the upper-minus-lower vacuum steps are 2.849682, 3.163173 and 3.180818 eV.

For the established frontier diagnostic, `b=NELECT/2`: b=62/124/186 for 2L/4L/6L. Use `max_k ε_b` and `min_k ε_(b+1)`, retaining partially occupied states. EIGENVAL/XML arrays agree within printed precision and weighted occupations reproduce NELECT.

| Aligned extremum (eV) | 2L | 4L | 6L | 2L→4L | 4L→6L |
| --- | --- | --- | --- | --- | --- |
| Band-b maximum, lower vacuum | -3.821925 | -3.758227 | -3.782734 | +0.063698 | -0.024508 |
| Band-(b+1) minimum, lower vacuum | -3.968151 | -3.730159 | -3.738394 | +0.237992 | -0.008236 |
| Band-b maximum, upper vacuum | -6.671607 | -6.921399 | -6.963553 | -0.249792 | -0.042153 |
| Band-(b+1) minimum, upper vacuum | -6.817833 | -6.893331 | -6.919213 | -0.075498 | -0.025881 |

The largest absolute 4L→6L change is **0.042153 eV**, below 0.05 eV (margin 0.007847 eV); 2L→4L fails. The sampled frontier separations are −0.146226, +0.028068 and +0.044340 eV. Partial occupations persist; these values are not validated semiconductor gaps, rigorous bulk VBM/CBM values or evidence of identified surface-state localization. The pass is for these sampled, vacuum-aligned screening descriptors at the stated smearing.

## Structural and interior-electrostatic convergence

All accepted trajectories remain finite, preserve layer order and show no collapse or overlap. Minimum distances encountered are 2.28691/2.29448/2.29552 Å for 2L/4L/6L; largest ionic movements are 0.023573/0.017307/0.017222 Å. First-shell cation–S graphs retain their initial identities after periodic-image registration: 28/56/84 links at a 3 Å cutoff. The common observed surface topology and termination identities persist within the present lateral cell. Free energies decrease across evaluated geometries, while force norms show finite rebounds; monotonic force reduction is not required.

| Slab | Relaxed interlayer gaps, lower→upper (Å) | Actual final vacuum (Å) |
| --- | --- | --- |
| 2L | 2.841470 | 19.833551 |
| 4L | 2.758836, 2.757570, 2.752280 | 19.939621 |
| 6L | 2.756110, 2.749106, 2.745266, 2.743679, 2.748316 | 19.925823 |

The 2L gap is distinctly wider. From 4L to 6L the outer gaps change by only −0.002727 and −0.003965 Å; the central gap changes by −0.012304 Å. This supports improving structural consistency, without assigning an untested universal bond/displacement convergence tolerance. Starting vacuum is 20 Å; full relaxation slightly changes the actual atom-free separation while the cell remains fixed.

The interior analysis was recomputed from native grids using the existing convention: periodic piecewise-linear planar potentials, a **12.1484331333 Å box average**, registration at the central interlayer-gap midpoint, and the common local window `u ∈ [−L/4,+L/4]`. The averaging footprint remains inside every slab. Direct native-knot quadrature checks the integrals; no slope removal or fitted potential shift is applied.

| Interior descriptor | 2L | 4L | 6L |
| --- | --- | --- | --- |
| Central mean, lower vacuum (eV) | -8.738397 | -8.624581 | -8.574700 |
| Central mean, upper vacuum (eV) | -11.588079 | -11.787754 | -11.755518 |
| Residual slope (eV/Å) | 0.063051 | 0.029048 | 0.013620 |
| Central potential range (eV) | 0.345784 | 0.178816 | 0.083573 |
| Central-gap interstitial density (e/Å³) | 0.098655 | 0.106420 | 0.107312 |

The 4L→6L central means change by +0.049881/+0.032236 eV under lower/upper-vacuum alignment. Central-gap interstitial potentials change by +0.042862/+0.025216 eV, and density changes by **+0.838%**. These local measures support partial interior similarity. The full central profiles retain maximum differences of **0.097462/0.079816 eV**, with RMS differences **0.056747/0.042086 eV**. Slope and range decrease by about half. In 6L, the three inner gap-centre repeat averages span **0.346920 eV over 24.297451 Å**. Thus a broad thickness-invariant interior profile has not been demonstrated.

Interior convergence concerns **thickness invariance**; it does not require the field in an asymmetric slab to vanish. The measured profile changes limit interpretation of the interior field and bulk-like band edges, while the paired excess, work functions and aligned frontier descriptors meet their own refinement criteria. No separate slope threshold was predefined.

## Relaxed-region convergence

J37.1 fixes the central two complete layers (14 atoms) while relaxing the outer two; J34.1 relaxes all 28 atoms. The starting geometry, lattice, k mesh and PAW datasets match, and executed numerical settings agree apart from the descriptive SYSTEM label. The constraint flags are the intended intervention. The key quantities below were independently regenerated from the native endpoints rather than copied only from the earlier report.

| FIXED − FULL diagnostic | Verified result | Interpretation |
| --- | --- | --- |
| ΔE₀ / ΔΓ_pair | +0.00152675 eV/cell / +0.117208 meV/Å² | PASS: bulk term cancels at fixed composition and area |
| ΔΓ from F / without entropy | +0.110817 / +0.123600 meV/Å² | Same energetic PASS |
| Lower / upper ΔΦ | +0.000953 / −0.006344 eV | PASS: both <0.05 eV |
| Largest aligned frontier change | 0.032970 eV | PASS: <0.05 eV |
| Largest surface bond-length change | 0.002241 Å | Small change; no first-shell topology change |
| Largest outer-atom position difference | 0.052379 Å | Mostly outward layer translations, not identical positions |
| FIXED central-atom max / RMS force | 0.100972 / 0.068259 eV/Å | Real suppressed relaxation; not force-free |
| Largest aligned central-profile difference | 0.018612 eV | Similar profile with measurable perturbation |
| Central slope / range changes | +0.002817 eV/Å / +0.014841 eV | FIXED is slightly steeper, not more bulk-like |
| Runtime FULL / FIXED | 12568.158 / 17874.445 s | FIXED took 42.2% longer |

J37.1 has **28/28 genuinely converged electronic cycles**, explicit constrained ionic convergence and maximum free-atom force **0.01275239 eV/Å**. The nonzero fixed-atom forces are present in native OUTCAR/XML; they are not zeroed artifacts and are not subject to the free-atom stopping test. Fixed atoms remain fixed to output precision. The 56 cation–S links are unchanged at 2.8, 3.0 and 3.2 Å cutoffs, with no gained/lost surface neighbor or altered termination. No claim of an absent larger-cell reconstruction follows.

**Relaxed-region verdict: PASS for the intended measured screening observables.** This is not force equivalence, zero strain or proof that FIXED is more bulk-like. FULL allows the measured interior relaxation and had no runtime disadvantage, so keep **full ionic relaxation** for production.

## Selected production model

| Item | Selected setting verified in J34.1 |
| --- | --- |
| Surface / thickness | β(001), existing β001_t07-derived inequivalent face pair; four complete layers, Zn₄In₈S₁₆ |
| Lateral cell / area | Accepted β bulk in-plane cell; A=13.025965219827203 Å² |
| Vacuum | Nominal initial 20 Å; final atom-free gap 19.939621 Å; c=65.8651576100 Å |
| Relaxation | All atoms free; IBRION=2; ISIF=2; POTIM=0.5; NSW=100 |
| Force criterion | EDIFFG=−0.02 eV/Å; final maximum norm 0.01802549 eV/Å |
| Functional / PAWs | PBE+D3(BJ), IVDW=12; matching PAW-PBE S/In_d/Zn repeats |
| Cutoff / precision | ENCUT=500 eV; PREC=Accurate; LREAL=False; LASPH=True; ADDGRID=False |
| Sampling / symmetry | Γ-centred 10×10×1; ISYM=0 |
| SCF | ALGO=Normal; AMIN=0.01; AMIX=0.40; BMIX=1.00; EDIFF=10⁻⁶ eV |
| Electronic budget / start | NELM=120; NELMIN=2; fresh ISTART=0 / ICHARG=2 initially |
| Smearing / spin | ISMEAR=0; SIGMA=0.05 eV; ISPIN=1 non-spin-polarized basal reference |
| Electrostatics | LDIPOL=True; IDIPOL=3; DIPOL=(0.5,0.5,0.5); LVHAR=True |
| Resources / budget | 8 MPI ranks × 1 thread; NCORE=4; KPAR=1; existing 48-hour managed budget |

**J34.1 already embodies this exact combined production recipe.** Its numerical and ionic convergence are verified above. The independent J39.1 bulk energy now supplies the missing reference; J36.1 supplies thickness refinement and J37.1 the relaxed-region control. These existing calculations jointly constitute the combined verification; another calculation merely to populate a stage directory would add no changed scientific variable.

The earlier [k-point](../01_kpoints_and_encut/KPOINTS_ANALYSIS.md), [cutoff](../01_kpoints_and_encut/ENCUT_ANALYSIS.md) and [vacuum](../02_vacuum/VACUUM_ANALYSIS.md) studies selected 10×10×1, 500 eV and approximately 20 Å on the fixed basal seed. Their final native energies were rechecked here: k10→k12 changes E₀/A by approximately +0.0980 meV/Å²; 500→600 eV by −0.048947 meV/Å²; 20→25→30 Å successive changes are +0.069497 and +0.040477 meV/Å². Reported maximum work-function changes are 0.024501, 0.001838 and 0.002676 eV, respectively. All support the selected screening settings.

Those tests established sampling, cutoff and vacuum sensitivity on a single-layer geometry. They did not jointly refine those settings on the final 4L geometry or include a matched bulk calculation at each k-point/cutoff refinement. Consequently, they support the selected screening settings without establishing the full systematic error in absolute Γ_pair. The final endpoints verify the combined recipe, vacuum-reference extraction and the specified thickness/constraint comparisons.

## Need for additional calculations

- **Duplicate combined-verification VASP run: NO.** J34.1 contains the selected model and J39.1 supplies its independent bulk reference. There is no identified recipe mismatch that a duplicate would fix.
- **8L slab: NO for the current Stage 02 screening decision.** The 4L→6L paired excess changes by only 0.246638 meV/Å²; both work functions and all four aligned extrema meet 0.05 eV; qualitative structure is stable and the relaxed-region control passes. **No unresolved Stage 02 acceptance criterion requires an 8L slab.**

This does not claim that 8L contains no additional information or that the infinite-thickness limit has been proved. If later work requires a thickness-invariant interior field, bulk-like band edges or resolution finer than the present targets, a thicker series can become useful. Likewise, an accuracy question dominated by bulk k sampling or final-geometry vacuum sensitivity should be addressed by the relevant controlled numerical refinement, not automatically by increasing slab thickness. No additional run is required to resolve the present production-model selection from the supplied evidence.

## Final Stage 02 decision

**VERDICT: PASS WITH DOCUMENTED LIMITATION.**

Stage 02 is closed **for the intended clean β(001) screening scope** using **4L FULL**. The validated comparisons are:

1. genuine electronic/ionic completion of the accepted endpoints;
2. 4L→6L convergence of paired surface excess at the independent matched bulk reference, within 1 meV/Å²;
3. face-specific work functions and the defined vacuum-aligned frontier descriptors within 0.05 eV;
4. stable observed layer/coordination/termination topology in this lateral cell;
5. insensitivity of these screening observables to the tested central-layer constraint, with FULL retained.

The full interior electrostatic profile remains thickness-dependent. That limitation does not prevent progression for the measured surface-screening observables, but prevents calling the entire slab electrostatically or bulk-interior converged. The 2L failures support selecting 4L for these observables.

## Remaining scope limitations

- Rigorous bulk-like VBM/CBM, layer-resolved state localization and a semiconductor gap are unvalidated; the frontier separation is an index-based diagnostic with partial occupations.
- The basal reference assumes ISPIN=1. No spin-state competition is established by these nonmagnetic runs.
- Absolute Γ_pair inherits finite smearing, bulk/basis sampling, vacuum and relaxation tolerances. The reported E₀/F/without-entropy robustness is not a full systematic-error bound, and individual inequivalent-face surface energies are not resolved.
- Adsorption energies, reaction intermediates, reaction barriers and adsorption coverage require their own convergence evidence.
- Interfaces/cocatalysts, larger-cell reconstruction, other facets/terminations and production hybrid-functional slabs are outside this validation.
- Earlier setup-validation endpoint audits do not automatically transfer to these independently relaxed production endpoints; no tighter independent force audit of the present endpoints is claimed.

Only this Markdown report was written inside the repository. Native outputs and inputs were read without modification; temporary work stayed outside the repository. CMW records were read without submitting, cancelling, releasing or changing a job. No calculation, staging, commit or push was performed.
