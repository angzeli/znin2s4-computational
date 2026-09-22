# WP2 Stage 02 — β(001) Relaxed-Region Convergence

## Executive conclusion

**J37.1 passes numerical convergence, and the relaxed-region test is a PASS for the measured WP2 screening observables.** Fixing the central two complete layers preserves the surface coordination topology, changes surface bond lengths by no more than 0.002241 Å, and satisfies both predefined project targets: the paired excess-energy change is **+0.117208 meV/Å²** using E₀, versus a 1 meV/Å² target; the largest face-specific work-function change is **0.006344 eV** and the largest aligned frontier-extremum change is **0.032970 eV**, versus a 0.05 eV target. The interior potential remains similar, although it is not identical.

This is a screening-level agreement, not a claim of force equivalence or an unstrained fixed interior. The outer layers move outward relative to FULL by approximately 0.044–0.046 Å, and the fixed atoms carry residual forces up to **0.100972 eV/Å**. Those forces diagnose real suppressed relaxation, but are not accompanied by a reconstruction change, large bond distortion, energetic failure or major electrostatic change in the tested observables. **Retain FULL relaxation for production:** it imposes fewer assumptions and, in these runs, completed sooner. FIXED is an informative convergence control, not an improvement in the production model.

## Compared calculations and provenance

| Model | Job | Native calculation directory | Ionic degrees of freedom |
| --- | --- | --- | --- |
| FULL | J34.1 | [4L thickness reference](../../../calculation/02_numerical_convergence/03_slab_thickness/01_relaxation/beta001_4L) | All 28 atoms free |
| FIXED | J37.1 | [4L relaxed-region test](../../../calculation/02_numerical_convergence/04_relaxed_region/beta001_4L) | Outer 14 atoms free; central 14 fixed |

Source POSCAR coordinates and lattice vectors are exactly equal when parsed, with identical ordering: four repeated `S₄ In₂ Zn₁` units, totaling **Zn₄In₈S₁₆**. The substantive input difference is only Selective Dynamics. POSCAR title and INCAR `SYSTEM` labels differ harmlessly. A complete parsed INCAR comparison and the actual XML parameter comparison find no other setting difference. KPOINTS and POTCAR are byte-identical, verified by SHA-256; the PAW sequence is `PAW_PBE S 06Sep2000`, `In_d 06Sep2000`, `Zn 06Sep2000`, repeated four times to match the atom ordering.

The common cell vectors, in Å, are `(3.8782842923, 0, 0)`, `(−1.9391421462, 3.3586927203, 0)` and `(0, 0, 65.8651576100)`. The one-face area is **13.0259652198 Å²**, and the volume is **857.9572522266 Å³**. Both cells remain fixed during relaxation.

The matched model is VASP 6.6.1, PBE (`GGA=PE`), D3(BJ) (`IVDW=12`), `ENCUT=500 eV`, Γ-centred 10×10×1 sampling, `PREC=Accurate`, `LREAL=False`, `LASPH=True`, `ISPIN=1`, `ISYM=0`, Gaussian `ISMEAR=0`, `SIGMA=0.05 eV`, and `LDIPOL=True`, `IDIPOL=3`, `DIPOL=(0.5,0.5,0.5)`. Both use `ALGO=Normal` (executed `IALGO=38`), `AMIN=0.01`, effective `AMIX=0.40`, `BMIX=1.00`, `EDIFF=10⁻⁶ eV`, `NELM=120`, `NELMIN=2`, and fresh `ISTART=0`, `ICHARG=2`. Ionic controls are `IBRION=2`, `ISIF=2`, effective `POTIM=0.5`, `NSW=100`, and `EDIFFG=−0.02 eV/Å`. Output settings (`LCHARG`, `LVHAR`, `LORBIT=11`, `LWAVE=False`) also match. Both ran with eight MPI ranks × one thread, `NCORE=4`, `KPAR=1`, under the same managed launcher and VASP binary identity recorded in RUN_METADATA.

Energies, forces and evaluated geometries below come from complete `vasprun.xml`, checked against OUTCAR. CONTCAR supplies high-precision endpoint coordinates; LOCPOT/CHGCAR headers match those endpoints within their printing precision. EIGENVAL supplies frontier eigenvalues, cross-checked against XML. Runtime records identify J34.1 and the J37.1 attempt `ac28454671a94f51a592371e10442d2f`. No calculation or job was altered for this analysis.

## Numerical convergence

Both runs have an explicit EDIFF termination at every ionic evaluation, with both final `|dE|` and `|d eps| < 10⁻⁶ eV`. No cycle reaches its 120-step ceiling. Each XML electronic-cycle count agrees with the parsed DAV sequence; force-block counts agree with ionic-evaluation counts. J37.1 has **28/28 converged cycles**, a genuine structural-convergence marker and a normal timing footer. It stopped by the force criterion, before NSW or walltime exhaustion. All free force components satisfy 0.02 eV/Å; the stricter maximum vector norm is **0.012752 eV/Å**. The 14 constrained atoms are excluded from that stopping test.

OUTCAR, standard output and standard error contain no BRMIX, non-Hermitian, NaN/Inf, EDDDAV/diagonalisation failure, charge-sloshing warning, internal error or abnormal-termination indication. VASP's generic LREAL efficiency advice is not a numerical failure; the approved reciprocal-space projection was retained.

| Quantity | J34.1 FULL | J37.1 FIXED |
| --- | ---: | ---: |
| Ionic evaluations | 14 | 28 |
| First SCF iterations | 75 | 76 |
| Subsequent SCFs: range; median | 4–19; 6 | 4–17; 6 |
| Total electronic iterations | 169 | 263 |
| VASP elapsed time | 3 h 29 min 28 s | 4 h 57 min 54 s |
| Initial maximum force on active atoms (eV/Å) | 0.358748 | 0.359412 |
| Final maximum force on active atoms (eV/Å) | 0.018025 | 0.012752 |
| Largest single ionic displacement (Å) | 0.017307 | 0.017339 |

The constrained calculation does **not** relax faster: it uses twice as many ionic evaluations, 55.6% more electronic iterations and 42.2% more elapsed time. The costly first SCF is followed by much cheaper propagated SCFs in both models; their subsequent medians are both six iterations. The same starting geometry nevertheless takes 75 versus 76 initial iterations, with an initial free-energy difference of 0.00001127 eV and maximum force-vector difference of 0.000783 eV/Å. These small numerical differences are not evidence of a changed physical input model or exact trajectory reproducibility.

Both free-energy trajectories decrease at every evaluated geometry. Forces are nonmonotonic: FIXED has rebounds around evaluations 11, 18 and 22–27, ending at **0.028677→0.033065→0.012752 eV/Å** over its last three evaluations. FULL ends at **0.033430→0.028789→0.018025 eV/Å**. These are converged electronic force evaluations, not NELM-limited optimization steps. Runtime is descriptive and is not a convergence criterion.

### J37.1 electronic audit by ionic evaluation

Every row below has an explicit OUTCAR EDIFF marker. `rms` is the terminal DAV value; `rms(c)` is the last available value in that cycle because the terminating row omits it. Residuals retain VASP's native normalization. Neither an exit code nor a force block is used as a substitute for the electronic termination check.

| Ionic evaluation | SCF steps | Final dE (eV) | Final d eps (eV) | Final rms | Late rms(c) | Max free force (eV/Å) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 76 | 2.22180e-07 | -9.38270e-07 | 1.160e-03 | 2.780e-03 | 0.359412 |
| 2 | 13 | 1.75110e-07 | -1.66340e-07 | 3.780e-04 | 1.160e-02 | 0.113622 |
| 3 | 17 | 8.35160e-07 | -7.44770e-08 | 8.870e-04 | 5.460e-04 | 0.123796 |
| 4 | 10 | 6.18300e-07 | -5.32190e-08 | 2.460e-04 | 1.570e-04 | 0.079730 |
| 5 | 10 | 3.03790e-07 | -1.49920e-07 | 3.560e-04 | 2.750e-04 | 0.075732 |
| 6 | 6 | -8.30750e-07 | -5.04690e-07 | 6.610e-04 | 1.170e-03 | 0.040692 |
| 7 | 6 | -3.47350e-07 | -3.63940e-07 | 6.050e-04 | 9.010e-04 | 0.039486 |
| 8 | 6 | 4.38830e-07 | -2.95310e-07 | 7.220e-04 | 9.790e-04 | 0.020741 |
| 9 | 5 | 9.02280e-07 | -9.02190e-08 | 5.410e-04 | 9.310e-04 | 0.025750 |
| 10 | 5 | 2.17930e-07 | -3.73970e-07 | 6.950e-04 | 1.940e-03 | 0.024334 |
| 11 | 6 | 1.82020e-07 | -4.09170e-07 | 5.150e-04 | 9.150e-04 | 0.039543 |
| 12 | 4 | 3.75950e-07 | -3.09310e-07 | 3.970e-04 | 1.700e-03 | 0.035084 |
| 13 | 6 | -6.91960e-07 | -5.06340e-07 | 5.120e-04 | 1.560e-03 | 0.023136 |
| 14 | 7 | -5.64370e-07 | -4.54900e-07 | 5.900e-04 | 9.360e-04 | 0.026490 |
| 15 | 5 | -1.28800e-07 | -1.60880e-07 | 3.850e-04 | 2.010e-03 | 0.022084 |
| 16 | 7 | -7.54160e-07 | -5.53250e-07 | 5.550e-04 | 6.930e-04 | 0.027374 |
| 17 | 7 | 1.09860e-07 | -1.91690e-07 | 3.740e-04 | 3.060e-04 | 0.025464 |
| 18 | 7 | 1.42360e-07 | -1.73450e-07 | 3.470e-04 | 2.970e-04 | 0.048021 |
| 19 | 4 | 5.73550e-07 | -3.14680e-07 | 4.300e-04 | 1.730e-03 | 0.045851 |
| 20 | 6 | 1.44620e-07 | -4.67550e-07 | 6.100e-04 | 1.890e-03 | 0.042137 |
| 21 | 7 | -3.29040e-07 | -5.49290e-07 | 6.580e-04 | 6.270e-04 | 0.035744 |
| 22 | 8 | 3.43690e-07 | -1.52720e-07 | 4.360e-04 | 1.780e-03 | 0.041317 |
| 23 | 7 | -9.82310e-07 | -9.09100e-07 | 7.210e-04 | 1.050e-03 | 0.025999 |
| 24 | 4 | 7.09890e-07 | -3.98370e-07 | 5.510e-04 | 1.950e-03 | 0.038088 |
| 25 | 7 | -8.42870e-07 | -1.43600e-07 | 4.050e-04 | 7.210e-04 | 0.029205 |
| 26 | 6 | 1.31370e-07 | -4.84830e-07 | 5.860e-04 | 1.290e-03 | 0.028677 |
| 27 | 4 | 8.22380e-07 | -3.77030e-07 | 4.430e-04 | 1.770e-03 | 0.033065 |
| 28 | 7 | -3.89000e-07 | -5.63750e-07 | 5.210e-04 | 9.700e-04 | 0.012752 |

CONTCAR matches the last evaluated J37.1 geometry within **3.3×10⁻⁷ Å**. Fixed-atom positions remain within **2.5×10⁻⁷ Å** of their originals throughout the XML trajectory, consistent with coordinate printing precision. The smallest pair distance encountered is **2.295437 Å** (FULL: **2.294485 Å**), and the largest single-step movement is only 0.017339 Å. Layer ordering is preserved at every evaluation; no overlap, slab collapse or transfer between layer units is observed.

## Constraint definition and residual forces

Layers are numbered by increasing z. Each comprises seven atoms with ordering `S, S, S, S, In, In, Zn`. There are no partial-direction constraints.

| Layer | Atom indices (1-based) | Starting z range (Å) | J37.1 constraint |
| --- | ---: | ---: | ---: |
| L1 | 1–7 | 10.000000–19.419858 | T T T |
| L2 | 8–14 | 22.148433–31.568291 | F F F |
| L3 | 15–21 | 34.296866–43.716724 | F F F |
| L4 | 22–28 | 46.445299–55.865158 | T T T |

The following are **physical, unmasked forces**: raw OUTCAR TOTAL-FORCE agrees with XML to printed precision, including the nonzero constrained forces. They are not the masked optimization vector. Values below are in eV/Å; “mean norm” means the average magnitude, and “mean Fz” retains its sign.

| J37.1 fixed group | Mean norm | RMS norm | Maximum norm | max \|Fx\| | max \|Fy\| | max \|Fz\| | Mean Fz |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| L2 | 0.071136 | 0.074243 | 0.100972 | 2.768e-05 | 1.562e-05 | 0.100972 | −0.003102 |
| L3 | 0.056809 | 0.061696 | 0.082133 | 2.767e-05 | 1.555e-05 | 0.082133 | +0.002865 |
| All 14 | 0.063972 | 0.068259 | 0.100972 | 2.768e-05 | 1.562e-05 | 0.100972 | −0.000119 |

L2 has a small net force toward lower z (**−0.021717 eV/Å summed over the layer**); L3 has a similar opposing net force (**+0.020056 eV/Å**). Individual z forces alternate in sign and are much larger than these means, while lateral components are tiny. The constraint therefore suppresses internal relaxation as well as a small tendency for the fixed layers to separate; the fixed region is not “fully force-free.”

For the corresponding 14 central atoms in FULL, the maximum and RMS force norms are only **0.010749 and 0.004826 eV/Å**, versus **0.100972 and 0.068259 eV/Å** in FIXED. This is a real mechanical difference, not a failed constrained stopping criterion. It makes a strain-free/bulk-anchor interpretation inappropriate. Its observable consequences are assessed below rather than inferred from a force threshold that applies only to free coordinates.

## Structural comparison

Endpoint displacements are **FIXED minus FULL**, using periodic minimum-image vectors in the common triclinic cell. No rigid translation, rotation or fitted z shift is removed. For the 14 outer atoms, RMS and maximum positional differences are **0.045495 and 0.052379 Å**; RMS and maximum absolute z differences are **0.045495 and 0.052379 Å**. Almost all of the difference is slab-normal.

| Face / layer | RMS position difference (Å) | Maximum difference (Å) | Mean signed Δz (Å) | Differential z-rumpling RMS (Å) |
| --- | ---: | ---: | ---: | ---: |
| Lower / L1 | 0.046539 | 0.049530 | −0.046352 | 0.004169 |
| Upper / L4 | 0.044426 | 0.052379 | +0.044176 | 0.004710 |

“Differential z-rumpling RMS” is the standard deviation of the seven atomwise Δz values after subtracting their layer mean, used only to distinguish internal distortion from layer translation. It is not a rigid fit applied to the positional comparison. Each atomic plane has one atom in this lateral cell, so an independent in-plane corrugation amplitude is not resolved. The opposite signed surface shifts represent a slightly thicker FIXED slab, not a removable common translation.

| Geometric quantity (Å) | FULL | FIXED | FIXED − FULL |
| --- | ---: | ---: | ---: |
| L1 seven-plane thickness | 9.447246 | 9.459722 | +0.012476 |
| L2 seven-plane thickness | 9.387962 | 9.419858 | +0.031897 |
| L3 seven-plane thickness | 9.390270 | 9.419858 | +0.029588 |
| L4 seven-plane thickness | 9.431372 | 9.445509 | +0.014137 |
| L1–L2 interlayer gap | 2.758836 | 2.789247 | +0.030410 |
| L2–L3 interlayer gap | 2.757570 | 2.728575 | −0.028995 |
| L3–L4 interlayer gap | 2.752280 | 2.764676 | +0.012396 |
| Total slab thickness | 45.925536 | 46.027445 | +0.101909 |
| Atom-free vacuum gap | 19.939621 | 19.837712 | −0.101909 |

Layer thickness is the difference between its outermost atomic z coordinates, and each interlayer gap is the separation between adjacent layer envelopes. The largest gap effects are the **+0.030410 Å L1–L2 widening** and **−0.028995 Å L2–L3 narrowing**. Total thickness increases by 0.101909 Å; at the fixed cell height, the vacuum gap decreases by the same amount. These changes are part of the response to the constraint, not an independently changed cell setting.

### Bonds, coordination and reconstruction identity

Periodic S–In/Zn neighbors were matched by atom identity and periodic image after minimum-image registration to the common starting structure. A 3.0 Å distance cutoff gives the same **56 cation–S links per cell** in the original structure, FULL and FIXED; this identity persists at both 2.8 and 3.2 Å. Thus the result is insensitive to those cutoff choices. The table contains intralayer bonds in each outer layer; repeated periodic neighbors are counted separately.

| Surface bond family | Bonds per outer layer | FULL range (Å) | FIXED range (Å) | Mean Δlength (Å) | Largest \|Δlength\| (Å) |
| --- | ---: | ---: | ---: | ---: | ---: |
| L1: In-S | 10 | 2.437845–2.726087 | 2.437708–2.726617 | +0.000492 | 0.001022 |
| L1: Zn-S | 4 | 2.294485–2.571293 | 2.296121–2.573137 | +0.001689 | 0.001844 |
| L4: In-S | 10 | 2.440303–2.716003 | 2.442544–2.716365 | +0.001069 | 0.002241 |
| L4: Zn-S | 4 | 2.298840–2.584817 | 2.299723–2.587044 | +0.001218 | 0.002227 |

Nearest-neighbor distances are **2.294485→2.296121 Å** on L1 and **2.298840→2.299723 Å** on L4 (FULL→FIXED). Both outer layers retain the neighbor-count sequence **3, 4, 4, 3, 4, 6, 4** in `S₄In₂Zn` atom order. There is no gained/lost first-shell cation–S link, atom exchange, altered termination identity, bond healing/breaking or major internal rumpling. The surface layers primarily shift outward with small internal distortions. This establishes preservation of the observed topology within the present lateral cell, not the absence of larger-cell reconstructions or a universal chemical coordination assignment.

## Energetic comparison

Endpoint energy definitions are kept separate: F is finite-smearing TOTEN, E_without_entropy omits the entropy contribution, and E₀ is VASP's `energy(sigma→0)` estimate. OUTCAR and XML agree within their printed precision. Positive differences mean FIXED is higher in energy.

| Quantity | FULL (eV/cell) | FIXED (eV/cell) | ΔE (eV/cell) | ΔΓ_pair (meV/Å²) |
| --- | ---: | ---: | ---: | ---: |
| F / TOTEN | −114.45163599 | −114.45019249 | +0.00144350 | +0.110817 |
| E_without_entropy | −114.45034546 | −114.44873545 | +0.00161001 | +0.123600 |
| E₀ (primary) | −114.45099072 | −114.44946397 | +0.00152675 | +0.117208 |

With identical composition, parent phase, numerical model and cell, the **same bulk reference cancels** in this matched difference:

`ΔΓ_pair = (E_FIXED − E_FULL) / A`, with `A = |a × b| = 13.025965219827203 Å²`.

The denominator is the **one-face area A**, not 2A: this is the paired surface excess, not a separate energy assigned to either inequivalent face. Using E₀ gives **+0.117208 meV/Å²**, comfortably meeting the predefined `|ΔΓ_pair| ≤ 1 meV/Å²` refinement target in [ROADMAP §5.6.6](../../../ROADMAP.md#566-convergence-acceptance). F and E_without_entropy give the same pass. The sign is consistent with suppression of available ionic relaxation; no claim that the constrained state is a lower variational minimum is made. This matched cancellation does not supply the independent bulk reference missing from the separate slab-thickness surface-energy analysis.

## Vacuum-referenced electronic comparison

### Face-specific vacuum levels and work functions

Both final LOCPOTs contain the same LVHAR Hartree-plus-ionic potential. Planar means use native z grids, with **6 Å setback from each outermost atom and 1 Å exclusion from each cell boundary**, matching the [thickness analysis](../03_slab_thickness/THICKNESS_ANALYSIS.md). Native grids are 56×56×960 and match CHGCAR. Integrated charge reproduces 248 electrons to better than 10⁻³ electrons; the electron density is obtained by dividing CHGCAR by cell volume.

The established vacuum screening requires width ≥2 Å, |slope| ≤0.005 eV/Å, maximum detrended residual ≤0.02 eV, potential range ≤0.03 eV and maximum absolute density ≤10⁻⁵ e/Å³. All four windows pass. Their positions and diagnostics are:

| Run / side | Actual z interval (Å) | Slope (eV/Å) | Potential range (eV) | Max detrended residual (eV) | Max \|density\| (e/Å³) |
| --- | ---: | ---: | ---: | ---: | ---: |
| FULL / lower | 1.029143–3.910744 | -0.000344077 | 0.000987015 | 9.4972e-05 | 7.74804e-07 |
| FULL / upper | 61.885804–64.836015 | 0.000393492 | 0.00118307 | 0.000125663 | 9.05583e-07 |
| FIXED / lower | 1.029143–3.842134 | -0.000202286 | 0.000596256 | 7.53262e-05 | 5.47346e-07 |
| FIXED / upper | 61.954414–64.836015 | -5.00382e-05 | 0.00024004 | 0.000112768 | 1.79193e-06 |

No abrupt potential jump lies inside a selected interval: the largest adjacent-grid change among them is 0.000075 eV, well below the discontinuity across the periodic vacuum boundary. Changing the setback from 6 to 5 Å changes each reference by at most **0.000254 eV**; both faces still pass screening. The chosen windows are depleted, flat and separated from the dipole discontinuity. These checks validate reference extraction, not every possible periodic-image energy effect.

The raw vacuum and Fermi values below document **same-run references only**; differences between raw absolute offsets are not interpreted as physical shifts. Work functions use `Φ_side = V_vac,side − E_F` independently for each face.

| Quantity (eV) | FULL | FIXED | Meaningful FIXED − FULL difference |
| --- | ---: | ---: | ---: |
| Lower vacuum, raw reference | 5.976712 | 5.980947 | — |
| Upper vacuum, raw reference | 9.139884 | 9.136823 | — |
| Fermi level, raw reference | 2.209001 | 2.212284 | — |
| Lower-face Φ | 3.767710 | 3.768663 | +0.000953 |
| Upper-face Φ | 6.930883 | 6.924539 | −0.006344 |
| Upper-minus-lower vacuum step | 3.163173 | 3.155876 | −0.007296 |

**Both face-specific changes pass 0.05 eV.** The references are not averaged. The vacuum-step change is only −0.007296 eV, so the large face asymmetry persists with a small constraint-induced perturbation. As a consistent supporting diagnostic, OUTCAR's final z dipole moment is **0.227717→0.227123 e·Å**, a change of −0.000594 e·Å under the same cell/origin convention. This is a slab-dipole diagnostic, not a Berry-phase bulk polarization calculation.

### Frontier-state screening diagnostic

EIGENVAL and XML arrays agree within the XML's four-decimal printing precision. For ISPIN=1, these files use occupations in [0,1]; including spin degeneracy, `2 Σ_kb w_k f_kb` reproduces 248 electrons within 2.3×10⁻⁷ electrons. Both endpoints retain partial occupations. The common convention is **b = NELECT/2 = 124**, with the maximum of band 124 and minimum of band 125 taken over the sampled k points. This is an index-based screening diagnostic, not a layer-resolved state assignment.

| Frontier diagnostic (eV) | FULL | FIXED | FIXED − FULL |
| --- | ---: | ---: | ---: |
| Frontier separation min ε₁₂₅ − max ε₁₂₄ | 0.028068 | 0.056293 | +0.028225 |
| Band-124 maximum, relative to lower vacuum | −3.758227 | −3.791197 | −0.032970 |
| Band-125 minimum, relative to lower vacuum | −3.730159 | −3.734904 | −0.004745 |
| Band-124 maximum, relative to upper vacuum | −6.921399 | −6.947073 | −0.025674 |
| Band-125 minimum, relative to upper vacuum | −6.893331 | −6.890780 | +0.002551 |

The largest change among the four aligned extrema is **0.032970 eV**, within the predefined project-specific 0.05 eV working target. The frontier separation increases by 0.028225 eV, so relative changes of this small quantity should not be hidden by the absolute screening pass. Its 0.028068/0.056293 eV values are comparable to the 0.05 eV Gaussian smearing scale, and partial occupations persist. Neither value is a validated semiconductor band gap or a rigorous bulk VBM/CBM separation.

## Interior electrostatics

The method matches the slab-thickness audit. For each run, construct `V̄(z)=mean_xy V(x,y,z)` and a periodic piecewise-linear profile; its box average is `V_M(z) = (1/L) ∫[z−L/2,z+L/2] V̄(s) ds`, with the **fixed original repeat L = 12.1484331333 Å**. Compare the same local window `u=z−g ∈ [−L/4,+L/4]`, width **6.074217 Å**, where g is the central interlayer-gap midpoint: **32.921823 Å in FULL** and **32.932579 Å in FIXED**. The entire averaging footprint lies inside both slabs. Structural registration aligns equivalent gaps, not fitted potential curves or fitted slopes; the actual geometric shifts are retained in the structural analysis above.

Each profile is aligned separately to its own lower or upper vacuum. Central means and RMS differences are integrated on 2001 common u points, slopes are least-squares fits, and ranges are maximum minus minimum. Box integrals are independently checked by direct trapezoidal quadrature through the native knots. No slope is subtracted or artificial flattening applied. The central-interstitial descriptor averages the **inner half of the central gap**, excluding the quarters nearest its bounding planes; it does not sample ionic cores.

| Interior descriptor | FULL | FIXED | FIXED − FULL |
| --- | ---: | ---: | ---: |
| Central mean relative to lower vacuum (eV) | −8.624581 | −8.634568 | −0.009987 |
| Central mean relative to upper vacuum (eV) | −11.787754 | −11.790444 | −0.002690 |
| Central slope (eV/Å) | 0.029048 | 0.031864 | +0.002817 |
| Central potential range (eV) | 0.178816 | 0.193657 | +0.014841 |
| Interstitial mean relative to lower vacuum (eV) | −0.401007 | −0.437586 | −0.036579 |
| Interstitial mean relative to upper vacuum (eV) | −3.564180 | −3.593462 | −0.029282 |
| Central interstitial density (e/Å³) | 0.106420 | 0.108662 | +0.002242 |

| Aligned central-profile comparison | Lower-vacuum alignment | Upper-vacuum alignment |
| --- | ---: | ---: |
| Maximum pointwise \|difference\| (eV) | 0.018612 | 0.011315 |
| RMS profile difference (eV) | 0.011211 | 0.005760 |

The central slope increases by **0.002817 eV/Å (9.7%)**, and the central range increases by **0.014841 eV**. The fixed model is slightly steeper, not flatter or more bulk-like. The largest aligned central-profile discrepancy is 0.018612 eV, while central-interstitial density increases by **2.106%**. Local interstitial potentials change by −0.036579/−0.029282 eV under lower/upper alignment. Thus there is a measurable internal response, but no reversal of the gradient, major shift in the face-to-face step or large change of the aligned central environment. The two gauge choices lead to the same qualitative assessment.

The 0.05 eV criterion was predefined for screening-level electronic descriptors, not specifically calibrated for every macroscopic-potential statistic or an interior-slope test. The profile values are therefore supporting evidence, not an invented independent formal PASS threshold. A nonzero field in this asymmetric slab need not be a finite-size artifact. This comparison measures sensitivity to the imposed constraint; it does not establish the field's bulk limit, prove a bulk-like plateau, or show that suppressing interior relaxation makes the electrostatic solution more physical.

## Relaxed-region convergence decision

**VERDICT: PASS — for the tested WP2 screening observables at the specified numerical model and tolerances.**

| Decision component | Evidence | Assessment |
| --- | --- | --- |
| Numerical validity | 28/28 J37.1 SCFs reach EDIFF; genuine constrained ionic convergence | Pass |
| Qualitative surface structure | Same layer identities and first-shell coordination; no bond-topology change | Pass within this lateral cell |
| Surface geometry | ≤0.052379 Å atomwise difference, mainly layer shifts; ≤0.002241 Å bond change | Small geometric response, not identical structures |
| Paired excess energy | +0.117208 meV/Å² in E₀ | Meets 1 meV/Å² target |
| Face-specific electronic references | maximum absolute ΔΦ = 0.006344 eV; max aligned frontier change = 0.032970 eV | Meets 0.05 eV target |
| Interior solution | max aligned central-profile difference 0.018612 eV; slope +9.7% | Similar with measurable perturbation |
| Fixed-region mechanics | max/RMS residual force 0.100972/0.068259 eV/Å | Real imposed strain; not force-equivalent or strain-free |

The force diagnostic deserves weight: the constraint prevents the small central-layer displacements that reduce those physical forces in FULL. Nevertheless, central-atom endpoint differences are at most **0.025846 Å**, both surfaces retain their topology, and the measured energetic/electronic effects are well within the screening targets. Taken together, the evidence supports a limited screening PASS rather than severe constraint-induced distortion of the tested model. It does not qualify FIXED for studies whose observable is the interior force field, strain response or a finer-than-tested structural displacement.

**Production recommendation: keep all four layers relaxed.** The constrained control provides no observed cost advantage, retains avoidable residual forces and offers no evidence of a more physical interior. This test resolves sensitivity to this particular relaxed-region choice for the measured observables; it does not mandate a central fixed region or require a new calculation. No 6L constrained test was launched.

## Remaining caveats

- The structural pass is supported by small bond changes and unchanged topology, not a previously specified universal positional-displacement tolerance. Surface translations of approximately 0.045 Å and a 0.101909 Å thickness change remain real.
- The independent tighter endpoint audit from the earlier setup-validation campaign does not automatically validate either of these endpoints. Both present runs use EDIFF=10⁻⁶ eV and EDIFFG=−0.02 eV/Å; the small energy difference is not claimed more accurately than their demonstrated numerical model permits.
- The constraint comparison does not resolve the 4L→6L interior-field dependence from the thickness study. Neither FULL nor FIXED has thereby acquired a proven bulk-like interior.
- The matched bulk term cancels here, but the independent matched-bulk reference for thickness-dependent surface energetics remains a separate requirement. No absolute or individual-face surface energy is assigned.
- Frontier extrema remain screening diagnostics; rigorous bulk-like VBM/CBM and a semiconductor gap are unvalidated. Smearing, localization and larger-cell reconstruction questions are not resolved by this comparison.
- Adsorption energies, reaction energetics and interface/cocatalyst models require their own convergence checks. Agreement of clean-slab observables cannot be transferred automatically.

All differences were calculated from unrounded source values and rounded only for display. Native calculation files were read without modification. Temporary analysis files remained outside the repository; only this report was written in the repository. No VASP calculation or CMW job was launched or altered, and no staging, commit or push was performed.
