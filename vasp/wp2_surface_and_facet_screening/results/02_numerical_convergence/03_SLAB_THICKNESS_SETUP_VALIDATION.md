# Setup Validation for β-ZnIn₂S₄ Slab Calculations

## 1. Purpose and evidence scope

**The tested 4L electronic and relaxation protocol passes setup validation:** `ALGO=Normal`, `AMIN=0.01`, and `EDIFF=10⁻⁶ eV` produce electronically converged force evaluations throughout a full relaxation; a fresh static audit at `EDIFF=10⁻⁷ eV` independently retains the maximum-force criterion of **0.02 eV/Å**. This qualifies the protocol for the forthcoming thickness study. It does not establish that 4L is a converged slab thickness.

Process completion, electronic self-consistency and structural convergence are assessed separately. A normal exit alone does not qualify energies or forces, and an electronically unconverged force evaluation cannot support production ionic motion. The objective is a reproducible, usable protocol rather than termination alone.

The [setup-validation archive](../../calculation/02_numerical_convergence/03_slab_thickness/00_setup_validation) contains eight input sets. Native results are present for **E3, E2x, R4L, P4L and F4L**. **E0, E1 and E2 contain inputs only:** no OUTCAR, OSZICAR, completed-run metadata, parsed trajectories or analysis records are present in those three directories. No existing figures or trajectory CSVs are present anywhere in this archive. Consequently, the later validation chain can be reconstructed quantitatively, but the requested retrospective E0/E1/E2 performance comparison cannot. Their historical convergence claims and E1→E2 improvement factors are not treated as independently verified results here.

## 2. Starting problem: canonical slab convergence

The archived baseline inputs specify static fresh-start calculations on a 14-atom 2L slab (E0) and a 28-atom 4L slab (E1), with `NELM=60`, `ALGO=Normal`, `EDIFF=10⁻⁶ eV`, `ISTART=0` and `ICHARG=2`. They omit explicit AMIN, AMIX and BMIX overrides. The intended canonical mixing is AMIN=0.10, AMIX=0.40 and BMIX=1.0; execution of that exact default configuration in E0/E1 cannot be checked without their outputs. AMIN=0.10 is explicitly recorded for the archived E3 calculation, which also leaves AMIN unspecified, while AMIX=0.40 and BMIX=1.0 are recorded for the validated AMIN=0.01 runs.

The baseline cell heights are 41.56829 Å for 2L and 65.86516 Å for 4L. The archived E3 output explicitly warns that a lattice vector exceeds 50 Å and that a large AMIN may permit charge sloshing; it suggests trying 0.01 if convergence is problematic. That warning identifies a numerical risk, not a diagnosis proving that slab thickness or dipole correction caused the instability. No 6L result is present in this archive, and no quantitative 6L comparison is made.

AMIN controls the lower bound of the Kerker initial mixing approximation used by the Broyden/Pulay scheme. Reducing that bound is a physically motivated numerical intervention for slowly converging long-wavelength density modes. VASP also identifies smaller AMIN as a possible remedy for slow dipole-corrected slab convergence. These considerations motivate the test; they do not establish the mechanism of this particular trajectory. See the official [AMIN definition](https://vasp.at/wiki/AMIN) and [dipole-correction guidance](https://vasp.at/wiki/LDIPOL).

The original baseline diagnosis—oscillatory, bounded-budget nonconvergence—cannot be quantified from the retained E0/E1 inputs. The available results nevertheless show why extending a budget alone is not a sufficient qualification rule: E3 terminates normally with severely amplified residuals, whereas the selected recipe genuinely reaches EDIFF in E2x and remains usable after ionic motion. More iterations are justified by the behavior of a candidate recipe, not by an assumption that any finite trajectory will eventually converge.

## 3. Controlled SCF qualification

Input comparisons establish the experimental design exactly. E1→E2 adds only `AMIN=0.01`; E2→E2x changes only `NELM=60` to `120`; E1→E3 adds `AMIX=0.2` and `BMIX=0.0001`. POSCAR, KPOINTS and the locally retained POTCAR are byte-identical within each of these comparisons. E3 is an alternative mixing recipe, not a one-parameter causal comparison with E2.

| Test | Model | Intervention / budget | Electronic result supported by this archive | Decision |
| --- | --- | --- | --- | --- |
| [E0](../../calculation/02_numerical_convergence/03_slab_thickness/00_setup_validation/E0_r1) | 2L, 14 atoms | Canonical inputs; 60-step ceiling | Outputs absent; actual iteration count, residuals and EDIFF status unverified | Baseline design retained; no convergence claim |
| [E1](../../calculation/02_numerical_convergence/03_slab_thickness/00_setup_validation/E1_r1) | 4L, 28 atoms | Canonical inputs; 60-step ceiling | Outputs absent; stronger oscillation relative to E0 cannot be quantified | Baseline design retained; no convergence claim |
| [E2](../../calculation/02_numerical_convergence/03_slab_thickness/00_setup_validation/E2_r1) | Same fixed 4L geometry | AMIN=0.01; 60-step ceiling | Outputs absent; reported near-threshold behavior and initial nonconvergence cannot be rechecked | Candidate recipe subsequently tested by E2x |
| [E3](../../calculation/02_numerical_convergence/03_slab_thickness/00_setup_validation/E3_r1/OUTCAR) | Same fixed 4L geometry | AMIX=0.2, BMIX=0.0001; effective AMIN=0.10; 60 steps | 60 iterations, no EDIFF termination; strong late amplification and numerical warnings | Rejected |
| [E2x](../../calculation/02_numerical_convergence/03_slab_thickness/00_setup_validation/E2x_r1/OUTCAR) | Same fixed 4L geometry | Exact E2 recipe; ceiling raised to 120 | Genuine EDIFF at iteration **69** | Fixed-geometry qualification passed |

### E0–E2: controlled design, incomplete historical result archive

Neither E0 nor E1 is described as converged. Their inputs cannot establish failure either: the absent native trajectories prevent verification of the reported 60-step outcomes, comparison of 2L and 4L residuals, or assignment of a specific oscillatory regime. Likewise, no numerical E1→E2 reduction factors are recoverable from these records. The first 60 steps of E2x are not substituted for the missing E2 trajectory; agreement between two fresh-start nonlinear SCF trajectories must be measured, not assumed.

Selection of AMIN=0.01 is therefore supported here by its **demonstrated sufficiency in E2x, R4L, P4L and F4L**, together with rejection of the archived E3 alternative. Its claimed orders-of-magnitude improvement over E1 remains an archival evidence gap rather than a quantified result of this document.

### E3: severe amplification, despite normal process completion

The [E3 electronic trajectory](../../calculation/02_numerical_convergence/03_slab_thickness/00_setup_validation/E3_r1/OSZICAR) ends after 60 iterations without an EDIFF marker. The final `dE` is −55,767 eV and `d eps` is −25,772 eV; the final wavefunction residual is 73.4, and the last printed charge residual is 38.7 at iteration 59. These are failed SCF diagnostics, not physically interpretable slab energetics.

The following window statistics use absolute values for medians and signed values for peak-to-peak amplitudes. Energy changes are in eV; `rms` and `rms(c)` retain their native VASP diagnostic normalization.

| E3 electronic steps | Median \|dE\| (eV) | Median \|d eps\| (eV) | Median rms | Median rms(c) | Peak-to-peak dE (eV) | Peak-to-peak d eps (eV) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 31–40 | 153.38 | 11.2809 | 2.205 | 14.20 | 1,180.78 | 37.5909 |
| 41–50 | 298.845 | 84.6115 | 4.210 | 15.85 | 1,190.64 | 479.057 |
| 51–60 | 11,799.5 | 5,954.15 | 27.25 | 20.60 | 146,261 | 39,432.1 |

The final window contains nine reported charge residuals because iteration 60 omits `rms(c)`. The [stdout log](../../calculation/02_numerical_convergence/03_slab_thickness/00_setup_validation/E3_r1/stdout.log) contains BRMIX charge-density inconsistency messages, repeated non-Hermitian subspace-matrix warnings, and an explicit NELM nonconvergence warning. Thus the trajectory exhibits severe late amplification rather than a stable near-threshold plateau. A normal timing footer and launcher exit status zero do not change its scientific failure.

### E2x: genuine fixed-geometry convergence of the selected recipe

E2x reaches the explicit OUTCAR EDIFF termination at electronic iteration **69**, below its 120-step ceiling. Its final `dE=−1.0896×10⁻⁸ eV` and `d eps=−1.6641×10⁻⁷ eV` both satisfy `EDIFF=10⁻⁶ eV`. Final `rms=1.15×10⁻³`; the last printed `rms(c)=2.81×10⁻³` belongs to iteration 68, not 69. VASP elapsed time is **4,522.753 s (75.38 min)**.

| E2x electronic steps | Median \|dE\| (eV) | Median \|d eps\| (eV) | Median rms | Median rms(c) | Peak-to-peak dE (eV) | Peak-to-peak d eps (eV) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 31–40 | 2.038×10⁻⁴ | 1.312×10⁻⁵ | 0.006510 | 0.04175 | 1.240×10⁻³ | 5.602×10⁻⁵ |
| 41–50 | 1.652×10⁻⁴ | 9.146×10⁻⁶ | 0.006865 | 0.02430 | 1.609×10⁻³ | 2.157×10⁻⁵ |
| 51–60 | 1.323×10⁻⁴ | 1.449×10⁻⁵ | 0.007555 | 0.01975 | 1.029×10⁻³ | 4.184×10⁻⁵ |
| 61–69 | 1.867×10⁻⁵ | 2.737×10⁻⁶ | 0.002510 | 0.004970 | 8.928×10⁻⁵ | 8.153×10⁻⁶ |

The last row contains eight charge residuals. The trajectory is **oscillatory with late damping**, not monotonic: iteration 61 increases `|dE|` relative to 60, and some intermediate windows show renewed wavefunction or band-energy residual growth. Nevertheless, the final window improves jointly and the actual stopping condition is reached. An increase in NELM permitted completion of this particular recipe without changing its physical model. The resulting maximum force, 0.358811 eV/Å, is far above the structural threshold; static electronic success is not a relaxed structure.

The final available default-mixing–dielectric spectra also differ substantially:

| Run | Spectrum size | Printed mean | Eigenvalue range | Standard deviation |
| --- | ---: | ---: | --- | ---: |
| E3 | 44 | 0.0896 | 0.0016–0.2644 | 0.0793 |
| E2x | 43 | 1.2486 | 0.1345–3.3759 | 0.9247 |
| F4L | 43 | 1.2669 | 0.0590–2.9138 | 0.7977 |

These are the last printed spectra, not necessarily quantities evaluated on the terminating iteration. E3's spectrum lies far below unity and accompanies a demonstrably unstable trajectory; the accepted recipe's spectra remain broad. Neither the mean nor spectral width is a convergence certificate or proof of a unique failure mechanism.

## 4. Moving-geometry validation: R4L

A converged SCF solution at one geometry does not establish that force evaluations remain reliable after atomic motion. [R4L](../../calculation/02_numerical_convergence/03_slab_thickness/00_setup_validation/R4L_r1/OUTCAR) tests this explicitly using the E2x initial geometry and electronic recipe, with `IBRION=2`, `NSW=3`, fixed cell and `NELM=120` per electronic cycle. All three cycles have their own EDIFF termination marker; XML and OSZICAR iteration counts agree.

| Ionic evaluation | Initialization | Electronic iterations | Final dE (eV) | Final d eps (eV) | Maximum force (eV/Å) | Largest move from preceding evaluated geometry (Å) |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| 1 | Fresh electronic start | 97 | −3.7003×10⁻⁷ | −7.8882×10⁻⁷ | 0.360704 | — |
| 2 | After first ionic move | 19 | −1.2053×10⁻⁸ | −2.4635×10⁻⁸ | 0.112020 | 0.017401 |
| 3 | After second ionic move | 8 | −2.5389×10⁻⁷ | −1.7571×10⁻⁷ | 0.141429 | 0.003702 |

The expensive fresh-start cycle is followed by substantially cheaper propagated cycles. Two moved geometries establish short-run usability, although they are insufficient to predict long-run iteration costs. The force increase at evaluation 3 is not electronic failure: that cycle converges, while free energy decreases from −114.43259816 to −114.44378564 to −114.44405487 eV. Nearest periodic interatomic distances remain at least 2.30419 Å, and the cell is unchanged. No BRMIX, non-Hermitian, long-cell/charge-sloshing or other identified serious numerical warning appears in the retained native logs.

R4L therefore passes **short moving-geometry validation**, not structural relaxation: the final force remains above 0.02 eV/Å and there is no ionic-convergence marker. Its final evaluated geometry is the byte-identical source of P4L's POSCAR. The 97-iteration fresh cycle differs from E2x's 69 iterations on the same geometry, and their maximum atomic force-vector difference is **0.007953 eV/Å**. This measured sensitivity motivates an independent endpoint check rather than assuming that every EDIFF-converged initialization yields indistinguishable forces.

## 5. Full 4L relaxation: P4L

[P4L](../../calculation/02_numerical_convergence/03_slab_thickness/00_setup_validation/P4L_r1/OUTCAR) starts from R4L's final force-evaluated geometry with a fresh electronic initialization. It retains AMIN=0.01, `EDIFF=10⁻⁶ eV`, `NELM=120`, `IBRION=2`, `ISIF=2` and `EDIFFG=−0.02 eV/Å`, with `NSW=100` as its ceiling.

**All 16 ionic evaluations reach EDIFF**, and VASP explicitly reports `reached required accuracy - stopping structural energy minimisation`. The electronic iteration counts are:

`104, 7, 8, 7, 10, 3, 13, 3, 4, 4, 6, 6, 5, 6, 8, 6`.

This is **200 electronic iterations in total**: 104 for the initial fresh start and 96 across the 15 moved geometries. Propagated cycles require **3–13 iterations, median 6**, supporting the pattern of an expensive initial SCF followed by much cheaper electronic propagation. No cycle exhausts NELM, and there is no progressive increase in electronic cost. Native VASP elapsed time is **12,959.982 s (3.600 h)**; the local launch record includes **12,960.166 s**, a small difference attributable to the measurement boundary, not a conflicting scientific result.

| Ionic evaluation | Electronic iterations | Free energy (eV/cell) | Maximum force (eV/Å) |
| --- | ---: | ---: | ---: |
| 1 | 104 | −114.44410225 | 0.141328 |
| 7 | 13 | −114.45085738 | 0.034877 |
| 11 | 6 | −114.45105425 | 0.022923 |
| 15 | 8 | −114.45165473 | 0.029360 |
| 16 | 6 | −114.45183967 | **0.01818259** |

The final maximum force is below the retained threshold by **0.00181741 eV/Å**. The final energy without entropy is **−114.45052152 eV/cell**. Free energy decreases overall by **0.00773742 eV/cell** from the first P4L evaluation, with a small intermediate increase of 0.00002150 eV at evaluation 9. Forces need not decrease monotonically, and the actual ionic-convergence criterion is decisive.

The largest displacement between consecutive evaluated geometries is **0.011808 Å**. The minimum periodic interatomic distance throughout P4L is **2.29502 Å**. The outer-atom slab extent changes modestly from **45.88059 to 45.94657 Å**, with a fixed cell. These checks show finite, small ionic moves and no evidence of collapse, overlaps or a force/energy explosion along the sampled path. They do not establish global structural stability. No serious numerical warning is identified in the retained P4L output. CONTCAR agrees with the last force-evaluated XML geometry within 3.2×10⁻⁷ Å per Cartesian component.

**P4L passes structural relaxation under the selected recipe.** Its narrow force margin and the measured same-geometry sensitivity require the independent audit below before accepting the endpoint.

## 6. Independent endpoint force audit: F4L

The archived [F4L input](../../calculation/02_numerical_convergence/03_slab_thickness/00_setup_validation/F4L/INCAR) corresponds to the revised audit budget: `NELM=180`. Its POSCAR is byte-identical to P4L's final CONTCAR; KPOINTS and POTCAR are unchanged. Relative to P4L's INCAR, only `EDIFF: 10⁻⁶→10⁻⁷ eV`, `NELM: 120→180`, `NSW: 100→0` and `IBRION: 2→−1` change. The audit retains AMIN=0.01 and uses fresh `ISTART=0`, `ICHARG=2`, with no inherited electronic restart. The force-evaluated endpoint positions in the two XML records are identical at printed precision.

[F4L OUTCAR](../../calculation/02_numerical_convergence/03_slab_thickness/00_setup_validation/F4L/OUTCAR) explicitly reaches EDIFF at **iteration 139**. Final `dE=+8.0319×10⁻⁸ eV` and `d eps=−4.1843×10⁻⁹ eV` satisfy the tighter criterion. VASP elapsed time is **10,076.510 s (2.799 h)**. The 180-iteration limit is a ceiling, not a requirement to consume all iterations.

Differences below are F4L minus P4L on the same geometry. Forces use full atomic vector norms; component RMSE is over all 84 Cartesian components of the 28-atom slab.

| Endpoint quantity | P4L | F4L / comparison |
| --- | ---: | ---: |
| Electronic tolerance (eV) | 10⁻⁶ | 10⁻⁷ |
| Electronic iterations at endpoint | 6, propagated | 139, fresh |
| Maximum force (eV/Å) | **0.01818259** | **0.01850780** |
| Margin below 0.02 eV/Å | 0.00181741 | 0.00149220 |
| Free energy (eV/cell) | −114.45183967 | −114.45190037 |
| Free-energy difference (eV/cell) | — | **−0.00006070** |
| Energy-without-entropy difference (eV/cell) | — | −0.00005114 |
| Maximum atomic force-vector difference (eV/Å) | — | **0.00225570** |
| Force-component RMSE (eV/Å) | — | **0.00056516** |
| Cosine similarity of complete force arrays | — | 0.99237 |
| Maximum occupation difference, OUTCAR 0–2 scale | — | **0.04533** |

All audited force norms remain below 0.02 eV/Å. The stricter **0.001 eV/Å maximum force-vector reproducibility target is not met**: ten atoms exceed that difference, with the largest at atom 25. This target is a strong consistency test, not an automatic structural failure. The small component RMSE, closely aligned overall force field, small same-geometry energy difference and independently satisfied force threshold support the same structural conclusion. They do not imply mathematically identical electronic states.

Occupation comparison uses the same 52 k points and 148 energy-ordered bands per k point. The largest OUTCAR difference occurs at k point 1, band 124: **0.75111→0.79644**. Six states differ by more than 0.01; none differs by more than 0.05 on this scale. EIGENVAL reports occupations normalized to 0–1 for this non-spin-polarized calculation: its corresponding difference is **0.022665**, which becomes 0.045330 after multiplying by two. XML's more coarsely printed occupations give approximately 0.0452 after the same normalization. These are representation/precision differences, not contradictory electronic results. Nearby bands 122/123 at that k point are nearly degenerate, so band-index matching is descriptive rather than proof of orbital identity; no wavefunction-overlap matching or uniqueness claim is made.

The audit's approach to EDIFF is oscillatory, with intermittent amplification. In steps 131–139, median `|dE|` is 8.0761×10⁻⁶ eV and the largest `|dE|` is 2.1567×10⁻⁴ eV before the genuine final termination. The last printed charge residual is **0.00895 at iteration 138**, compared with **0.00110 at P4L endpoint iteration 5**. Thus tighter energy convergence does not establish uniformly smaller density residuals. Final F4L `rms=7.38×10⁻⁵`, and no BRMIX or non-Hermitian warning is identified. The audit passes the specified **electronic and endpoint-force criteria**, with finite residual sensitivity retained as a limitation; it is not evidence of a unique electronic ground state.

## 7. Validated production protocol

### Production relaxation

The accepted settings are those of [P4L INCAR](../../calculation/02_numerical_convergence/03_slab_thickness/00_setup_validation/P4L_r1/INCAR), checked against the executed OUTCAR/XML configuration. AMIN=0.01 is the accepted intervention; AMIX and BMIX remain their recorded default values. The validated physical model is non-spin-polarized PBE with D3(BJ), using PAW-PBE S, In_d and Zn datasets dated 06Sep2000 in the POSCAR species-block order. Species blocks must retain their corresponding potential ordering.

| Setting | Accepted value / interpretation |
| --- | --- |
| Electronic algorithm | `ALGO=Normal` (`IALGO=38` in executed output) |
| Mixing | **`AMIN=0.01`**; effective `AMIX=0.40`, `BMIX=1.0`, `IMIX=4`, `INIMIX=1` |
| Electronic tolerance | **`EDIFF=10⁻⁶ eV`** during relaxation |
| Force criterion | **`EDIFFG=−0.02 eV/Å`**; require every final atomic force norm below 0.02 eV/Å |
| Plane-wave and grid setup | `ENCUT=500 eV`, `PREC=Accurate`, `LREAL=False`, `LASPH=True` |
| Functional / dispersion | `GGA=PE`, `IVDW=12` |
| Occupations | `ISMEAR=0`, `SIGMA=0.05 eV` |
| Spin / symmetry | `ISPIN=1`, `ISYM=0` |
| Ionic optimizer / cell | `IBRION=2`, `ISIF=2`; effective `POTIM=0.5`; cell fixed |
| Electronic iteration controls | `NELM=120`, `NELMIN=2`; effective initial `NELMDL=−5` |
| Ionic budget | `NSW=100` ceiling in the validated full relaxation |
| Dipole correction | `LDIPOL=True`, `IDIPOL=3`, `DIPOL=(0.5,0.5,0.5)` in fractional coordinates |
| Sampling | Γ-centered **10×10×1**, zero shift; 52 sampled k points in the executed 4L calculation |
| Initial electronic state | Fresh `ISTART=0`, `ICHARG=2`; no inherited WAVECAR/CHGCAR/STOPCAR |
| Parallel resources | 8 MPI ranks × 1 thread; `NCORE=4`, `KPAR=1` |
| Output controls | `LCHARG=True`, `LWAVE=False`, `LVHAR=True`, `LORBIT=11` |
| Executable | VASP 6.6.1, standard complex executable recorded in native output |

`NELM=120` limits each electronic cycle; it does not require 120 iterations. Fresh initialization applies at launch; subsequent ionic evaluations propagate the electronic state within the same relaxation. Every force evaluation used for ionic motion must genuinely satisfy EDIFF. Exhausting the electronic budget is not an acceptable substitute. Likewise, reaching NSW or a walltime limit without ionic convergence leaves an incomplete relaxation rather than a successful structure.

These settings are accepted for the same basal slab construction and physical model. The center used for the dipole correction, potential ordering, mesh interpretation and fixed-cell geometry must remain consistent when preparing other thicknesses. Setup validation does not replace checking their actual electronic and structural outcomes.

### Independent endpoint audit

Use the final **force-evaluated** geometry with fresh electronic initialization and the same physical model, mixing, sampling and resources. The demonstrated audit uses **`NSW=0`, `IBRION=−1`, `EDIFF=10⁻⁷ eV`, `NELM=180`**. It requires genuine EDIFF termination, independent maximum force below 0.02 eV/Å, and quantified force/energy consistency. A maximum force-vector difference ≤0.001 eV/Å remains a strong reproducibility target; exceeding it requires assessment of the overall force field and force margin rather than automatic rejection. Occupation differences require consistent normalization and attention to band matching. The tighter audit tolerance is distinct from the validated 10⁻⁶ eV production relaxation tolerance.

## 8. Final decision and limitations

**4L setup validation: CLOSED / PASSED.** The AMIN=0.01 recipe genuinely converges the tested fixed 4L geometry, remains electronically convergent across ionic motion, produces a structurally converged full 4L relaxation, and independently reproduces the structural force criterion under tighter electronic convergence. Acceptance rests on that chain of evidence, not on one fortuitous EDIFF crossing or normal process termination.

Two qualifications remain: the endpoint does not meet the strongest force-vector reproducibility target, and the repository lacks the E0/E1/E2 outputs needed to substantiate the complete historical comparison. The latter limits retrospective claims about superiority over canonical mixing, but does not negate the directly verified later validation results. No additional calculation or automatic extension follows from this report.

The campaign does **not** establish 2L/4L/6L thickness convergence, that 4L is already sufficiently thick, surface-energy or work-function convergence, adsorption energetics, a unique electronic ground state, or transferability to unrelated facets or interfaces. Absolute total energies of different thicknesses are not compared here. The next scientific step is the formal **2L/4L/6L thickness-convergence study using this validated protocol**, with electronic, structural and property-convergence checks performed separately.

## Data provenance

All campaign conclusions and numerical results in this analysis were derived from the locally archived calculations under [calculation/02_numerical_convergence/03_slab_thickness/00_setup_validation](../../calculation/02_numerical_convergence/03_slab_thickness/00_setup_validation). Native VASP OUTCAR, OSZICAR, stdout, structures and electronic-state records were treated as primary evidence. XML was used for higher-precision forces and evaluated geometries, cross-checked against OUTCAR; force components agree within its 5×10⁻⁷ eV/Å printing precision. Local execution metadata was used only for resource, initialization and process checks. Its process-success fields were not substituted for scientific convergence evidence.

Window statistics, force differences and periodic geometry checks were recomputed from native records. A missing terminating-step `rms(c)` was left missing rather than replaced with zero. No parsed tables or plots were available in this directory, so none was linked or duplicated. E0/E1/E2 result omissions are explicitly retained as evidence gaps. Ignored local XML, execution metadata and licensed potential files support verification but are not assumed to be distributed with a Git checkout; linked evidence uses the archived repository paths. The external VASP references explain parameter semantics only and supply no campaign measurements.
