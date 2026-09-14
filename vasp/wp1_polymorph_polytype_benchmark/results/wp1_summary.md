# WP1 — ZnIn2S4 Polymorph / Polytype Benchmark

**Status: COMPLETE for the agreed FYP scope, 14 September 2026.** The five-phase PBE/XRD benchmark, targeted β HSE06 gap/near-edge validation and phase-selection handoff are concluded. β Stage 05-02 terminated normally at 14:13:53 China Standard Time (06:13:53 UTC). Its gap and near-edge results are accepted with the explicit high-unoccupied-band and finite-sampling caveats below; this is not an all-polymorph hybrid or complete optical-spectrum validation.

## 1. Executive Summary

WP1 established a consistent bulk benchmark for spinel, α₁ (`alpha1`), β (`beta`), IIa′ (`IIa_prime`) and IIb ZnIn2S4. All five relaxed structures passed geometry validation, and the subsequent static calculations supplied an accepted PBE+D3(BJ) energy ranking. The layered order is β < IIb < α₁ < IIa′, spanning only 29.594 meV per formula unit. IIb lies 3.513 meV/f.u. above β, making it a useful stacking-dependent comparison rather than a clearly separated energetic alternative.

Spinel is structurally and electronically distinct. Its calculated energy is 462.789 meV/f.u. below β, but this does not establish it as the ambient-pressure thin-film ground state. The experimental high-pressure association of spinel remains an unresolved interpretation boundary. The result is retained, not discarded or translated into an unsupported experimental phase assignment.

The completed PBE electronic comparison finds small layered gaps of approximately 0.279–0.304 eV and a larger spinel gap of 1.318 eV. Sulfur-p states dominate every valence edge; the layered conduction edges share mixed S-p/In-s/S-s character. Ideal powder XRD separates the overall spinel fingerprint from the layered patterns, while exposing substantial overlap among layered basal reflections and useful differences elsewhere.

The completed β HSE06 calculation gives a sampled gap of approximately **1.091 eV**, versus **0.303 eV at PBE**, an increase of approximately **0.788 eV** on the same accepted relaxed geometry. Both methods retain a Γ conduction minimum and an almost flat Γ–A valence edge, with the sampled valence maximum at or near A. The minimum direct gap exceeds the fundamental gap by only about 1.2–1.4 meV, so β remains **near-degenerate rather than unambiguously indirect**. HSE06 preserves this qualitative near-edge picture but does not act as an exact rigid shift across the full path.

WP1 is **complete for the agreed five-phase PBE/XRD benchmark plus targeted β HSE06 validation**. β remains the primary layered reference and WP2 parent. No WP1 calculation remains active or mandatory for this closure. Broader hybrid coverage and quantitative high-energy empty-band analysis are explicitly deferred; neither is silently claimed as completed. WP2 retains its PBE+D3(BJ) structural-screening workflow, with targeted electronic validation only where a later surface or reaction conclusion requires it.

## 2. Scope and Computational Protocol

The workflow separates geometry, energy ranking, electronic sampling and higher-level validation. Stage 01 used PBE with D3(BJ) dispersion, PAW-PBE Zn/In_d/S datasets, a 500 eV plane-wave cutoff and full cell/position relaxation. Electronic and force criteria were 10⁻⁶ eV and 0.01 eV Å⁻¹, respectively. All calculations were non-spin-polarised. Dataset ordering followed each structure's species ordering; no PAW contents are reproduced here.

The [k-point audit](02_kpoint_convergence/KPOINT_CONVERGENCE_AUDIT.md) directly tested spinel and β. Their accepted meshes differed from the densest tested meshes by 0.011139 and 0.118590 meV/f.u., respectively, comfortably below the 1 meV/f.u. criterion. Transfer to other layered cells used the actual reciprocal-vector intervals, including 2π, with a maximum spacing of approximately 0.156 Å⁻¹.

| Stage | Spinel sampling | Layered sampling | Role |
| --- | --- | --- | --- |
| 01 geometry | 4 × 4 × 4 | β: 12 × 12 × 4; others: 12 × 12 × 2 | Relaxed bulk structures |
| 03 static SCF, actually completed | 4 × 4 × 4 | All four: 12 × 12 × 4 | Authoritative energies and charge densities |
| 04 DOS/PDOS | 6 × 6 × 6 | All four: 18 × 18 × 6 | Dense fixed-density electronic sampling |
| 04 bands | cF2 path | α₁: hR1; others: hP2 | 40 points per segment |
| 05 β HSE validation | Not propagated | β: 12 × 12 × 4; optional hP2 path | Regular SCF, then 20-point/segment path |

Stage 02 recommended 12 × 12 × 2 for α₁, IIa′ and IIb, but the completed Stage 03 calculations used 12 × 12 × 4. The [static audit](03_static_scf/STATIC_SCF_ANALYSIS.md) explicitly accepts this as conservative c-axis oversampling. The summary records the executed protocol, not merely the earlier recommendation.

Stage 03 used fixed geometries, tetrahedron occupations and EDIFF = 10⁻⁷ eV; its final converged `free energy TOTEN`, divided by formula-unit count, defines the energy comparison. Stage 04 instead used `ICHARG=11` with the accepted Stage 03 charge density and unchanged structure. These were fixed-density evaluations, not new dense-mesh self-consistent charge calculations. D3(BJ) shaped the geometry; it is not a band-gap correction.

Version provenance is stage-specific: Stages 01–04 and the five Stage 05-00 PBE bootstrap donors used VASP `5.4.4.18Apr17-6-g9f103f2a35`, built in February 2022. All ten current Stage 04 OUTCAR headers confirm that build. Stage 05 local production used the validated native arm64 VASP 6.6.1 build dated 7 September 2026, launched using the validated synthetic OpenMPI topology. Native describes the executable architecture, not `--mpi-mode native`; production used `--mpi-mode synthetic`. A single executable did not produce the entire dataset.

## 3. Structural Benchmark

The authoritative structures are Stage 01 final CONTCARs, retained numerically in Stage 03. The [geometry audit](01_geometry_optimisation/FINAL_GEOMETRY_OPT_AUDIT.md) documents normal termination, convergence, composition, contacts and periodic bonding for each. The XRD [structure summary](simulated_xrd/wp1_simulated_xrd_summary.csv) independently records the same cells; densities below come from the [electronic summary](04_electronic_structure/electronic_structure_summary.csv).

| Phase | Space group | Atoms in cell / f.u. per cell | a (Å) | c (Å) | V/f.u. (Å³) | Density (g cm⁻³) | Maximum force (eV Å⁻¹) |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Spinel | Fd-3m, 227 | 56 / 8 | 10.61750 | 10.61750 | 149.6155 | 4.6976 | 0.000772 |
| α₁ | R3m, 160 | 21 / 3 | 3.88340 | 36.49448 | 158.8768 | 4.4238 | 0.005699 |
| β | P3m1, 156 | 7 / 1 | 3.87828 | 12.14843 | 158.2451 | 4.4414 | 0.007821 |
| IIa′ | P-3m1, 164 | 14 / 2 | 3.88965 | 24.56060 | 160.9016 | 4.3681 | 0.009849 |
| IIb | P6₃mc, 186 | 14 / 2 | 3.88335 | 24.35177 | 159.0170 | 4.4199 | 0.006224 |

Spinel is cubic; the other cells have hexagonal metrics. These are detected symmetries of the accepted calculated structures, not newly refined experimental assignments. The geometry audit found the labels stable over 10⁻⁵–0.05 Å symmetry tolerances without moving atoms. Its comparison against the published PBE+D3(BJ) reference cells gives a largest lattice-parameter deviation of 0.1772%; this is computational-reference agreement, not a new experimental fit.

β, IIb and α₁ contain one, two and three complete ZnIn2S4 septuple-layer networks per calculation cell. They retain tetrahedral Zn and tetrahedral/octahedral In environments, with no cation–S bond crossing their identified van der Waals gaps. Spinel instead retains a connected three-dimensional framework of tetrahedral Zn and octahedral In.

The legacy IIa starting model relaxed to **IIa′**, not simply a slightly distorted IIa cell. Periodic connectivity changed from two septuple networks to one S–Zn–Zn–S quadruple network and two S–In–S–In–S quintuple networks. Explicit Zn–S bond exchange supports that identification despite the unchanged P-3m1 label. The accepted source therefore remains [IIa_to_IIa_prime/CONTCAR](../calculation/01_geometry_optimisation/IIa_to_IIa_prime/CONTCAR), consistently described as IIa′ downstream.

## 4. Relative Energetics and Stability

The [Stage 03 energy CSV](03_static_scf/static_scf_energies.csv), not the relaxation energies or bootstrap energies, is authoritative. Every static run converged, terminated normally and agreed between OUTCAR and XML. All cells have the same composition per formula unit, allowing direct comparison under the common Hamiltonian.

| Phase | Energy (eV/f.u.) | ΔE versus spinel (meV/f.u.) | ΔE versus β (meV/f.u.) |
| --- | ---: | ---: | ---: |
| Spinel | −29.217802 | 0.000 | −462.789 |
| α₁ | −28.748661 | 469.141 | 6.352 |
| β | −28.755013 | 462.789 | 0.000 |
| IIa′ | −28.725418 | 492.384 | 29.594 |
| IIb | −28.751500 | 466.302 | 3.513 |

β is the lowest-energy layered reference. The 2.839 meV/f.u. α₁–IIb separation is small, and even IIa′ remains within approximately 30 meV/f.u. of β. These results motivate a selective comparison of stacking and layer topology; they do not determine equilibrium film fractions or exclude metastable products.

Spinel's substantially lower calculated energy survived the protocol, composition, normalisation, k-point and residual-stress checks. The accepted analysis also records comparable spinel stabilisation in its literature comparison. Nevertheless, a zero-temperature internal-energy ordering is not a pressure-dependent or finite-temperature phase diagram. Vibrational free energies, pressure-dependent enthalpies, kinetics and synthesis conditions were not resolved. The tension between the low computed energy and the experimental high-pressure association is scientifically interesting but outside the minimum FYP closure scope.

No standalone relative-energy diagram is present in the inspected result set. The validated numerical ranking is complete; an additional diagram would be a presentation asset, not missing scientific evidence or a closure blocker.

## 5. Simulated XRD

The frozen [powder-XRD comparison](simulated_xrd/wp1_simulated_xrd_comparison.pdf) uses the exact accepted relaxed cells without standardisation, primitive reduction or symmetrisation before diffraction. The reproducible [script](../script/simulate_wp1_xrd.py) uses pymatgen's CuKa convention, wavelength 1.54184 Å, over 5–80° in 2θ. Cu Kα is a configurable simulation assumption, not a confirmed experimental instrument specification. Each phase is normalised to a maximum of 100 within that interval; the figure contains unbroadened sticks with display-only offsets.

Numerical records are the [summary](simulated_xrd/wp1_simulated_xrd_summary.csv) and five reflection CSVs: [spinel](simulated_xrd/spinel_xrd_reflections.csv), [α₁](simulated_xrd/alpha1_xrd_reflections.csv), [β](simulated_xrd/beta_xrd_reflections.csv), [IIa′](simulated_xrd/IIa_prime_xrd_reflections.csv) and [IIb](simulated_xrd/IIb_xrd_reflections.csv). They retain d-spacings, intensities, contributing families and multiplicities. Hexagonal CSVs retain four-index families; the figure and following table omit only the redundant third index for compact three-index labels.

| Phase | Three strongest peaks: 2θ in degrees (hkl; relative intensity) |
| --- | --- |
| Spinel | 27.869 (311; 100.0), 48.502 (440; 71.1), 33.768 (400; 47.4) |
| α₁ | 28.282 (104; 100.0), 21.919 (009; 97.8), 46.786 ($2\bar{1}0$; 82.1) |
| β | 27.554 (101; 100.0), 21.949 (003; 52.2), 30.411 (102; 52.1) |
| IIa′ | 27.456 (102; 100.0), 46.706 ($2\bar{1}0$; 52.3), 36.994 (107; 45.6) |
| IIb | 21.899 (006; 100.0), 28.733 (103; 95.1), 46.786 ($2\bar{1}0$; 83.7) |

Spinel is distinguishable through the combined cubic fingerprint, not through a claim that every individual peak is unique. Among α₁, β and IIb, the strong basal peaks near 21.9° differ by only about 0.05°, and their low-angle basal peaks also nearly coincide. Those reflections alone are poor stacking identifiers. The in-plane $(2\bar{1}0)$ peaks likewise cluster near 46.8°.

The full layered patterns are not identical. The 26–31° region contains different peak sequences: α₁ has its maximum at 28.282°, β at 27.554°, and IIb has strong peaks at 27.515° and 28.733°. Some contrasts depend strongly on intensity, including β's 30.411° peak versus the weaker nearby IIb peak at 30.363°. IIa′ offers a more distinctive 18.059° (005) peak with relative intensity 41.1; none of the other frozen patterns contains a nearby reflection. This is a candidate discriminator within this five-phase set, not a universal exclusion of substrates or impurities.

Thus ideal multi-peak powder fingerprints can differentiate the models, but a readily resolved experimental distinction among the closely related layered stackings is not established. Random orientation underlies the simulated intensities. Texture, strain, finite crystallite size, mixed phases and instrumental effects can alter thin-film patterns; strong experimental (00l) intensity must not be converted directly into phase abundance. Positions and reflection combinations are the primary comparison targets, with possible DFT-lattice, temperature and substrate-induced shifts retained as caveats.

## 6. PBE Electronic Structure

The [accepted Stage 04 analysis](04_electronic_structure/ELECTRONIC_STRUCTURE_ANALYSIS.md), [electronic summary](04_electronic_structure/electronic_structure_summary.csv) and [edge-coordinate table](04_electronic_structure/band_edges.csv) support the following finite-mesh assignments. The fundamental gap is the sampled minimum conduction energy minus the sampled maximum valence energy; the direct gap minimises that separation at the same sampled k point.

| Phase | Fundamental gap (eV) | Minimum direct gap (eV) | Direct minus fundamental (meV) | Sampled classification | Sampled VBM | Sampled CBM |
| --- | ---: | ---: | ---: | --- | --- | --- |
| Spinel | 1.318 | 1.428 | 110.2 | Indirect | (0.5, 0.5, 0) | Γ |
| α₁ | 0.304 | 0.304 | 0.0 | Direct | Γ and (0, 0, z), z = 1/6, 1/3, 1/2 | Γ |
| β | 0.303 | 0.304 | 1.4 | Ambiguous / near-degenerate | A = (0, 0, 0.5) | Γ |
| IIa′ | 0.304 | 0.346 | 42.0 | Indirect | (0.111111, 0.055556, 0) | Γ |
| IIb | 0.279 | 0.279 | 0.0 | Direct | Γ | Γ |

Coordinates refer to the calculation-cell reciprocal basis and symmetry-reduced representatives. The β separation is only 1.380 meV before rounding, so calling it unambiguously indirect would exceed the accepted evidence. IIa′ retains an off-Γ valence maximum and a larger 41.959 meV separation. IIb is sampled direct, while α₁ includes Γ among several degenerate valence representatives.

The [PDOS edge table](04_electronic_structure/pdos_band_edge_character.csv) gives the complementary orbital comparison. Percentages below use the existing 0.50 eV edge windows, normalised within the PAW-projected weight, not to exact chemical populations.

| Phase | Leading VBM projection | Main CBM projections |
| --- | --- | --- |
| Spinel | S-p 89.1% | In-s 39.6%, S-s 21.9%, S-p 17.5%, Zn-s 15.9% |
| α₁ | S-p 76.1% | S-p 42.5%, In-s 32.7%, S-s 19.1% |
| β | S-p 76.9% | S-p 43.9%, In-s 31.7%, S-s 18.4% |
| IIa′ | S-p 78.8% | S-p 44.1%, In-s 30.0%, S-s 19.8% |
| IIb | S-p 75.8% | S-p 42.8%, In-s 32.5%, S-s 18.9% |

S-p dominance at the valence edge is common to all five phases. Layered conduction edges are hybridised rather than pure orbital states; the leading channel can depend on the integration window. Spinel has a distinct In-s/S-s-rich conduction edge and appreciable Zn-s contribution. The narrower 0.20 eV analysis preserves these qualitative distinctions. Similar orbital identities do not mean stacking has no electronic effect or establish carrier-separation dynamics.

Canonical band PDFs: [spinel](04_electronic_structure/figures/spinel_pbe_band_structure.pdf), [α₁](04_electronic_structure/figures/alpha1_pbe_band_structure.pdf), [β](04_electronic_structure/figures/beta_pbe_band_structure.pdf), [IIa′](04_electronic_structure/figures/IIa_prime_pbe_band_structure.pdf), [IIb](04_electronic_structure/figures/IIb_pbe_band_structure.pdf). Canonical near-edge selected-orbital PDOS: [spinel](04_electronic_structure/figures/spinel_pbe_selected_orbital_pdos.pdf), [α₁](04_electronic_structure/figures/alpha1_pbe_selected_orbital_pdos.pdf), [β](04_electronic_structure/figures/beta_pbe_selected_orbital_pdos.pdf), [IIa′](04_electronic_structure/figures/IIa_prime_pbe_selected_orbital_pdos.pdf), [IIb](04_electronic_structure/figures/IIb_pbe_selected_orbital_pdos.pdf). These frozen assets were inspected, not regenerated.

## 7. Folded-Cell and Electronic-Structure Caveats

Spinel and α₁ band plots retain calculation cells containing four and three primitive cells, respectively. Their bands are folded; the canonical primitive-path labels are represented in the original cell basis, not an unfolded primitive-Brillouin-zone spectral representation. The [path justification](04_electronic_structure/KPOINTS_SELECTION_JUSTIFICATION.md) documents the coordinate transformation and charge-density compatibility.

Finite meshes and high-symmetry lines do not prove continuous-zone global extrema. The path samples give slightly smaller gaps than the uniform meshes for spinel and IIa′, by 3.440 and 4.580 meV, respectively; both records are retained. β remains near-degenerate at PBE. Absolute eigenvalues from different cells are not vacuum-aligned band offsets, and the PBE gaps are not quantitative experimental optical gaps. Dispersion ranges are not effective-mass calculations.

## 8. Stage 05: Targeted HSE06 Validation

### 8.1 Completed β Stage 05-01

The [Stage 05-00 bootstrap audit](05_hybrid_validation/05_00_PBE_WAVECAR_BOOTSTRAP_AUDIT.md) accepts all five WAVECARs as PBE orbital donors for Stage 05-01. Later symmetry/restart-provenance work resolved the current Stage 05-01/05-02 restart route. The production sequence is PBE bootstrap → regular-mesh HSE06 SCF → converged HSE WAVECAR → β HSE KPOINTS_OPT band validation. The [hybrid sampling justification](05_hybrid_validation/KPOINTS_SELECTION_JUSTIFICATION.md) preserves full crystal-symmetry/time-reversal orbit coverage and multiplicity checks: different equivalent irreducible representatives are not a sampling mismatch.

β Stage 05-01 is a converged HSE reference and donor. The archived [OUTCAR](../calculation/05_hybrid_validation/01_hse06_scf/beta/OUTCAR), [OSZICAR](../calculation/05_hybrid_validation/01_hse06_scf/beta/OSZICAR) and [stdout](../calculation/05_hybrid_validation/01_hse06_scf/beta/stdout.log), supplemented by local runtime metadata, establish convergence and normal termination independently of an orchestration status label.

| Evidence | Verified β 05-01 result |
| --- | --- |
| Method / executable | VASP 6.6.1; HSE06, AEXX = 0.25, HFSCREEN = 0.20 Å⁻¹ |
| Algorithm / symmetry | Davidson (`ALGO=Normal`), ACE enabled, ISYM = 3 |
| Regular sampling / bands | Γ-centred 12 × 12 × 4; NKPTS = 69, NBANDS = 40 |
| Restart | PBE WAVECAR accepted; ISTART = 1, ICHARG = 0 |
| Convergence | 18 electronic iterations; explicit EDIFF-reached marker |
| Final electronic energy change | +5.641476 × 10⁻⁸ eV; below EDIFF = 10⁻⁷ eV |
| Final total energy | −32.73892725 eV per seven-atom cell |
| VASP elapsed wall time | 125566.593 s = 34.880 h |
| Termination / restart output | Normal footer, complete XML, launcher exit 0; newly written 91679104-byte WAVECAR |

XML independently contains 18 electronic steps and the same final energy. The WAVECAR headers in the completed Scratch donor and repository copy contain 69 k points, 40 bands and a 500 eV cutoff. The later band-run staging record identifies this HSE donor, whose successful reading supplies additional restart-use evidence. The HSE energy is not inserted into the five-phase PBE+D3(BJ) ranking: it is a different electronic method on the accepted relaxed geometry, not a sixth comparable thermodynamic datum.

The generic launcher's `scientific_convergence: NOT ASSESSED` field is not a failed convergence verdict; raw VASP records establish convergence here. Conversely, older Stage 05 reports saying production HSE convergence had not yet been demonstrated describe their earlier checkpoint and are now superseded for β only. The β 05-01 results were archived and pushed in commit `5cffc22`; machine-specific metadata and ignored restart/XML files remain local.

### 8.2 Completed β Stage 05-02: Gap and Near-Edge Validation

**PASS WITH CAVEAT for the sampled gap and near-edge comparison.** Attempt `beta-20260911-101155` completed normally on 14 September 2026. The archived [OUTCAR](../calculation/05_hybrid_validation/02_hse06_band/beta/OUTCAR), [OSZICAR](../calculation/05_hybrid_validation/02_hse06_band/beta/OSZICAR) and [stdout](../calculation/05_hybrid_validation/02_hse06_band/beta/stdout.log) record the completed calculation; results were archived and pushed in commit `1590692`. Launcher and CMW payload exit codes were both zero. Runtime completion and the numerical acceptance tests below are separate evidence.

#### Protocol, provenance and convergence

The resolved method is VASP 6.6.1 HSE06: AEXX = 0.25, HFSCREEN = 0.20 Å⁻¹, 500 eV, ALGO = Normal, ACE enabled, ISYM = 3, ISPIN = 1 and no SOC. The accepted seven-atom Stage 01 β CONTCAR is numerically identical, within 10⁻¹⁰ in lattice-vector components and periodic fractional coordinates, to the Stage 04 references, HSE input, HSE donor geometry and final HSE CONTCAR. The HSE and Stage 04 PAW files are byte-identical. No geometry was relaxed or transformed for this comparison. The launcher added only NCORE = 1 and KPAR = 1 to the repository INCAR; it used eight MPI ranks.

The [regular KPOINTS](../calculation/05_hybrid_validation/02_hse06_band/beta/KPOINTS) specify Γ-centred 12 × 12 × 4 sampling, with 69 irreducible points. The converged 05-01 WAVECAR was read successfully with ISTART = 1 and ICHARG = 0, resolving 40 bands and 62 electrons without an ISTART = 0 fallback. HFRCUT = −1 was specified. The [KPOINTS_OPT](../calculation/05_hybrid_validation/02_hse06_band/beta/KPOINTS_OPT) hP2 path is Γ–M–K–Γ–A–L–H–A | L–M | H–K, with 20 points per segment: 180 entries in total. It matches the Stage 04 path and reciprocal basis, but Stage 04 has 40 points per segment. The optional-path evaluation follows regular-mesh self-consistency; it does not use the PBE fixed-density ICHARG = 11 shortcut. This distinction follows the [VASP hybrid-band workflow](https://vasp.at/wiki/Band-structure_calculation_using_hybrid_functionals).

| Completion evidence | Verified β 05-02 result |
| --- | --- |
| Regular SCF | Eight iterations; explicit EDIFF-reached marker; final energy −32.73892722 eV/cell |
| Optional batch 1–69 | 11 iterations; final absolute energy change 6.1798 × 10⁻⁸ eV |
| Optional batch 70–138 | 11 iterations; final absolute energy change 6.7375 × 10⁻⁸ eV |
| Optional batch 139–180 | 11 iterations; final absolute energy change 3.8024 × 10⁻⁸ eV |
| Stopping target | EDIFF = 10⁻⁷ eV; all three batch changes below target |
| VASP elapsed time | 273605.651 s = 76.002 h |
| Termination | Normal OUTCAR footer, closed/parseable XML, launcher/payload exit 0 |

The optional-loop energy is a band-eigenvalue convergence sum, not the cell total energy or a band gap. The archived MPI TCP-socket warning did not prevent subsequent convergence and normal termination. A NUL byte in stdout and an overflowed XML `totalsc` timing field are retained as formatting anomalies: clean OSZICAR iteration records and the OUTCAR elapsed time provide cross-checks. Neither anomaly enters the eigenvalue analysis.

#### Numerical source and validation boundary

The 180-point spectrum was read from `calculation/05_hybrid_validation/02_hse06_band/beta/vasprun.xml`, specifically `modeling/calculation/eigenvalues_kpoints_opt`. The ordinary [EIGENVAL](../calculation/05_hybrid_validation/02_hse06_band/beta/EIGENVAL) contains **69 regular-mesh points, not the optional band path**. The local XML is complete and byte-identical to the completed runtime copy. It remains ignored under the repository's existing data policy; the [portable eigenvalue CSV](05_hybrid_validation/beta_hse06_kpoints_opt_eigenvalues.csv) now preserves all **180 × 40 = 7200 optional-path eigenvalues**, including bands 38–40 unchanged despite their documented high-empty-state convergence limitation. The [provenance metadata](05_hybrid_validation/beta_hse06_kpoints_opt_metadata.json) records source hashes, geometry, sampling, the common path-VBM zero and repeated-point diagnostics. The [export/plot script](../script/export_wp1_hse06_band.py) can regenerate the [standalone β HSE06 band figure](05_hybrid_validation/figures/beta_hse06_band_structure.pdf) from the CSV alone using the frozen Stage 04 plotting primitives. Its accepted near-edge **−4 to +4 eV** presentation range excludes the problematic upper states from view but does not waive their convergence caveat. No PBE/HSE overlay is claimed.

Analysis used the existing dependency-free structure/EIGENVAL readers with explicit paths and a standard-library XML parser, trimming XML parameter-name whitespace. All 180 path entries and 40 eigenvalues per entry are present, finite and band-ordered. Coordinates reproduce the requested segment order to within 5 × 10⁻⁹ in fractional reciprocal coordinates. With NELECT = 62 and ISPIN = 1, bands 31 and 32 (one-based) define the valence and conduction edges; regular-mesh occupations confirm the occupied/empty assignment. The optional XML stores eigenvalues to 0.0001 eV, so extra decimal places are not inferred.

The regular-mesh eigenvalues agree independently between EIGENVAL, XML and OUTCAR to their printed precision (maximum difference 0.00005 eV). All 18 repeated-path-point pairs agree for bands 1–37 at the optional XML precision; the VBM/CBM bands also agree with coincident regular-mesh points within 0.000064 eV. These checks support the gap and near-edge comparison.

**The entire 40-band spectrum is not uniformly validated.** At repeated path points, the maximum discrepancies in bands 38, 39 and 40 are 0.0005, 0.0048 and 0.1421 eV, respectively; band 40 differs from its coincident regular-mesh value by up to 0.14412 eV. These bands lie at least approximately 5.50 eV above the path VBM. The differences are consistent with insufficient convergence of the highest empty states, which is not independently guaranteed by an occupied-band-energy stopping test. They do not affect the extracted band-31/32 gap, but prohibit treating the upper spectrum or optical transitions involving those states as quantitatively converged. Additional empty-state convergence work is a follow-up only if those observables become necessary.

#### HSE06 versus the frozen Stage 04 PBE reference

For each sampling set, the sampled fundamental gap is `min(E32) − max(E31)` and the minimum direct gap is `min(E32 − E31)` at the same k point. Uniform-mesh and path results are retained separately; neither proves continuous-zone global extrema.

| β dataset | Sampling | Sampled fundamental gap (eV) | Minimum direct gap (eV) | Direct minus fundamental (meV) |
| --- | --- | ---: | ---: | ---: |
| Stage 04 PBE uniform | 18 × 18 × 6; 202 irreducible points | 0.303118 | 0.304498 | 1.380 |
| Stage 04 PBE path | 40 points/segment; 360 entries | 0.303118 | 0.304498 | 1.380 |
| Stage 05-02 HSE06 uniform | 12 × 12 × 4; 69 irreducible points | 1.090677 | 1.091920 | 1.243 |
| Stage 05-02 HSE06 path | 20 points/segment; 180 entries | 1.0907 | 1.0919 | 1.2 |

The recalculated PBE numbers reproduce the frozen Stage 04 summary. Both uniform meshes place the VBM at A and CBM at Γ, with the minimum direct gap at Γ. At the optional XML precision, the HSE path VBM is shared by points near A along Γ–A (fractional kz approximately 0.447–0.500). The HSE path and its own uniform mesh agree in gap to better than 0.0001 eV. The original 05-01 donor's regular gap was 1.091399 eV, 0.000722 eV above the restarted 05-02 value; this small restart sensitivity reinforces reporting the HSE gap as approximately 1.091 eV, not as a sub-meV accuracy claim.

HSE06 increases the sampled β gap by approximately **0.788 eV** while retaining the **near-degenerate A/Γ edge arrangement**. Under the existing Stage 04 10-meV near-degeneracy convention, the 1.2–1.4 meV direct–fundamental separation does not support a categorical indirect-gap claim. Neither functional establishes a distinct, experimentally resolved change in directness.

| Band-edge range along a segment (eV) | PBE path | HSE06 path |
| --- | ---: | ---: |
| Top valence, Γ–M | 0.604 | 0.662 |
| Bottom conduction, Γ–M | 1.232 | 1.283 |
| Top valence, Γ–A | 0.00138 | 0.0012 |
| Bottom conduction, Γ–A | 0.674 | 0.609 |

The strong contrast between in-plane valence dispersion and an almost flat Γ–A valence edge survives HSE06. These are sampled band ranges, not effective masses or carrier mobilities, and the two paths have different point densities. An independent comparison at exactly shared special points also rules out an exact uniform scissor correction across the path: the HSE–PBE direct-gap increases are approximately 0.7874 eV at Γ, 0.7225 eV at A and 0.8604 eV at K. Thus the main change is gap opening with modest, non-rigid near-edge dispersion changes, not wholesale replacement of the PBE band-edge picture.

The comparison uses identical geometry, PAW datasets and cutoff, but different VASP versions, electronic-density protocols and sampling densities; it is not a perfectly isolated same-build functional-only benchmark. No HSE orbital-character validation is claimed: LORBIT = 0 and no HSE projected-band dataset was produced. The PBE PDOS assignments in Section 6 remain PBE-only. The approximately 1.091 eV HSE result is a calculated electronic gap, not a quantitative experimental optical-gap prediction, vacuum-referenced band alignment or evidence of catalytic activity.

### 8.3 Why the Remaining Polymorph HSE Matrix Was Not Completed

The final scope is an intentional resource-gated decision, not a failed calculation campaign. β is the lowest-energy layered reference and primary downstream surface parent; its 34.880 h SCF demonstrated that dense local HSE was feasible but expensive. IIb was then used as a bounded cost benchmark before releasing the remaining phases.

After **40480.214 s (11.2445 h)** of wall time, IIb had not completed its first electronic iteration in the surviving record. This alone exceeded the agreed approximately 3 h/iteration resource gate, so the calculation was deliberately cancelled. Stdout contains no completed DAV line and OUTCAR no completed LOOP timing; no exact IIb per-iteration duration is claimed. This is bounded resource evidence supporting an intentional scope decision, not a completed scientific HSE result.

The preserved attempt is `20260910-225541-2ffe65f1/IIb` under `/Volumes/Scratch/wp1_stage05_01_hse06_scf_runs/`. Its metadata reports accepted WAVECAR reading, SIGTERM, no normal footer and an unchanged staged donor WAVECAR. Read-only CMW history in state `wp1-IIb-05-01` records J1.1 cancellation requested and confirmed on 11 September at 02:10:32–33 UTC.

**IIb, IIa′, α₁ and spinel HSE production are deferred for publication-level follow-up.** IIb was stopped after sufficient cost evidence; the other three were not propagated into production. The internally consistent five-phase PBE comparison and targeted β HSE anchor offer greater immediate FYP value than filling the entire hybrid matrix at the expense of surface/reaction work. Later HPC execution remains appropriate if a publication question requires all-phase HSE gaps or projections. One β result does not prove a rigid correction for every phase.

## 9. Reconciliation with the Original ROADMAP Deliverables

The [original WP1 roadmap](../../ROADMAP.md) is retained as planning provenance. Its deliverables are reconciled here rather than silently rewritten.

| Original deliverable | Final WP1 interpretation and status |
| --- | --- |
| Optimised structural parameters | **COMPLETE:** five accepted structures and structural table above. |
| Relative-energy diagram | **SCIENTIFIC RESULT COMPLETE:** the validated Stage 03 numerical energy ranking is complete; a standalone relative-energy diagram was not generated and is a presentation asset rather than a WP1 closure requirement. |
| Simulated XRD comparison | **COMPLETE AND FROZEN:** canonical powder PDF and reflection CSVs. |
| Band-gap comparison | **COMPLETE:** five-phase PBE comparison plus accepted β HSE06 gap/near-edge validation with explicit upper-empty-band caveat; no all-phase HSE table is claimed. |
| VBM/CBM character comparison | **COMPLETE at PBE:** projected edge windows and canonical PDOS; not all-phase HSE projections. |
| WP2 shortlist | **HANDOFF DEFINED:** β primary, IIb provisional stacking control, spinel bounded feasibility contrast; WP2 chooses the final surfaces. |

The old generic “2–3 structures” instruction becomes a phase hierarchy, not an obligation to propagate three bulk phases. Routine multi-structure HSE, spin–orbit coupling (SOC) and electron/hole effective-mass calculations are no longer mandatory WP1 closure criteria. They remain targeted follow-up options when a specific mechanism or publication claim needs them. The roadmap also mentioned band-decomposed charge density; no completed WP1 charge-localisation deliverable is asserted here. Selected state-localisation analysis is part of the surface-focused WP2 plan, not something silently inferred from bulk PDOS.

## 10. WP2 Handoff

The current [WP2 roadmap](../../wp2_surface_and_facet_screening/ROADMAP.md) will select **reaction-facing surfaces**, not merely bulk phases. β supplies the primary basal and prismatic starting geometries. Its lowest layered energy and HSE donor justify priority, but do not establish the experimentally dominant film termination. IIb is a provisional stacking-dependent control because its energy is close and its stacking differs. It should expand only where a bounded comparison adds information.

Spinel offers a three-dimensional structural/electronic contrast, conditional on constructing chemically defensible slabs. Its low bulk energy is not sufficient experimental justification. α₁ and IIa′ remain meaningful bulk benchmarks but do not automatically enter the initial surface matrix.

WP2 [Stage 0](../../wp2_surface_and_facet_screening/results/00_surface_registry/SURFACE_GENERATION_AUDIT.md) subsequently deferred the tested unreconstructed spinel (001)/(110) cleavage-only branches pending a chemically justified compensation/reconstruction treatment; this does not reinterpret the WP1 bulk result.

WP2 inherits exact relaxed parents, the validated sampling logic, electronic hypotheses and comparison caveats. It must still establish termination identity, stoichiometry, polarity, thickness/vacuum convergence and surface-specific electronic behaviour. Its planned VASP 6.6.1 surface energies require matched parent-phase bulk static references, not blind subtraction of historical energies from another executable or the lowest polymorph. Bulk gaps and XRD fingerprints alone do not rank catalytic activity.

The completed β comparison does **not** trigger automatic reruns, HSE slab relaxations or expansion of the polymorph HSE matrix. WP2 retains PBE+D3(BJ) slab generation, relaxation and structural/numerical screening. Its consequential electronic conclusions may receive separately justified and budgeted higher-level checks. A bulk gap opening must not be applied automatically to other polytypes, surface/defect states, adsorption energies or reaction barriers. Relative VBM-zero band comparisons also do not partition the correction into absolute VBM and CBM shifts; vacuum/potential alignment remains a separate surface-specific task.

## 11. Remaining Limitations and Publication Follow-Up

No unresolved **FYP blocker** remains for the agreed five-phase PBE/XRD benchmark, targeted β HSE06 gap/near-edge validation or WP2 handoff. β 05-02 is finished and its accepted numerical conclusions are recorded above. The high-empty-band limitation is not waived: it limits the observables claimed, rather than requiring an unrequested full-spectrum calculation to close this scoped benchmark.

**Publication follow-up opportunities** include all-polymorph HSE comparisons, convergence of high unoccupied bands where required, an optional combined PBE/HSE comparison figure, question-driven SOC/effective-mass calculations, finer extremum sampling, appropriate folded-band interpretation, experimental textured-film diffraction comparison, and pressure/temperature-dependent interpretation of spinel. HSE mesh convergence is not guaranteed by PBE energy convergence, and 20 points per segment is a bounded sampling choice rather than a path-convergence study. These boundaries constrain later claims without invalidating the accepted gap/near-edge result.

Provenance limitations remain explicit: historical Stage 05 audits predate the completed β runs, and the surviving IIb logs do not resolve an exact iteration duration. The optional HSE path is now preserved in the portable CSV and standalone figure, while the ignored native XML remains the source-of-truth rather than the regular-mesh EIGENVAL. This summary uses raw completion evidence for β and cancellation/elapsed-time evidence for IIb without editing those historical records. Native XML exists both in the repository outside Scratch and in the byte-identical Scratch runtime copy; no independent backup has been verified. Verifying an independent native-source backup is the sole remaining archival follow-up. Local CMW metadata remains supporting local evidence, not a portable repository deliverable.

## 12. WP1 Closure Statement

WP1 established a validated five-phase structural dataset, relative PBE+D3(BJ) stability, ideal powder-XRD fingerprints, a five-phase PBE electronic comparison, PBE band-edge orbital trends, a completed targeted β HSE06 gap/near-edge benchmark and a scientifically bounded handoff to WP2. The β HSE result supports retaining the qualitative PBE near-edge picture while increasing the predicted electronic gap for this reference; it does not validate universal gap shifts, absolute redox alignment or catalytic performance.

**WP1 is complete for the agreed FYP benchmark and phase-selection scope as of 14 September 2026. No mandatory WP1 calculation remains.** Acceptance is explicitly limited to the demonstrated observables, with upper-empty-band convergence, all-phase HSE and the other publication follow-ups above left outside this closure. Stage 04 figures were frozen in commit `5a1dea3`; the XRD script and seven outputs in `45259cf`; β 05-01 and 05-02 results were archived in `5cffc22` and `1590692`. This archival pass adds the portable HSE export, provenance metadata, standalone figure and focused script/tests, and updates the present summary for human review. No calculation, input, native result, frozen Stage 04 plotting workflow or WP2 file was changed, and no staging, commit or push was performed.
