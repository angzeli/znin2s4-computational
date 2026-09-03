# WP1 ZnIn2S4 Geometry Optimisation — Final Audit

## Executive Summary

All five completed bulk geometry optimisations pass the final scientific audit. Each calculation terminated normally, satisfied the requested ionic force criterion, preserved its composition, and produced a crystallographically and chemically defensible final geometry. The five optimised structures are suitable to freeze as the geometry set for the next workflow stage.

The IIa calculation shows a clear IIa → IIa-prime reconstruction. Its final periodic cation–sulfur bond network consists of one S–Zn–Zn–S quadruple layer and two S–In–S–In–S quintuple layers, rather than two conventional S–Zn–S–In–S–In–S septuple layers. This conclusion is based on explicit periodic bonding and layer connectivity, not on the unchanged space-group label alone.

**Overall decision: 5 PASS, 0 FAIL.** No structure requires another geometry optimisation before the planned static-SCF stage. Final static-SCF calculations and k-point convergence checks remain necessary before using energy differences for a definitive polymorph ranking.

## Calculation Protocol

All five calculations used VASP `5.4.4.18Apr17-6-g9f103f2a35` (build 10 February 2022) with one common optimisation protocol:

- PBE (`GGA = PE`) with D3(BJ) dispersion (`IVDW = 12`)
- `ENCUT = 500 eV`, `PREC = Accurate`, `LREAL = .FALSE.`, and `LASPH = .TRUE.`
- `EDIFF = 1E-6 eV` and force-based `EDIFFG = -0.01 eV Å^-1`
- conjugate-gradient ionic optimisation (`IBRION = 2`) with full cell relaxation (`ISIF = 3`) and `NSW = 200`
- non-spin-polarised calculations (`ISPIN = 1`)
- Gaussian smearing (`ISMEAR = 0`, `SIGMA = 0.05 eV`)
- `LWAVE = .FALSE.` and `LCHARG = .FALSE.`

The PAW datasets were consistent across the complete set: `PAW_PBE S 06Sep2000`, `PAW_PBE Zn 06Sep2000`, and `PAW_PBE In_d 06Sep2000`. Their concatenation order matched the species order in the corresponding POSCAR in every calculation.

All KPOINTS files used Gamma-centred, unshifted meshes. The meshes were `4 × 4 × 4` for spinel, `12 × 12 × 2` for alpha1, `12 × 12 × 4` for beta, and `12 × 12 × 2` for IIa and IIb. The largest reciprocal-space interval was approximately 0.148 Å^-1 for spinel and 0.156–0.157 Å^-1 for the layered phases; alpha1 was deliberately denser along its long c direction. No meaningful protocol deviation was found.

## Final Structure Summary

| Structure | Converged | Ionic Steps | Max Force (eV Å^-1) | Final Energy (eV) | Energy / f.u. (eV) | Final Space Group | a (Å) | c (Å) | Volume / f.u. (Å³) | Lee a dev. | Lee c dev. | Verdict |
|---|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|---|
| Spinel | Yes | 9 | 0.000772 | -233.739012 | -29.217377 | Fd-3m (No. 227) | 10.61750 | 10.61750 | 149.6155 | -0.0236% | n/a | PASS |
| Alpha1 | Yes | 25 | 0.005699 | -86.247502 | -28.749167 | R3m (No. 160) | 3.88340 | 36.49448 | 158.8768 | +0.0876% | -0.0151% | PASS |
| Beta | Yes | 26 | 0.007821 | -28.755010 | -28.755010 | P3m1 (No. 156) | 3.87828 | 12.14843 | 158.2451 | -0.0442% | -0.1772% | PASS |
| IIa → IIa-prime | Yes | 81 | 0.009849 | -57.449033 | -28.724516 | P-3m1 (No. 164) | 3.88965 | 24.56060 | 160.9016 | -0.0090% | +0.0024% | PASS |
| IIb | Yes | 25 | 0.006224 | -57.502438 | -28.751219 | P6_3mc (No. 186) | 3.88335 | 24.35177 | 159.0170 | +0.0863% | -0.0748% | PASS |

The reported energies are final relaxation `TOTEN` values. They are useful for provenance and gross consistency checks, but they are **not** an authoritative polymorph-energy ranking because the cells and k-point meshes differ. A consistently converged static-SCF energy set is required for that purpose.

## Structure-by-Structure Audit

### Spinel

- **Identity:** S32Zn8In16, 56 atoms, corresponding to 8 ZnIn2S4 formula units; reduced composition ZnIn2S4.
- **Convergence:** normal VASP timing footer, explicit “reached required accuracy” termination, 9 ionic steps, and 8 electronic iterations in the final ionic step. The final maximum force is 0.000772 eV Å^-1 (RMS 0.000584 eV Å^-1), well below 0.01 eV Å^-1. The complete `vasprun.xml` contains all 9 ionic calculations and agrees with the OUTCAR final energy.
- **Final cell:** a = b = c = 10.617497 Å; α = β = γ = 90°; V = 1196.9238 Å³; V/f.u. = 149.6155 Å³.
- **Start-to-final change:** Δa = Δb = Δc = +0.2597%; ΔV = +0.7810%.
- **Symmetry:** Fd-3m, No. 227, for both starting and final structures; unchanged from `symprec = 10^-5` to `5 × 10^-2 Å`.
- **Chemistry and topology:** all 8 Zn are tetrahedral ZnS4 with Zn–S = 2.3682 Å; all 16 In are octahedral InS6 with In–S = 2.6149 Å. The full 56-atom periodic framework remains connected and no duplicate sites or implausible contacts occur. The shortest contact is 2.3682 Å (Zn–S).
- **Benchmark:** Lee et al. PBE+D3(BJ) a = 10.62 Å; deviation = -0.0236%.
- **Decision:** PASS.

### Alpha1

- **Identity:** Zn3In6S12, 21 atoms, corresponding to 3 ZnIn2S4 formula units; reduced composition ZnIn2S4.
- **Convergence:** normal termination, explicit ionic-accuracy marker, 25 ionic steps, and 6 electronic iterations in the final step. The final maximum force is 0.005699 eV Å^-1 (RMS 0.003761 eV Å^-1). The complete XML record agrees with OUTCAR.
- **Final cell:** a = b = 3.883399 Å, c = 36.494483 Å; α = β = 90°, γ = 120°; V = 476.6305 Å³; V/f.u. = 158.8768 Å³.
- **Start-to-final change:** Δa = Δb = +0.2737%; Δc = -1.5430%; ΔV = -1.0033%.
- **Symmetry:** R3m, No. 160, for both starting and final structures across `symprec = 10^-5` to `5 × 10^-2 Å`.
- **Chemistry and topology:** the final cell retains three disconnected periodic 7-atom Zn1In2S4 septuple networks. Zn is four-coordinate (2.3079–2.5380 Å); the two In environments are four-coordinate (2.4530–2.4849 Å) and six-coordinate (2.5730–2.7112 Å). The shortest contact is 2.3079 Å. Three equivalent outer-S van der Waals gaps of 2.7686 Å separate the septuple layers; no cation–S bond crosses these gaps.
- **Benchmark:** Lee et al. a = 3.88 Å and c = 36.50 Å; deviations = +0.0876% and -0.0151%, respectively.
- **Decision:** PASS.

### Beta

- **Identity:** Zn1In2S4, 7 atoms and 1 formula unit; reduced composition ZnIn2S4.
- **Convergence:** normal termination, explicit ionic-accuracy marker, 26 ionic steps, and 5 electronic iterations in the final step. The final maximum force is 0.007821 eV Å^-1 (RMS 0.004596 eV Å^-1). The complete XML record agrees with OUTCAR.
- **Final cell:** a = b = 3.878284 Å, c = 12.148433 Å; α = β = 90°, γ = 120°; V = V/f.u. = 158.2451 Å³.
- **Start-to-final change:** Δa = Δb = +0.7347%; Δc = -1.5524%; ΔV = -0.1006%.
- **Symmetry:** P3m1, No. 156, for both starting and final structures across `symprec = 10^-5` to `5 × 10^-2 Å`.
- **Chemistry and topology:** one complete periodic 7-atom Zn1In2S4 septuple network is retained. Zn is four-coordinate (2.3067–2.5373 Å); In occurs in four-coordinate (2.4510–2.4889 Å) and six-coordinate (2.5719–2.7092 Å) environments. The shortest contact is 2.3067 Å. The outer-S van der Waals gap is 2.7286 Å, with no cation–S bond across it.
- **Benchmark:** Lee et al. a = 3.88 Å and c = 12.17 Å; deviations = -0.0442% and -0.1772%, respectively.
- **Decision:** PASS.

### IIa → IIa-prime

- **Identity:** Zn2In4S8, 14 atoms, corresponding to 2 ZnIn2S4 formula units; reduced composition ZnIn2S4.
- **Convergence:** normal termination, explicit ionic-accuracy marker, 81 ionic steps, and 6 electronic iterations in the final step. The final maximum force is 0.009849 eV Å^-1 (RMS 0.005739 eV Å^-1), below the requested 0.01 eV Å^-1 threshold. The complete XML record agrees with OUTCAR.
- **Final cell:** a = b = 3.889649 Å, c = 24.560599 Å; α = β = 90°, γ = 120°; V = 321.8032 Å³; V/f.u. = 160.9016 Å³.
- **Start-to-final change:** Δa = Δb = +1.0298%; Δc = -0.4838%; ΔV = +1.5765%.
- **Symmetry:** P-3m1, No. 164, for both starting and final structures across `symprec = 10^-5` to `5 × 10^-2 Å`. The unchanged label does not conceal the bond-topology change described below.
- **Starting topology:** at a periodic cation–S cutoff of 2.90 Å, the input consists of two disconnected 7-atom `{Zn1 In2 S4}` septuple networks. Between the two van der Waals gaps, their plane order is the conventional S–Zn–S–In–S–In–S sequence, or its c-axis reverse.
- **Final topology:** the same periodic analysis gives one 4-atom `{Zn2 S2}` network and two 5-atom `{In2 S3}` networks. With the unit-cell origin treated cyclically, the final plane sequence is

  `S–In–S–In–S | S–Zn–Zn–S | S–In–S–In–S`,

  where the bars are the three outer-S van der Waals-gap boundaries. This is exactly two In–S quintuple layers plus one Zn–S quadruple layer, i.e. the IIa-prime topology.
- **Explicit Zn bond exchange:** Zn1 loses its starting apical S5 bond (2.4557 Å) and forms an S7 bond (2.5729 Å), while its three S8 bonds contract from 2.3548 to 2.3061 Å. Zn2 analogously loses S6 (2.4557 Å), forms S8 (2.5729 Å), and contracts its three S7 bonds from 2.3548 to 2.3061 Å. Thus the topology assignment is supported by actual bond exchange, not merely by a visual layer description.
- **Chemistry and contacts:** both Zn remain four-coordinate (2.3061–2.5729 Å). The two In environments remain four-coordinate (2.3934–2.5214 Å) and six-coordinate (2.5487–2.7643 Å). The shortest contact is 2.3061 Å. Final outer-S gap widths are 2.9063, 2.9063, and 2.9365 Å, and no cation–S bond crosses a disconnected-layer gap.
- **Classification:** **A. CLEAR IIa → IIa-prime RESTRUCTURING.**
- **Benchmark:** Lee et al. IIa-prime a = 3.89 Å and c = 24.56 Å; deviations = -0.0090% and +0.0024%, respectively.
- **Decision:** PASS.

### IIb

- **Identity:** Zn2In4S8, 14 atoms, corresponding to 2 ZnIn2S4 formula units; reduced composition ZnIn2S4.
- **Convergence:** normal termination, explicit ionic-accuracy marker, 25 ionic steps, and 7 electronic iterations in the final step. The final maximum force is 0.006224 eV Å^-1 (RMS 0.003619 eV Å^-1). The complete XML record agrees with OUTCAR.
- **Final cell:** a = b = 3.883347 Å, c = 24.351774 Å; α = β = 90°, γ = 120°; V = 318.0340 Å³; V/f.u. = 159.0170 Å³.
- **Start-to-final change:** Δa = Δb = +0.8661%; Δc = -1.3299%; ΔV = +0.3867%.
- **Symmetry:** P6_3mc, No. 186, for both starting and final structures across `symprec = 10^-5` to `5 × 10^-2 Å`.
- **Chemistry and topology:** two disconnected periodic 7-atom Zn1In2S4 septuple networks are retained. Zn is four-coordinate (2.3074–2.5413 Å); In occurs in four-coordinate (2.4511–2.4853 Å) and six-coordinate (2.5724–2.7117 Å) environments. The shortest contact is 2.3074 Å. Two equivalent outer-S van der Waals gaps of 2.7792 Å separate the layers, with no cation–S bond across either gap.
- **Benchmark:** Lee et al. a = 3.88 Å and c = 24.37 Å; deviations = +0.0863% and -0.0748%, respectively.
- **Decision:** PASS.

## Input Consistency

The INCAR settings, VASP executable version, PAW dataset identities, and element-to-POTCAR ordering are consistent across the five calculations. Every calculation used the same electronic and ionic convergence settings and the same dispersion treatment. The chosen meshes also provide comparable maximum reciprocal-space intervals:

- spinel `4 × 4 × 4`: approximately 0.148, 0.148, 0.148 Å^-1;
- alpha1 `12 × 12 × 2`: approximately 0.156, 0.156, 0.085 Å^-1;
- beta `12 × 12 × 4`: approximately 0.157, 0.157, 0.127 Å^-1;
- IIa and IIb `12 × 12 × 2`: approximately 0.157, 0.157, 0.127 Å^-1.

No input mismatch or isolated protocol deviation invalidates a cross-structure geometry comparison. These checks establish protocol consistency, but not fully converged energy differences; a dedicated k-point convergence test is still required for energetic ranking.

## Structural Sanity Checks

All final CONTCAR files parse completely, retain the expected integer composition and atom count, and contain neither partial occupancies nor duplicate periodic sites. The shortest final interatomic contacts lie between 2.3061 and 2.3682 Å and are chemically assigned cation–S bonds, not overlaps. With a 2.90 Å periodic cation–S cutoff, all Zn sites are four-coordinate and all In sites have the expected four- or six-coordinate environments for their crystallographic layers; spinel instead has the expected tetrahedral Zn and octahedral In framework.

For alpha1, beta, and IIb, the number and composition of the disconnected septuple-layer networks are preserved from input to output. Their outer-S plane separations remain meaningful van der Waals gaps, and no cation–S bond bridges a gap. IIa is the sole topology-changing case, and its reconstructed quadruple-plus-quintuple network is internally coordinated without unphysical contacts.

Symmetry was independently redetected from both starting and final Cartesian structures using spglib 2.7.0. The assigned groups were stable at `symprec` values of `10^-5`, `10^-4`, `10^-3`, `10^-2`, and `5 × 10^-2 Å`; no symmetry standardisation or atomic displacement was applied during this audit.

## Lee et al. Benchmark Comparison

The final lattice parameters agree closely with the PBE+D3(BJ) values reported by J. Lee et al., *Chemistry of Materials* **2019**, 31, 9148–9155, DOI `10.1021/acs.chemmater.9b03539`. Across all values for which a comparison is available, the largest absolute deviation is 0.1772% (the beta c parameter). Spinel, alpha1, beta, the reconstructed IIa-prime structure, and IIb are therefore all consistent with the corresponding published relaxed-cell benchmarks. The IIa-prime final cell is especially close to the reported a = 3.89 Å and c = 24.56 Å values.

## Acceptance Decision

- **Spinel:** accepted; freeze the converged geometry.
- **Alpha1:** accepted; freeze the converged geometry.
- **Beta:** accepted; freeze the converged geometry.
- **IIa → IIa-prime:** accepted as a clear IIa-prime reconstruction; freeze the converged geometry under the IIa-prime identity.
- **IIb:** accepted; freeze the converged geometry.

No geometry requires follow-up relaxation or manual structural repair based on the available calculation records. Overall acceptance is **5/5 PASS**.

## Recommended Next Step

1. Freeze the accepted final geometries without further ionic relaxation.
2. Perform and document a k-point convergence study suitable for energy differences between polymorphs.
3. Run consistently configured final static-SCF calculations on the frozen cells.
4. Use only those converged static-SCF energies, with identical energy-extraction conventions, for the final polymorph ranking.
