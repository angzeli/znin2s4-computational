# WP1 ZnIn2S4 K-Point Convergence — Final Audit

## Executive Summary

- **Recommended spinel mesh:** Gamma-centred `4 × 4 × 4`.
- **Recommended beta mesh:** Gamma-centred `12 × 12 × 4`.
- **Production density rule:** use an unshifted Gamma-centred mesh with maximum directional reciprocal-space interval **Δk_max ≤ 0.156 Å^-1** (approximately 0.16 Å^-1).
- **Primary 1 meV/f.u. criterion:** satisfied. The spinel `4 × 4 × 4 → 5 × 5 × 5` change is 0.001481 meV/f.u.; the beta `12 × 12 × 4 → 15 × 15 × 5` change is 0.102590 meV/f.u.
- **Additional convergence jobs:** **NO ADDITIONAL K-POINT TESTING REQUIRED**.
- **Transferability:** **A. STRONGLY SUPPORTED — no additional phase-specific convergence test needed.**

The six completed calculations form two valid convergence sequences. Within each structure series, only the KPOINTS mesh changes; the fixed geometry, electronic Hamiltonian, occupation method, PAW datasets, and all other calculation settings are consistent. The lowest tested mesh in each series is already converged by a wide margin against the next denser tested mesh, so selecting a denser production mesh would add cost without a demonstrated scientific need.

## Calculation and Input Validation

All six calculations completed normally with VASP `5.4.4.18Apr17-6-g9f103f2a35` (build 10 February 2022). Each OUTCAR contains both the explicit `EDIFF`-reached marker and the normal timing/accounting footer. Each `vasprun.xml` parses completely, contains exactly one static calculation, and reproduces the final OUTCAR `free energy TOTEN` exactly at the reported precision. The final electronic cycle required 16–19 iterations, safely below `NELM = 150`; no matched fatal, internal-error, or warning marker was found.

All six INCAR files are byte-identical. The resolved common protocol is:

- PBE (`GGA = PE`) with D3(BJ) (`IVDW = 12`);
- `ENCUT = 500 eV`, `PREC = Accurate`, `EDIFF = 1E-7 eV`, and `ALGO = Normal`;
- tetrahedron occupations (`ISMEAR = -5`; resolved `SIGMA = 0.20 eV` in every OUTCAR);
- non-spin-polarised (`ISPIN = 1`);
- fixed geometry (`IBRION = -1`, `NSW = 0`);
- `LREAL = .FALSE.`, `LASPH = .TRUE.`, `ADDGRID = .TRUE.`, and `ISYM = 2`;
- `ISTART = 0`, `ICHARG = 2`, `LWAVE = .FALSE.`, and `LCHARG = .FALSE.`.

The three spinel POSCAR files are byte-identical to one another, and the three beta POSCAR files are byte-identical to one another. VASP's corresponding CONTCAR files retain the same cells exactly and the same fractional coordinates to within `1.9 × 10^-15`; their textual differences are only signed-zero and last-digit floating-point rewriting. The structures are therefore fixed in both intent and numerical outcome.

The spinel cell contains S32Zn8In16 = 8 ZnIn2S4 formula units. The beta cell contains S4In2Zn1 = 1 ZnIn2S4 formula unit. These counts were read from the actual POSCAR files and are used for every normalized energy below.

Within each series the POTCAR files are byte-identical. Both series use the same per-element PAW datasets: `PAW_PBE S 06Sep2000`, `PAW_PBE Zn 06Sep2000`, and `PAW_PBE In_d 06Sep2000`. Their concatenation order differs only because the POSCAR species orders differ (`S, Zn, In` for spinel and `S, In, Zn` for beta); in every run the POTCAR order matches its POSCAR. This is required bookkeeping, not a Hamiltonian inconsistency.

Every KPOINTS file is Gamma-centred and has zero shift. The mesh in each file agrees with its directory name. No scientifically meaningful input discrepancy was found, and the six calculations genuinely differ only in k-point sampling within their respective structure series.

## Energy Convention

For every run, the energy is the final converged static-calculation `free energy TOTEN` from OUTCAR, independently cross-checked against `vasprun.xml`. The normalized energy is

`E_per_f.u. = TOTEN / N(ZnIn2S4 formula units)`.

For each structure, with meshes ordered from coarsest to densest:

- `ΔE_to_dense = E(mesh) - E(densest tested mesh)`;
- `ΔE_step = E(current mesh) - E(previous mesh)`.

The tables report the absolute magnitudes of these differences in meV per ZnIn2S4 formula unit. No raw total-energy difference is compared across cells with different formula-unit counts.

## Spinel Convergence

| Mesh | Max Δk (Å^-1) | Total Energy (eV) | Energy / f.u. (eV) | ΔE to Densest (meV/f.u.) | ΔE from Previous (meV/f.u.) |
|---|---:|---:|---:|---:|---:|
| 4 × 4 × 4 | 0.147944 | -233.74241459 | -29.2178018238 | 0.011139 | — |
| 5 × 5 × 5 | 0.118355 | -233.74240274 | -29.2178003425 | 0.012620 | 0.001481 |
| 6 × 6 × 6 | 0.098629 | -233.74250370 | -29.2178129625 | 0.000000 | 0.012620 |

The `4 × 4 × 4` mesh differs from the next denser `5 × 5 × 5` mesh by only 0.001481 meV/f.u., and it differs from the densest tested `6 × 6 × 6` mesh by only 0.011139 meV/f.u. Both are negligible relative to the 1.0 meV/f.u. target. The minute non-monotonic variation at `5 × 5 × 5` is bounded within 0.013 meV/f.u. across the complete sequence and does not indicate unresolved convergence.

**Recommended spinel mesh: `4 × 4 × 4`.** A `7 × 7 × 7` test is not scientifically justified by these results.

## Beta Convergence

| Mesh | Max Δk (Å^-1) | Total Energy (eV) | Energy / f.u. (eV) | ΔE to Densest (meV/f.u.) | ΔE from Previous (meV/f.u.) |
|---|---:|---:|---:|---:|---:|
| 12 × 12 × 4 | 0.155894 | -28.75501253 | -28.7550125300 | 0.118590 | — |
| 15 × 15 × 5 | 0.124715 | -28.75511512 | -28.7551151200 | 0.016000 | 0.102590 |
| 18 × 18 × 6 | 0.103929 | -28.75513112 | -28.7551311200 | 0.000000 | 0.016000 |

The `12 × 12 × 4` mesh differs from the next denser `15 × 15 × 5` mesh by 0.102590 meV/f.u. and from the densest tested `18 × 18 × 6` mesh by 0.118590 meV/f.u. Both values are comfortably below 1.0 meV/f.u., with almost a tenfold margin even for the required next-mesh comparison.

**Recommended beta mesh: `12 × 12 × 4`.** A `21 × 21 × 7` test is not scientifically justified by these results.

## Reciprocal-Space Density and Final Static Meshes

Directional intervals were calculated as `|b_i| / N_i` from the actual reciprocal vectors of each relaxed cell, including the `2π` convention. They were not inferred from mesh labels alone. The two independently validated production meshes give Δk_max = 0.147944 Å^-1 for cubic spinel and 0.155894 Å^-1 for layered beta. A practical common production rule is therefore:

**Choose the smallest unshifted Gamma-centred integer mesh satisfying Δk_max ≤ 0.156 Å^-1.**

Applied to the actual accepted relaxed cells, this gives:

| Structure | Recommended Mesh | Δk1, Δk2, Δk3 (Å^-1) | Max Δk (Å^-1) | Rationale |
|---|---|---|---:|---|
| Spinel | 4 × 4 × 4 | 0.147944, 0.147944, 0.147944 | 0.147944 | Directly validated cubic production mesh. |
| Alpha1 | 12 × 12 × 2 | 0.155688, 0.155688, 0.086084 | 0.155688 | Matches the validated in-plane density; one c point would give 0.172168 Å^-1 and miss the target. |
| Beta | 12 × 12 × 4 | 0.155894, 0.155894, 0.129300 | 0.155894 | Directly validated layered production mesh. |
| IIa-prime | 12 × 12 × 2 | 0.155438, 0.155438, 0.127912 | 0.155438 | Approximately twice beta's c repeat, so two c points reproduce beta-like c-axis density. |
| IIb | 12 × 12 × 2 | 0.155690, 0.155690, 0.129009 | 0.155690 | Approximately twice beta's c repeat, so two c points reproduce beta-like c-axis density. |

This strategy deliberately does not copy beta's integer mesh along c to every layered structure. Alpha1's much longer c repeat requires only two c points to remain inside the density limit; IIa-prime and IIb likewise require two rather than four.

## Transferability Assessment

**A. STRONGLY SUPPORTED — no additional phase-specific convergence test needed.**

The beta calculation directly validates the common layered in-plane density at a = 3.878 Å. Alpha1, IIa-prime, and IIb have nearly identical in-plane lattice constants, so `12 × 12` gives Δk1 and Δk2 between 0.155438 and 0.155690 Å^-1. Their differing c repeats are explicitly handled by reciprocal-density matching rather than by assuming a common integer c mesh. The resulting c-axis intervals, 0.086084–0.129009 Å^-1, are no coarser than beta's validated 0.129300 Å^-1.

The transfer also has a substantial numerical safety margin: beta's coarsest-to-next change is 0.102590 meV/f.u., roughly one tenth of the acceptance threshold, and its coarsest-to-densest difference is only 0.118590 meV/f.u. The phases share the same chemistry and planned Hamiltonian, while differences in symmetry and stacking are accommodated by applying the density rule to each actual reciprocal lattice. Together with the independent 3D spinel validation, this is sufficient for the planned five-structure comparison at the stated 1 meV/f.u. standard.

## Additional-Testing Decision

**NO ADDITIONAL K-POINT TESTING REQUIRED**

The existing six calculations establish convergence for both a three-dimensional cubic representative and a layered representative. Neither medium-to-dense change exceeds 1 meV/f.u.; both coarsest meshes also pass against the next denser tested level. No `7 × 7 × 7` spinel, `21 × 21 × 7` beta, or separate alpha1/IIa-prime/IIb spot-check job is warranted by the evidence.

## Exact Strategy for `03_static_scf`

Use the five meshes listed above on the corresponding frozen relaxed structures:

- spinel: `4 × 4 × 4`;
- alpha1: `12 × 12 × 2`;
- beta: `12 × 12 × 4`;
- IIa-prime: `12 × 12 × 2`;
- IIb: `12 × 12 × 2`.

All should be Gamma-centred and unshifted. Use the same fixed-geometry PBE+D3(BJ), 500 eV, `EDIFF = 1E-7 eV`, `ISMEAR = -5`, non-spin-polarised static protocol across all five phases, with each POTCAR concatenated in the exact order of its POSCAR species list. Extract the same final `free energy TOTEN` convention and normalize every energy to one ZnIn2S4 formula unit.

The convergence-test energies in this report are not a final polymorph ranking because they cover only spinel and beta. Final polymorph energies must be taken only from the consistently configured five-structure `03_static_scf` calculations after those calculations have completed successfully.
