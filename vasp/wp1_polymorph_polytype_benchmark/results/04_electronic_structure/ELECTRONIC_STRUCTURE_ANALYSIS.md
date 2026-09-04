# WP1 ZnIn2S4 PBE Electronic-Structure Analysis

## Executive Summary

All 10 Stage 04 calculations (five dense DOS/PDOS branches and five high-symmetry band branches) passed termination, electronic-convergence, fixed-geometry, structure/charge-density-lineage, dimensional, and parseability checks.

The dense-mesh sampled PBE fundamental gaps span 0.279-1.318 eV. IIb has the smallest sampled gap (0.279 eV), and Spinel has the largest (1.318 eV).

Sampled classifications comprise 2 direct, 2 indirect, and 1 near-degenerate cases. S-p is the dominant VBM projection in every phase. The layered CBMs are strongly hybridized S-p/In-s/S-s states whose first-place channel is window-dependent; the spinel CBM instead has substantial In-s/S-s/Zn-s character.

Phase-level acceptance is 4 PASS, 1 PASS WITH CAVEAT, and 0 REQUIRES RECALCULATION. No Stage 04 recalculation is required before HSE06. Stage 04 can be frozen, with the finite uniform-mesh extremum and PBE gap-magnitude limitations carried forward explicitly into Stage 05.

All values below are **PBE electronic structures on PBE+D3(BJ)-relaxed geometries**. D3(BJ) determined the geometry and is not an electronic band-gap correction.

## Calculation Validation

Both branches use the accepted Stage 03 POSCAR and a byte-identical Stage 03 CHGCAR. `NSW = 0` and `IBRION = -1` freeze the geometry; output CONTCAR coordinates differ from POSCAR only at numerical round-off. Both DOS/PDOS and band branches are fixed-density non-self-consistent evaluations (`ICHARG = 11`). This is valid here because the charge-density/structure lineage is exact, but it must not be described as a new dense-mesh self-consistent charge calculation.

| Phase | Branch | k sampling | Irreducible/path k points | Final electronic iteration | Validation |
| --- | --- | --- | ---: | ---: | --- |
| Spinel | DOS/PDOS | Gamma 6 x 6 x 6 | 20 | 11 | PASS |
| Spinel | Band | cF2 Line-mode, 40/segment | 240 | 11 | PASS |
| Alpha1 | DOS/PDOS | Gamma 18 x 18 x 6 | 202 | 11 | PASS |
| Alpha1 | Band | hR1 Line-mode, 40/segment | 280 | 11 | PASS |
| Beta | DOS/PDOS | Gamma 18 x 18 x 6 | 202 | 10 | PASS |
| Beta | Band | hP2 Line-mode, 40/segment | 360 | 10 | PASS |
| IIa-prime | DOS/PDOS | Gamma 18 x 18 x 6 | 202 | 10 | PASS |
| IIa-prime | Band | hP2 Line-mode, 40/segment | 360 | 10 | PASS |
| IIb | DOS/PDOS | Gamma 18 x 18 x 6 | 148 | 10 | PASS |
| IIb | Band | hP2 Line-mode, 40/segment | 360 | 10 | PASS |

All INCARs resolve to PBE (`GGA = PE`), 500 eV, `ISPIN = 1`, `LORBIT = 11`, and `EDIFF <= 1e-7` eV. DOS/PDOS uses the stated unshifted dense Gamma mesh, `ISMEAR = -5`, `NEDOS = 3000`, and `ISYM = 2`; band calculations use explicit reciprocal-coordinate Line-mode, `ISMEAR = 0`, and `ISYM = 0`. Every OUTCAR contains both the EDIFF and normal timing markers, every vasprun.xml closes normally, and every EIGENVAL/PROCAR/DOSCAR dimension agrees with the structure. Edge energies and occupied/unoccupied assignments independently agree between EIGENVAL and PROCAR at every sampled k point.

## Band-Gap Summary

| Phase | Uniform-mesh fundamental gap (eV) | Minimum direct gap (eV) | Delta direct-indirect (meV) | Classification | Sampled VBM | Sampled CBM |
| --- | ---: | ---: | ---: | --- | --- | --- |
| Spinel | 1.318 | 1.428 | 110.2 | INDIRECT ON THE SAMPLED UNIFORM MESH | (0.500000, 0.500000, 0.000000) | GAMMA |
| Alpha1 | 0.304 | 0.304 | 0.0 | DIRECT ON THE SAMPLED UNIFORM MESH | GAMMA; (0.000000, 0.000000, 0.166667); (0.000000, 0.000000, 0.333333); (0.000000, 0.000000, 0.500000) | GAMMA |
| Beta | 0.303 | 0.304 | 1.4 | AMBIGUOUS / NEAR-DEGENERATE | A | GAMMA |
| IIa-prime | 0.304 | 0.346 | 42.0 | INDIRECT ON THE SAMPLED UNIFORM MESH | (0.111111, 0.055556, 0.000000) | GAMMA |
| IIb | 0.279 | 0.279 | 0.0 | DIRECT ON THE SAMPLED UNIFORM MESH | GAMMA | GAMMA |

Classification uses the dense uniform mesh. A direct-minus-fundamental difference of <= 10 meV is treated as near-degenerate unless the same sampled k point explicitly hosts both extrema. The CSV retains raw VASP eigenvalue references; absolute eigenvalues from separate cells are not vacuum-aligned band positions and must not be compared as such.

The path-minus-uniform gap differences (spinel, alpha1, beta, IIa-prime, IIb) are -3.4, -0.0, +0.0, -4.6, +0.0 meV. The small negative differences for spinel and IIa-prime mean that the finer Line-mode sampling intersects a slightly higher valence point than the uniform mesh; they are retained explicitly rather than treated as corrupt data.

## High-Symmetry Band Structures

### Spinel

Path: `GAMMA-X-U | K-GAMMA-L-W-X` (cF2; 40 points per segment).

The path-only gap is 1.315 eV, with the path VBM at K-GAMMA (t=0.256); K-GAMMA (t=0.410) and path CBM at GAMMA; X. The top-valence/bottom-conduction ranges are 0.021/0.859 eV on GAMMA-X and 0.080/1.219 eV on GAMMA-L.

### Alpha1

Path: `GAMMA-T-H_2 | H_0-L-GAMMA-S_0 | S_2-F-GAMMA` (hR1; 40 points per segment).

The path-only gap is 0.304 eV, with the path VBM at GAMMA-T (entire segment) and path CBM at GAMMA; GAMMA-T (t=0.667). The top-valence/bottom-conduction ranges are 0.000/0.195 eV on GAMMA-T and 0.627/1.203 eV on L-GAMMA; the rhombohedral path geometry prevents treating this as a simple Cartesian mass comparison.

### Beta

Path: `GAMMA-M-K-GAMMA-A-L-H-A | L-M | H-K` (hP2; 40 points per segment).

The path-only gap is 0.303 eV, with the path VBM at GAMMA-A (t=0.949-1.000, near A) and path CBM at GAMMA. On GAMMA-M, the top-valence/bottom-conduction band ranges are 0.604/1.232 eV; on GAMMA-A they are 0.001/0.674 eV. These path-wise ranges indicate anisotropy but are not effective masses.

### IIa-prime

Path: `GAMMA-M-K-GAMMA-A-L-H-A | L-M | H-K` (hP2; 40 points per segment).

The path-only gap is 0.299 eV, with the path VBM at K-GAMMA (t=0.744) and path CBM at GAMMA. On GAMMA-M, the top-valence/bottom-conduction band ranges are 0.688/1.330 eV; on GAMMA-A they are 0.000/0.170 eV. These path-wise ranges indicate anisotropy but are not effective masses.

### IIb

Path: `GAMMA-M-K-GAMMA-A-L-H-A | L-M | H-K` (hP2; 40 points per segment).

The path-only gap is 0.279 eV, with the path VBM at GAMMA-A (t=0.000-0.128, near GAMMA) and path CBM at GAMMA. On GAMMA-M, the top-valence/bottom-conduction band ranges are 0.596/1.162 eV; on GAMMA-A they are 0.000/0.366 eV. These path-wise ranges indicate anisotropy but are not effective masses.

## DOS / PDOS

DOS onsets are used only as a qualitative check on the eigenvalue-derived gaps. Band-edge character is the integrated LORBIT=11 PAW projection over 0.50 eV edge windows; the machine-readable table also includes 0.20 eV windows. Percentages are normalized within each projected window and are comparative descriptors, not exact chemical populations.

| Phase | E-fermi (eV) | DOS valence/conduction onset (eV) | Gross occupied DOS width (eV) | Conduction TDOS states/f.u. (0-0.5 / 0.5-1.0 eV) | Dominant VBM character | Dominant CBM character | Notes |
| --- | ---: | --- | ---: | --- | --- | --- | --- |
| Spinel | 4.182 | 4.172 / 5.583 | 14.13 | 0.02 / 0.37 | S-p 89.1%, In-d 6.1%, Zn-d 2.8% | In-s 39.6%, S-s 21.9%, S-p 17.5% | Clean sampled DOS gap; No individual same-element site dominates the 0.50 eV edge windows. |
| Alpha1 | 4.519 | 4.510 / 4.847 | 15.71 | 0.02 / 0.44 | S-p 76.1%, Zn-d 12.2%, Zn-p 6.5% | S-p 42.5%, In-s 32.7%, S-s 19.1% | Clean sampled DOS gap; VBM: S3 carries 21.8% of total projected edge weight. |
| Beta | 4.534 | 4.524 / 4.883 | 15.66 | 0.02 / 0.37 | S-p 76.9%, Zn-d 11.3%, Zn-p 5.8% | S-p 43.9%, In-s 31.7%, S-s 18.4% | Clean sampled DOS gap; VBM: S4 carries 62.4% of total projected edge weight; CBM: S1 carries 32.7% of total projected edge weight; CBM: In1, S4 carries 23.6% of total projected edge weight. |
| IIa-prime | 4.164 | 4.161 / 4.501 | 15.55 | 0.03 / 0.13 | S-p 78.8%, Zn-d 13.7%, Zn-p 4.0% | S-p 44.1%, In-s 30.0%, S-s 19.8% | Clean sampled DOS gap; VBM: S8 carries 34.6% of total projected edge weight. |
| IIb | 4.512 | 4.503 / 4.817 | 15.72 | 0.02 / 0.38 | S-p 75.8%, Zn-d 12.4%, Zn-p 6.6% | S-p 42.8%, In-s 32.5%, S-s 18.9% | Clean sampled DOS gap; VBM: S1 carries 32.0% of total projected edge weight. |

The DOS Fermi/reference energies agree between OUTCAR, vasprun.xml, and DOSCAR within 1 meV. Integrated DOS at the Fermi level recovers the actual electron count, supporting the occupied-band assignment `NELECT / 2` for these non-spin-polarized runs.

The narrower 0.20 eV integrations preserve the S-p-dominated VBM assignment. They also preserve the mixed S-p/In-s/S-s layered CBM and In-s/S-s/Zn-s spinel CBM descriptions, so the qualitative edge-character conclusions are not artifacts of choosing the 0.50 eV reporting window.

## Cross-Phase Electronic Trends

The sampled PBE gap ordering is IIb (0.279 eV) < Beta (0.303 eV) < IIa-prime (0.304 eV) < Alpha1 (0.304 eV) < Spinel (1.318 eV).

VBM character is broadly conserved: S-p is dominant in all five phases.

The layered CBMs consistently mix S-p, In-s, and S-s weight; S-p and In-s are close enough that the leading label should not be interpreted as a pure orbital assignment. Spinel is distinct in retaining a much larger Zn-s contribution near the CBM.

Among the layered stacking variants, beta and IIb differ by 0.024 eV, while alpha1 differs from beta by 0.001 eV. Similar dominant edge channels support a common local electronic motif, but the nonzero gap and dispersion changes show that stacking is not electronically invisible.

IIa-prime does not introduce a new leading edge-orbital species, but it is distinct through its off-GAMMA sampled VBM, 42.0 meV direct-minus-fundamental separation, only 0.170 eV bottom-conduction-band range on GAMMA-A, and concentration of 34.6% of the 0.50 eV VBM projection on S8. These are topology-associated signatures, not proof of a single causal mechanism.

The all-phase gap-density Pearson coefficient is 0.98, but it is dominated by the dense, wide-gap spinel. Within the four layered phases it is -0.14, so density does not provide a robust standalone explanation of their gap ordering. The sample is too small for a causal claim.

Spinel is retained as a qualitatively distinct three-dimensional reference. IIa-prime is treated separately because its reconstructed layer topology can change site weighting and dispersion even where the leading element/orbital labels remain the same.

## Literature Context

Lee et al. focused the accurate electronic-structure comparison of the layered polytypes on hybrid-functional results and proposed the revised IIa-prime structure. Accordingly, no direct numerical literature match is claimed for these present PBE gaps. The Stage 04 results are used to define hypotheses for HSE06, not as substitutes for the literature HSE06 values. See [Lee et al., Chemistry of Materials 31, 9148-9155 (2019)](https://doi.org/10.1021/acs.chemmater.9b03539).

## Direct vs Indirect Behaviour

The primary fundamental gap is `min_k E_CBM(k) - max_k E_VBM(k)` on the dense uniform mesh. The minimum sampled direct gap is `min_k [E_CBM(k) - E_VBM(k)]`. Their difference is reported explicitly. High-symmetry path extrema are analysed separately and never substituted for the full-mesh sampled extrema.

These are finite-k sampled assignments, not mathematical proofs of the continuous-zone global extrema. Symmetry reduction means each reported uniform coordinate is an irreducible representative; equivalent full-zone points are implicit. HSE06 should recheck the important direct/indirect assignments rather than assume a rigid scissor shift.

## Stage 05 HSE06 Validation Targets

| Phase | PBE gap (eV) | PBE sampled classification | Sampled VBM | Sampled CBM | Dominant VBM | Dominant CBM | Main HSE06 question |
| --- | ---: | --- | --- | --- | --- | --- | --- |
| Spinel | 1.318 | INDIRECT ON THE SAMPLED UNIFORM MESH | (0.500000, 0.500000, 0.000000) | GAMMA | S-p 89.1% | In-s 39.6%, S-s 21.9%, S-p 17.5% | Does the three-dimensional reference remain electronically distinct and retain its extremum assignment? |
| Alpha1 | 0.304 | DIRECT ON THE SAMPLED UNIFORM MESH | GAMMA; (0.000000, 0.000000, 0.166667); (0.000000, 0.000000, 0.333333); (0.000000, 0.000000, 0.500000) | GAMMA | S-p 76.1% | S-p 42.5%, In-s 32.7%, S-s 19.1% | Does HSE06 preserve the sampled direct/indirect class and layered-phase gap ordering? |
| Beta | 0.303 | AMBIGUOUS / NEAR-DEGENERATE | A | GAMMA | S-p 76.9% | S-p 43.9%, In-s 31.7%, S-s 18.4% | Does HSE06 resolve the near-degenerate direct/indirect assignment without moving an extremum off the sampled point? |
| IIa-prime | 0.304 | INDIRECT ON THE SAMPLED UNIFORM MESH | (0.111111, 0.055556, 0.000000) | GAMMA | S-p 78.8% | S-p 44.1%, In-s 30.0%, S-s 19.8% | Does the revised layer topology retain its gap, edge locations, and site-weighting signature? |
| IIb | 0.279 | DIRECT ON THE SAMPLED UNIFORM MESH | GAMMA | GAMMA | S-p 75.8% | S-p 42.8%, In-s 32.5%, S-s 18.9% | Does HSE06 preserve the sampled direct/indirect class and layered-phase gap ordering? |

Across the phase set, HSE06 should test gap ordering, sampled direct/indirect character, edge-location stability, and whether the leading orbital-projection picture survives. SOC remains deferred; Stage 04 supplies no specific result that changes that project decision.

## Acceptance Decision

| Phase | Decision | Basis |
| --- | --- | --- |
| Spinel | PASS | Complete, converged, frozen, lineage-consistent, parseable, semiconducting, and internally cross-checked. |
| Alpha1 | PASS | Complete, converged, frozen, lineage-consistent, parseable, semiconducting, and internally cross-checked. |
| Beta | PASS WITH CAVEAT | Complete, converged, frozen, lineage-consistent, parseable, semiconducting, and internally cross-checked. |
| IIa-prime | PASS | Complete, converged, frozen, lineage-consistent, parseable, semiconducting, and internally cross-checked. |
| IIb | PASS | Complete, converged, frozen, lineage-consistent, parseable, semiconducting, and internally cross-checked. |

**Stage 04 freeze decision: PASS.** All five phases can be frozen for the stated PBE scope.

## Recommended Next Step

Proceed to the already planned Stage 05 HSE06 validation. No new PBE calculation is required before HSE06. Preserve the sampled-language and fixed-density provenance when comparing the two stages, and do not interpret HSE06 as a guaranteed rigid correction to every PBE band edge.
