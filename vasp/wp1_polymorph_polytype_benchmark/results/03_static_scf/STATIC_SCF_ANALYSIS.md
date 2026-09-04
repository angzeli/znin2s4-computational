# WP1 ZnIn2S4 Static-SCF Analysis

## Executive Summary

All five fixed-geometry calculations completed normally, reached the requested electronic convergence threshold, retained their frozen relaxed structures, produced non-empty CHGCAR files, and gave matching final `free energy TOTEN` values in OUTCAR and vasprun.xml. Four structures receive **PASS** and spinel receives **PASS WITH INTERPRETIVE CAVEAT**; no calculation receives **REQUIRES RECALCULATION**.

Within this 0 K fixed-composition PBE+D3(BJ) dataset, the five-phase order is **spinel < beta < IIb < alpha1 < IIa-prime**. The layered-only order is **beta < IIb < alpha1 < IIa-prime**, with a total layered spread of 29.594 meV/f.u. Spinel is 462.789 meV/f.u. below the lowest layered structure, beta. This separation passes the computational consistency checks but requires later physical interpretation because the experimental spinel source is associated with a high-pressure phase.

The five static energies can be frozen as the authoritative WP1 fixed-geometry PBE+D3(BJ) total-energy dataset, with the stated spinel interpretation caveat. No static-SCF recalculation is required before proceeding to the electronic-structure stage.

## Calculation Consistency

The actual common protocol was VASP 5.4.4 (`vasp.5.4.4.18Apr17-6-g9f103f2a35`, build 10 February 2022), PAW-PBE, `GGA = PE`, `IVDW = 12` (D3(BJ)), `ENCUT = 500 eV`, `PREC = Accurate`, `EDIFF = 1E-7 eV`, `ALGO = Normal` (`IALGO = 38`), `ISMEAR = -5`, `ISPIN = 1`, `IBRION = -1`, `NSW = 0`, `LREAL = .FALSE.`, `LASPH = .TRUE.`, `ADDGRID = .TRUE.`, and `ISYM = 2`. All KPOINTS files specify Gamma-centred, unshifted meshes. `LCHARG = .TRUE.` and each CHGCAR is non-empty.

All calculations use the same per-element datasets: `PAW_PBE Zn 06Sep2000`, `PAW_PBE In_d 06Sep2000`, and `PAW_PBE S 06Sep2000`. The POSCAR species order varies, but each POTCAR concatenation order matches its corresponding POSCAR exactly; this is harmless bookkeeping rather than a Hamiltonian difference.

The spinel and beta runs used the expected `4 x 4 x 4` and `12 x 12 x 4` meshes. Alpha1, IIa-prime, and IIb used `12 x 12 x 4`, rather than the expected `12 x 12 x 2`. This is conservative c-axis oversampling, not undersampling: their maximum reciprocal-grid intervals remain controlled by the `12 x 12` in-plane mesh (0.155438-0.155690 A^-1), while their c-axis intervals are finer (0.043042-0.064504 A^-1). The deviation increases cost but does not reduce the accuracy or comparability of this completed energy set.

The actual POSCAR compositions and cells identify the structures independently of their directory labels: cubic spinel has a 10.6175 A cubic cell and 56 atoms (8 f.u.); the remaining cells have hexagonal metrics and contain 21 atoms for alpha1 (3 f.u.), 7 for beta (1 f.u.), and 14 each for IIa-prime and IIb (2 f.u.). Every static POSCAR is numerically identical to its final geometry-optimisation CONTCAR. Each static CONTCAR retains the input lattice exactly and fractional coordinates to within `3.6 x 10^-15`.

## Static-SCF Validation

The authoritative energy is the final converged OUTCAR `free energy TOTEN`; every value agrees with the corresponding fully parsed vasprun.xml `e_fr_energy` at the printed precision. Each XML contains one fixed-geometry calculation step.

| Structure | f.u. | KPOINTS | Electronic Iterations | TOTEN (eV) | E/f.u. (eV) | Delta E (meV/f.u.) | Residual Pressure/Stress | Verdict |
| --- | ---: | --- | ---: | ---: | ---: | ---: | --- | --- |
| Spinel | 8 | 4 x 4 x 4 | 18 | -233.74241459 | -29.217801824 | 0.000000 | `p_ext = +0.09 kB`; diagonal stress `(+0.092, +0.092, +0.092) kB` | PASS WITH INTERPRETIVE CAVEAT |
| Alpha1 | 3 | 12 x 12 x 4 | 19 | -86.24598186 | -28.748660620 | 469.141204 | `p_ext = +0.10 kB`; diagonal stress `(+0.097, +0.097, +0.099) kB` | PASS |
| Beta | 1 | 12 x 12 x 4 | 16 | -28.75501253 | -28.755012530 | 462.789294 | `p_ext = +0.09 kB`; diagonal stress `(-0.050, -0.050, +0.369) kB` | PASS |
| IIa-prime | 2 | 12 x 12 x 4 | 19 | -57.45083623 | -28.725418115 | 492.383709 | `p_ext = -0.10 kB`; diagonal stress `(-0.097, -0.097, -0.107) kB` | PASS |
| IIb | 2 | 12 x 12 x 4 | 19 | -57.50299979 | -28.751499895 | 466.301929 | `p_ext = -0.29 kB`; diagonal stress `(-0.291, -0.291, -0.295) kB` | PASS |

The formula-unit normalisations follow directly from S32Zn8In16 (spinel), Zn3In6S12 (alpha1), S4In2Zn1 (beta), In4S8Zn2 (IIa-prime), and In4Zn2S8 (IIb).

## Final Polymorph Energy Ranking

Relative to the lowest E/f.u. across all five structures:

1. **Spinel:** 0.000000 meV/f.u.
2. **Beta:** 462.789294 meV/f.u.
3. **IIb:** 466.301929 meV/f.u.
4. **Alpha1:** 469.141204 meV/f.u.
5. **IIa-prime:** 492.383709 meV/f.u.

## Layered Polytype Ranking

Relative to beta, the lowest-energy layered phase:

1. **Beta:** 0.000000 meV/f.u.
2. **IIb:** 3.512635 meV/f.u.
3. **Alpha1:** 6.351910 meV/f.u.
4. **IIa-prime:** 29.594415 meV/f.u.

Required pairwise differences are:

- Alpha1 - beta: **+6.351910 meV/f.u.**
- IIb - beta: **+3.512635 meV/f.u.**
- IIa-prime - beta: **+29.594415 meV/f.u.**
- Alpha1 - IIb: **+2.839275 meV/f.u.**

## Comparison with Lee et al. 2019

The computed layered ordering, `beta < IIb < alpha1 < IIa-prime`, is qualitatively consistent with the reported PBE+D3(BJ) trend `beta <= alpha1 approximately equal to IIb < IIa-prime` in J. Lee et al., *Chemistry of Materials* **2019**, 31, 9148-9155, DOI: 10.1021/acs.chemmater.9b03539. In both cases beta is the lowest layered structure, alpha1 and IIb are near-degenerate, and IIa-prime is higher, with the layered family contained within several tens of meV/f.u.

The quantities are not numerically identical: this work compares fixed-composition static total energies, whereas Lee et al. report elemental formation energies. No absolute formation-energy comparison is made here. Relative ordering is nevertheless meaningful because all structures in this dataset have the same ZnIn2S4 composition and use the same Hamiltonian and energy convention.

## Spinel Forensic Assessment

**Classification B: NUMERICALLY ROBUST BUT PHYSICALLY REQUIRES FOLLOW-UP INTERPRETATION.**

The 462.789 meV/f.u. lowering of spinel relative to beta is reproduced by the final static calculations and is not attributable to an identified bookkeeping or protocol error:

- all phases use the same VASP build, PBE functional, D3(BJ) treatment, 500 eV cutoff, spin treatment, precision settings, and energy convention;
- the same Zn, In_d, and S PAW datasets are used, and each POTCAR order matches its POSCAR;
- spinel has the correct Zn8In16S32 composition and is normalised by 8 ZnIn2S4 formula units;
- its OUTCAR and vasprun.xml energies agree, EDIFF was reached in 18 iterations, and the run terminated normally;
- its static input is exactly the relaxed spinel geometry and remained fixed;
- its directly validated `4 x 4 x 4` mesh has a maximum reciprocal interval of 0.147944 A^-1, comparable to or finer than the layered maximum intervals of 0.155438-0.155894 A^-1;
- its residual pressure and stress are only about 0.09 kB, so the energy is not associated with a large residual cell stress.

**Computational result:** Lee et al. independently report essentially the same spinel stabilisation in Supporting Information Table S3. Their PBE+D3(BJ) formation energies are approximately -4.59 eV/f.u. for spinel and -4.13 eV/f.u. for beta, giving `E_spinel - E_beta` of approximately **-460 meV/f.u.**; the present static-energy difference is **-462.789 meV/f.u.** This is quantitative agreement at the precision of the tabulated literature values. Their LDA, PBE, and optB86b results also favour spinel over beta by approximately 450, 300, and 430 meV/f.u., respectively. The large computed separation is therefore not anomalous relative to that independent study.

**Physical interpretation:** This computational agreement does **not** establish spinel as the experimentally stable ambient-pressure phase. The present calculation compares 0 K internal energies at separately relaxed cells; it does not include finite-temperature vibrational free energies, pressure-dependent enthalpies, kinetic accessibility, precursor-dependent synthesis effects, or experimental phase abundance. Because the experimental source spinel is associated with a high-pressure ZnIn2S4 phase, the relationship between the static DFT ordering and its experimentally described high-pressure character remains unresolved. Pressure-dependent and/or thermodynamic follow-up may be useful if this issue becomes publication-relevant, but the present data do not establish a unique physical explanation and none is inferred here.

## Residual Stress Check

The OUTCAR-reported final external pressures range from -0.29 to +0.10 kB, with `PSTRESS = 0.0` in every run. The largest individual normal stress magnitude is 0.369 kB (beta), and all reported shear components are zero at the displayed precision. These values are small numerical residual stresses (at most approximately 0.037 GPa), not evidence of a physically applied pressure or an inadequately relaxed cell. All frozen geometries remain effectively stress-relaxed for this static comparison.

## Acceptance Decision

- **Accepted structures (PASS):** alpha1, beta, IIa-prime, IIb.
- **Accepted with interpretive caveat (PASS WITH INTERPRETIVE CAVEAT):** spinel.
- **Structures requiring recalculation (REQUIRES RECALCULATION):** none.

The five energies can be frozen as the authoritative fixed-geometry PBE+D3(BJ) WP1 energy dataset. The spinel value is computationally accepted but must retain its pressure/thermodynamic interpretation caveat in subsequent scientific discussion.

## Recommended Next Step

Proceed to `04_electronic_structure` for DOS/PDOS and PBE band structures, followed by `05_hybrid_validation` for HSE06 electronic validation. No additional static-SCF calculation is required before that progression. Pressure-dependent and finite-temperature work may be considered later if the physical origin of the spinel ordering becomes part of the project scope; those calculations are not prerequisites for accepting this static dataset and were not performed here.
