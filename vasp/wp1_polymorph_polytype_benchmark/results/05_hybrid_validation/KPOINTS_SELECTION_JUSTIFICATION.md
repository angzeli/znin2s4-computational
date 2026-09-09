# WP1 HSE06 Band KPOINTS Selection Justification

## Purpose

Stage 05 validates PBE gap magnitudes, sampled direct/indirect assignments, VBM/CBM locations, and phase trends for spinel, alpha1, beta, IIa_prime, and IIb. Stage 05 production now targets the validated native VASP 6.6.1 build.

This document distinguishes the scientific selection (regular mesh, canonical reciprocal-space path, and 20 points per segment) from its file representation. The migration preserves that selection and the original frozen POSCARs. It does not establish completed production HSE06 results or independent HSE convergence.

## Current VASP 6.6.1 Strategy

The **VASP 6.6.1 production representation** uses two files in each `calculation/05_hybrid_validation/02_hse06_band/<phase>/` directory:

- `KPOINTS`: the automatic, unshifted Gamma-centered regular mesh for the HSE SCF.
- `KPOINTS_OPT`: the additional high-symmetry Line-mode path, in reciprocal coordinates, with 20 points per segment, evaluated after regular-mesh self-consistency.

VASP reads KPOINTS_OPT when present unless explicitly disabled. Its KPOINTS partner must supply a uniform mesh. See the [official KPOINTS_OPT documentation](https://vasp.at/wiki/KPOINTS_OPT).

The **VASP 5.4.4 legacy/reference representation** placed weighted irreducible regular points and zero-weight path points together in one explicit KPOINTS file. That method was valid and independently validated. It is superseded as the production representation because the target native build supports KPOINTS_OPT; it remains a useful reference in the provenance sections below. HSE06 remains orbital dependent, so ICHARG=11 is not used.

## Regular-Mesh Selection

| Phase | Regular mesh | Approx reciprocal spacings (Angstrom^-1) | Basis |
| --- | --- | --- | --- |
| spinel | 4 x 4 x 4 | 0.147944, 0.147944, 0.147944 | Accepted Stage 02/03 sampling; matches Stage 05 bootstrap/SCF |
| alpha1 | 12 x 12 x 4 | 0.155688, 0.155688, 0.043042 | Accepted Stage 02/03 sampling; matches Stage 05 bootstrap/SCF |
| beta | 12 x 12 x 4 | 0.155894, 0.155894, 0.129300 | Accepted Stage 02/03 sampling; matches Stage 05 bootstrap/SCF |
| IIa_prime | 12 x 12 x 4 | 0.155438, 0.155438, 0.063956 | Accepted Stage 02/03 sampling; matches Stage 05 bootstrap/SCF |
| IIb | 12 x 12 x 4 | 0.155690, 0.155690, 0.064504 | Accepted Stage 02/03 sampling; matches Stage 05 bootstrap/SCF |

All meshes are Gamma centered with shift `0 0 0`. The spacings are the three values norm(b_i)/N_i using reciprocal vectors with the 2*pi convention. They remain the previously validated values because the lattices and mesh dimensions are unchanged, and satisfy the accepted approximately 0.156 Angstrom^-1 PBE production criterion.

The Stage 02 convergence study and Stage 03 static-SCF framework established the regular sampling selection. Longer layered c repeats receive finer c-axis sampling with the same integer divisions. PBE energy convergence does not by itself establish convergence of HSE06 exchange, gaps, or band extrema.

The automatic meshes agree across Stage 05-00, 05-01, and 05-02, as do the exact frozen POSCARs and atom ordering within each phase. The completed PBE bootstrap retains ISYM=2; both VASP 6.6.1 HSE stages now explicitly use ISYM=3, the [documented hybrid default and orbital-based symmetry strategy](https://vasp.at/wiki/ISYM). Native initialization-only runs measured regular counts of spinel 10, alpha1 69, beta 69, IIa_prime 69, and IIb 57. Full crystal-symmetry/time-reversal orbit equivalence, complete grid coverage and exact multiplicity agreement with the donors were verified; identical representative coordinates or ordering are not required. These are measured initialization counts, not production HSE convergence results. The underlying full Gamma meshes and optional path coordinates are unchanged.

## High-Symmetry Paths

| Phase | Space group | SeeK-path type | Canonical path | KPOINTS_OPT points/segment |
| --- | --- | --- | --- | --- |
| spinel | Fd-3m (227) | cF2 | GAMMA-X-U &#124; K-GAMMA-L-W-X | 20 |
| alpha1 | R3m (160) | hR1 | GAMMA-T-H_2 &#124; H_0-L-GAMMA-S_0 &#124; S_2-F-GAMMA | 20 |
| beta | P3m1 (156) | hP2 | GAMMA-M-K-GAMMA-A-L-H-A &#124; L-M &#124; H-K | 20 |
| IIa_prime | P-3m1 (164) | hP2 | GAMMA-M-K-GAMMA-A-L-H-A &#124; L-M &#124; H-K | 20 |
| IIb | P63mc (186) | hP2 | GAMMA-M-K-GAMMA-A-L-H-A &#124; L-M &#124; H-K | 20 |

The authoritative coordinate sources are the five validated `calculation/04_electronic_structure/<phase>/band/KPOINTS` files. Each KPOINTS_OPT copies their endpoint coordinate rows, labels, segment pairing, and blank separators exactly. Only the descriptive first line and the points-per-segment value differ from those Stage 04 files.

Line-mode interpolates each endpoint pair separately. A branch break is preserved by ending one pair and beginning the next at the required distinct branch start; no segment bridges U to K for spinel, H_2 to H_0 or S_0 to S_2 for alpha1, or the two hP2 breaks.

Line-mode includes both endpoints of each segment, so connected endpoints appear twice in the expanded list. At 20 points per segment, the nominal expanded path counts are 120 for spinel, 140 for alpha1, and 180 for each hP2 phase. After removing only the consecutive shared endpoints, these correspond to the same 116, 136, and 174 sampled path points used by the legacy construction. Repeated coordinates are a representation difference, not a change in the path or its resolution.

## SeeK-path Validation

The endpoint selection retains the prior independent validation using official SeeK-path 2.2.1 and spglib 2.7.0, symprec=1e-5, with time reversal. The recorded maximum official Cartesian path mismatch was **1.385866231e-10 Angstrom^-1**, approximately 1.386e-10 and below the 1e-8 validation target.

This migration uses those same already-validated Stage 04 coordinates. No new crystallographic path selection or SeeK-path run was needed. Exact equality of all endpoint rows in each new KPOINTS_OPT and its Stage 04 source gives zero introduced fractional or Cartesian endpoint displacement for the unchanged lattice.

The following numerical evidence belongs to the prior explicit construction and remains its validation record:

| Phase | Max endpoint error vs official SeeK-path (Angstrom^-1) | Max path error vs official Cartesian interpolation (Angstrom^-1) | Max interpolation/serialization error vs Stage 04 (Angstrom^-1) |
| --- | --- | --- | --- |
| spinel | 2.419675e-16 | 4.904026e-15 | 4.855481e-15 |
| alpha1 | 1.385866e-10 | 1.385866e-10 | 1.537579e-14 |
| beta | 1.141949e-10 | 1.141950e-10 | 1.529594e-14 |
| IIa_prime | 1.138577e-10 | 1.138578e-10 | 1.531303e-14 |
| IIb | 1.140638e-10 | 1.140639e-10 | 1.526497e-14 |

The prior interpolation/serialization mismatch against Stage 04 was at most 1.537579e-14 Angstrom^-1; the larger official mismatch reflects finite decimal precision of the original endpoints.

Before replacement, the committed legacy KPOINTS files at `84c45af` were checked against their Stage 03 weighted rows and the 20-point Stage 04-derived interpolation. Counts, positive multiplicities, exact zero path weights, branch structure, and every path coordinate agreed. Maximum fractional rounding residual against decimal interpolation was 4.737e-15.

## Original POSCAR Reciprocal Basis

No POSCAR was standardized, rotated, or rewritten. For original real-space row-vector lattice A, B=2*pi*inverse(A).transpose() is the reciprocal lattice in Angstrom^-1, and a fractional reciprocal vector q corresponds to qB. KPOINTS_OPT uses that unchanged basis.

During the previous independent validation, SeeK-path constructed a standard primitive cell internally. Standard primitive reciprocal endpoint vectors were converted to Cartesian coordinates and rotated back using its returned rotation matrix before comparison in the original POSCAR frame. This transformation was part of validation, not a rewrite of any structure.

The spinel and alpha1 input cells contain respectively four and three primitive cells. Their canonical primitive-cell HPKOT labels are retained in the original conventional-cell reciprocal basis, where bands are folded and some coordinates lie beyond that cell's first zone. Coordinates were not wrapped or relabelled. The hP2 cells have primitive-cell volume ratio one.

## Why KPOINTS_OPT

In the legacy explicit hybrid workflow, the regular and zero-weight path states participate in the same hybrid iteration sequence. With KPOINTS_OPT, the regular HSE SCF uses KPOINTS and the additional path states are evaluated afterward. This avoids embedding the complete band path in every regular hybrid SCF iteration. It is a workflow change, without a claimed numerical speedup. See the [official hybrid-band workflow description](https://vasp.at/wiki/Band-structure_calculation_using_hybrid_functionals).

The project context supplied for this migration records successful HSE06, ACE/LFOCKACE, and reduced-beta HSE KPOINTS_OPT preflight on the native VASP 6.6.1 port. That runtime evidence was not rerun here. Full production-mesh performance remains a separate validation task.

## Relationship to the Legacy VASP 5.4.4 Representation

The following counts describe only the **VASP 5.4.4 legacy/reference representation**, preserved in Git commit `84c45af`:

| Phase | Legacy weighted IBZ points | Legacy zero-weight path points | Legacy total | Regular weight sum |
| --- | --- | --- | --- | --- |
| spinel | 10 | 116 | 126 | 64 |
| alpha1 | 69 | 136 | 205 | 576 |
| beta | 69 | 174 | 243 | 576 |
| IIa_prime | 69 | 174 | 243 | 576 |
| IIb | 57 | 174 | 231 | 576 |

All legacy weighted blocks came from `calculation/03_static_scf/<phase>/IBZKPT`. At the original construction, no Stage 05-00/01 IBZKPT was available. The Stage 03 and Stage 05 geometries, regular meshes, and SCF symmetry settings were verified to agree. Stage 03 OUTCARs recorded electronic convergence and normal completion.

The prior independent spglib 2.7.0 check reproduced every VASP irreducible orbit and multiplicity under the SCF ISYM=2, ISPIN=1, no-SOC/noncollinear, no-custom-MAGMOM settings, with symprec=1e-5 and time reversal. VASP representatives were mapped into the full grid and spglib orbits without requiring identical representative ordering. Weight sums were 64 or 576, with complete orbit coverage.

Legacy weighted rows were preserved verbatim, excluding tetrahedron metadata, and followed by 14-decimal zero-weight path coordinates. That construction included both endpoints of each branch's first segment and omitted shared starts on subsequent connected segments. With S segments and B branches, N_path=19S+B: 116 for 6/2, 136 for 7/3, and 174 for 9/3. Weighted/path coordinate coincidences were intentional.

For reproducibility, the legacy endpoint mapping is retained below. These are one-based indices into the old combined explicit list, not KPOINTS_OPT indices; old physical file line=index+3.

| Phase | Legacy endpoint index:label sequence |
| --- | --- |
| spinel | 11:GAMMA; 30:X; 49:U; 50:K; 69:GAMMA; 88:L; 107:W; 126:X |
| alpha1 | 70:GAMMA; 89:T; 108:H_2; 109:H_0; 128:L; 147:GAMMA; 166:S_0; 167:S_2; 186:F; 205:GAMMA |
| beta | 70:GAMMA; 89:M; 108:K; 127:GAMMA; 146:A; 165:L; 184:H; 203:A; 204:L; 223:M; 224:H; 243:K |
| IIa_prime | 70:GAMMA; 89:M; 108:K; 127:GAMMA; 146:A; 165:L; 184:H; 203:A; 204:L; 223:M; 224:H; 243:K |
| IIb | 58:GAMMA; 77:M; 96:K; 115:GAMMA; 134:A; 153:L; 172:H; 191:A; 192:L; 211:M; 212:H; 231:K |

The counts and index mappings are provenance cross-checks for the same scientific selection. They are no longer the current KPOINTS file format, and no legacy or backup copies have been added to the repository.

## Why 20 Points per Segment

Twenty points per segment remains the bounded Stage 05 validation-resolution choice. Stage 04 PBE uses 40 points per segment. This migration neither promotes the HSE path to 40 points nor claims that 20 is independently path converged.

A future 20-versus-40 KPOINTS_OPT comparison may be performed if results indicate insufficient resolution, especially for near-degenerate extrema. Finite meshes and high-symmetry paths do not prove the exact continuous-zone global extrema.

## VASP 6.6.1 Runtime Constraints

The project preflight context establishes KPOINTS_OPT functionality on the native VASP 6.6.1 build and requires NCORE=1 for the validated local hybrid KPOINTS_OPT workflow. Production HSE mesh convergence and runtime remain separate validation tasks.

The production dependency is **05-00 PBE WAVECAR → 05-01 regular-mesh HSE06 SCF → converged HSE WAVECAR → 05-02 HSE06 restart + KPOINTS_OPT**. The previous direct PBE-to-band route belonged to the valid VASP 5.4.4 explicit zero-weight implementation. It is no longer the configured production route. KPOINTS_OPT evaluates the unchanged canonical path after regular HSE self-consistency; splitting the stages is a provenance/reuse choice, not a claimed speed advantage.

Both HSE branches specify LHFCALC=.TRUE., GGA=PE, HFSCREEN=0.2, ISTART=1, ICHARG=0, ALGO=Normal, LFOCKACE=.TRUE., HFRCUT=-1, ISYM=3, ISPIN=1, and fixed geometry (NSW=0, IBRION=-1). TIME was removed with the Damped-to-Normal change. Explicit NBANDS is 304/112/40/80/80 for spinel/alpha1/beta/IIa_prime/IIb in both branches, matching the PBE donors. These counts support the current band-edge validation scope, not independently converged high-energy conduction states. All unrelated convergence and output controls remain unchanged; 05-01 retains LWAVE=.TRUE. to supply the next stage.

NCORE, NPAR, and KPAR remain absent from all ten canonical HSE INCARs. The frozen generic port launcher owns execution-layer parallelism; the controlled validation selects --ranks 8 --ncore 1 --kpar 1 --mpi-mode synthetic. Only staged execution copies receive those parallel tags. The local ignored `vasp5.4.sub` scripts remain legacy VASP 5.4.4 launchers and are not the native production entry point.

The initial KPOINTS-only migration did not edit INCARs. The subsequent restart alignment explicitly adopts Davidson/ACE, supported by the [official ACE documentation](https://vasp.at/wiki/LFOCKACE) and the port's bounded runtime evidence. HFRCUT=-1 follows the [official gapped hybrid-band recommendation](https://vasp.at/wiki/Band-structure_calculation_using_hybrid_functionals) and [HFRCUT definition](https://vasp.at/wiki/HFRCUT); applying it to both HSE branches keeps the finite-mesh Coulomb treatment consistent. It changes that numerical treatment from the previous default, not the intended HSE06 exchange fraction or screening. The frozen launcher was not edited. No HSE production convergence or speedup is established; the [bootstrap audit](05_00_PBE_WAVECAR_BOOTSTRAP_AUDIT.md) records the separate beta restart-control evidence and remaining runtime gates.

## ISYM=3 Runtime Provenance — 2026-09-09

All five updated 05-01 inputs passed the frozen launcher's eight-rank/NCORE=1/KPAR=1/synthetic-MPI dry-runs with original matching PBE WAVECAR donors. The unchanged Scratch orchestrator also passed its no-launch dry-run; its stale informational ISYM=2 caveat does not override the actual ISYM=3 input summaries.

Native initialization generated 10/69/69/69/57 points. Spinel and IIb directly match all donor representatives; alpha1, beta and IIa_prime each directly match 57/69, with all 69 mapping uniquely under crystal symmetry/time reversal. Orbit sets are disjoint, each multiplicity equals its orbit size, and full-grid coverage is 64/64 or 576/576 for both donor and target. The maximum symmetry-mapped coordinate residual is 1.001e-14. An otherwise-identical beta ISYM=2 control produced byte-identical IBZKPT to ISYM=3: the representative difference is not caused by the ISYM switch in that comparison.

The first plain-PBE/ISYM=3 control read WAVECAR but failed during full-k-point initialization; this failed attempt is preserved. A [VASP-developer-documented workaround](https://vasp.at/forum/viewtopic.php?t=19372), LHFCALC=T with AEXX=0, initializes ISYM=3 machinery while retaining zero exact exchange. It was applied only to a new private PBE control, never to production HSE inputs.

The corrected beta control directly passed VASP 5.4.4 to 6.6.1 WAVECAR reading, no ISTART=0 fallback, resolved NKPTS=69/NBANDS=40, and EDIFF convergence in 17 iterations within 73.937363 seconds of the 180-second cap. Resolved PBE coefficients were all one with AEXX=0. Its final -28.75501276 eV/f.u. differs from the PBE+D3(BJ) reference by 2e-8 eV/f.u., below the existing 1e-5 tolerance. The frozen launcher labels scientific convergence NOT ASSESSED for this LHFCALC=true control; independent raw OUTCAR/XML evidence establishes convergence. Stress was not evaluated (resolved ISIF=0). The [bootstrap audit](05_00_PBE_WAVECAR_BOOTSTRAP_AUDIT.md) records the complete runtime and integrity evidence.

The earlier ISYM=0 band-branch mismatch remains historical provenance: its k-to-minus-k-only reduction expanded the regular sets to 36/292/292/292/292 rather than donor counts 10/69/69/69/57. That was stale downstream policy, not a failed bootstrap. ISYM=3 preserves the approved sampling orbits without regenerating any path. All KPOINTS, KPOINTS_OPT, 20 points per segment and source structures are unchanged. NCORE=1 remains execution-layer launcher policy.

Only beta was directly runtime-tested under ISYM=3; spinel, alpha1, IIa_prime and IIb remain statically compatible. No HSE convergence, HSE timing advantage, or future HSE-to-KPOINTS_OPT restart acceptance is claimed.

## Final Selection

| Phase | KPOINTS regular mesh | KPOINTS_OPT path type | Points/segment | Verdict |
| --- | --- | --- | --- | --- |
| spinel | 4 x 4 x 4 | cF2 | 20 | PASS: VASP 6.6.1 KPOINTS representation |
| alpha1 | 12 x 12 x 4 | hR1 | 20 | PASS: VASP 6.6.1 KPOINTS representation |
| beta | 12 x 12 x 4 | hP2 | 20 | PASS: VASP 6.6.1 KPOINTS representation |
| IIa_prime | 12 x 12 x 4 | hP2 | 20 | PASS: VASP 6.6.1 KPOINTS representation |
| IIb | 12 x 12 x 4 | hP2 | 20 | PASS: VASP 6.6.1 KPOINTS representation |

All five input pairs pass the representation checks: automatic unshifted Gamma mesh; Line-mode/Reciprocal KPOINTS_OPT; exactly 20 points per segment; unchanged source endpoint coordinates, labels, and branch breaks. The original ten-file KPOINTS migration left INCARs unchanged and ran no VASP calculation. The subsequent restart alignment changes only HSE INCAR policy and these provenance documents, with a separate private beta PBE restart control; it does not change POSCAR, POTCAR, KPOINTS, KPOINTS_OPT, Stage 03/04, or any Stage 05-00 file.

## Caveats

- PBE k-point convergence alone does not prove HSE regular-mesh convergence.
- Twenty points per segment is not independently HSE path converged.
- Input selection and the reduced preflight do not establish production HSE runtime or memory requirements; no measured production speedup is claimed.
- HFRCUT=-1 is now explicit and consistent in both HSE branches; its production-mesh results still require HSE validation.
- KPOINTS_OPT execution requires the native VASP 6.6.1 runtime and appropriate parallel settings, including the validated local NCORE=1 constraint. The supplied legacy submission scripts do not establish that runtime.
