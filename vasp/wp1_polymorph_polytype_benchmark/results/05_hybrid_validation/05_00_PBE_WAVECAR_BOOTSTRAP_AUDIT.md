# WP1 Stage 05-00 — PBE WAVECAR Bootstrap Audit

Audit date: 2026-09-07. Scope: spinel, alpha1, beta, IIa_prime, and IIb only.

## Executive Verdict

**All five Stage 05-00 PBE WAVECAR bootstrap calculations PASS. Stage 05-00 freeze decision: PASS. No bootstrap recalculation was required or performed.**

- Donor decisions: **5 PASS, 0 PASS WITH CAVEAT, 0 REQUIRES RECALCULATION**.
- The finalized production dependency is **05-00 PBE WAVECAR → 05-01 HSE06 SCF → converged 05-01 HSE WAVECAR → 05-02 HSE06 restart + KPOINTS_OPT**. Direct 05-00 → 05-02 is no longer the production dependency.
- All five 05-01 targets are statically compatible with the existing donors. A private beta PBE control directly demonstrated VASP 5.4.4 → VASP 6.6.1 WAVECAR acceptance at the unchanged production mesh, without fresh-start fallback.
- The original audit correctly found a stale downstream ISYM=0 mismatch and held the broader two-branch freeze. That finding is retained below as history; the downstream settings and dependency have now been corrected without altering or rerunning 05-00.
- **Global evidence boundary:** only beta cross-version restart was exercised at runtime. Neither full-mesh HSE convergence nor an eventual 05-01 HSE WAVECAR → 05-02 restart has been demonstrated. These are downstream gates, not failures of the completed bootstrap contract.
- The scientific gate for the two requested local restart-alignment/provenance commits is satisfied. No Stage 05-00 calculation commit is warranted because its reproducible inputs are already committed and unchanged.

## Purpose and Scope

**WP1 Stage 05-00 — PBE WAVECAR Bootstrap** generates converged PBE orbitals on the accepted frozen WP1 structures and intended regular meshes. It is not geometry optimisation, not HSE06, and not a replacement thermodynamic or energy-ranking stage.

The production dependency is:

```text
Stage 03 frozen POSCAR / charge-density lineage
    → 05-00 PBE WAVECAR bootstrap
    → 05-01 regular-mesh HSE06 SCF
    → converged HSE WAVECAR
    → 05-02 regular-mesh HSE06 restart + KPOINTS_OPT path
```

The initial audit was read-only except for this report. This follow-up modifies only the five 05-01 INCARs, five 05-02 INCARs, this report, and KPOINTS_SELECTION_JUSTIFICATION.md. No Stage 05-00 file, POSCAR, POTCAR, regular/path sampling, Stage 03/04 file, launcher, or tracked port file was modified. No bootstrap or production HSE calculation was run.

One explicitly authorized, separately staged private beta PBE restart-control calculation was run through the frozen generic launcher, with a 180-second limit. It used the unchanged 12×12×4 mesh and donor WAVECAR. It is a restart-format/use test, not a new scientific production dataset. Parsing and integrity tables were in memory; no scripts, plots, JSON/YAML/CSV, or additional scientific-repository artifacts were created.

## Git Reconciliation

The initial audit and the follow-up checked status, branch, recent history, working/staged diffs, and remotes. The follow-up fetch succeeded after filesystem approval for Git metadata access. HEAD and origin/main remained at the following starting values; the two existing unpushed commits were not amended or pushed.

| Item | Starting state after fetch |
|---|---|
| Branch | `main` |
| HEAD | `893b22db696f7b870a46288aeabc644de0114cd5` |
| origin/main | `84c45af2bbca60e8e8df2124cb492582908e86f7` |
| Ahead / behind | 2 / 0 |
| Tracked working-tree changes | None |
| Staged changes | None |
| Untracked files | 30 Stage 05-00 outputs: CONTCAR, DOSCAR, EIGENVAL, OSZICAR, OUTCAR, XDATCAR in each phase; this audit report was additionally untracked at follow-up start |
| Unrelated visible changes | None |
| Remote | `https://github.com/angzeli/znin2s4-computational.git` |

The existing local commits are b1284e4 (calc: migrate WP1 HSE bands to KPOINTS_OPT) and 893b22d (docs: update WP1 HSE KPOINTS justification). All 15 bootstrap input files (INCAR, KPOINTS, POSCAR in each phase) were already committed and unchanged. The local `.sub` files and heavy/licensed restart files remain ignored. Existing output tracking in earlier stages does not override this task's stricter exclusion list.

## Inventory and Attempt Provenance

Every phase contains exactly the following 17 files: INCAR, POSCAR, POTCAR, KPOINTS, WAVECAR, CHGCAR, OUTCAR, OSZICAR, vasprun.xml, CONTCAR, EIGENVAL, IBZKPT, DOSCAR, XDATCAR, PCDAT, CHG, REPORT; each also contains `vasp5.4.sub` (18 files including the submission script). No separate scheduler stdout/stderr or job-exit record is present. There are no duplicate attempt-output sets or symlinked files in these directories.

All required inputs, convergence records, and WAVECAR files exist and are nonempty. CHG and REPORT are zero bytes in all five directories; these are not required products of this static, `LCHARG=F` bootstrap. The nonzero CHGCAR files are retained upstream inputs, not claimed new charge-density outputs. Each XML is well formed and closed, each EIGENVAL has the complete declared set of finite band records, and each XDATCAR contains one configuration. There is no evidence of a truncated required output or a conflicting failed attempt in the available files.

| Phase | OUTCAR bytes | OSZICAR bytes | XML bytes | CHGCAR bytes | Embedded run start, 2026-09-04 | Local OUTCAR mtime, 2026-09-07 | Local WAVECAR mtime, 2026-09-07 |
|---|---:|---:|---:|---:|---|---|---|
| spinel | 219416 | 2413 | 173212 | 74617440 | 20:20:58 | 18:10:46 | 18:13:46 |
| alpha1 | 361592 | 1483 | 389536 | 30847098 | 20:26:31 | 18:05:28 | 18:10:10 |
| beta | 177505 | 1204 | 205605 | 10282526 | 20:25:25 | 18:05:54 | 18:06:28 |
| IIa_prime | 289216 | 1762 | 309431 | 20564812 | 20:16:39 | 18:08:46 | 18:11:04 |
| IIb | 242412 | 1297 | 256647 | 20564812 | 20:21:45 | 18:06:40 | 18:08:34 |

OUTCAR and XML consistently identify VASP `5.4.4.18Apr17-6-g9f103f2a35`, complex parallel build dated 2022-02-10, with eight cores. The `.sub` scripts also target VASP 5.4.4 and request 24 hours. These are not VASP 6.6.1-generated donors.

Local filesystem dates are later than the embedded execution dates, consistent with a copied result set; the transfer itself is not independently logged. Filesystem mtimes therefore are not execution timestamps. WAVECAR is not suspiciously older than the companion local outputs. Its lattice, dimensions, k points, and stored eigenvalues independently tie it to each completed bootstrap result rather than relying on dates alone.

## Input Lineage

### Frozen structures and PAW sequence

For every phase, the 05-00 POSCAR is **byte-identical** to the Stage 03 POSCAR and to the corresponding 05-01 and 05-02 POSCARs. This establishes identical lattice representation, species names, counts, fractional coordinates, and atom ordering without a geometry transformation or CIF round-trip.

The Stage 03 and 05-00 CONTCARs retain the same cell and atom ordering. Comparing each 05-00 CONTCAR to its POSCAR gives zero lattice-vector difference and a maximum fractional-coordinate difference of at most `3.56e-15`. XML initial/final structures agree within their printed precision. XML's additional symmetry-derived primitive-cell record is not an ionic or cell update and was not substituted for the actual cell.

| Phase | POSCAR species order | Counts in that order | Total atoms | PAW order |
|---|---|---|---:|---|
| spinel | S, Zn, In | 32, 8, 16 | 56 | S, Zn, In_d |
| alpha1 | Zn, In, S | 3, 6, 12 | 21 | Zn, In_d, S |
| beta | S, In, Zn | 4, 2, 1 | 7 | S, In_d, Zn |
| IIa_prime | In, S, Zn | 4, 8, 2 | 14 | In_d, S, Zn |
| IIb | In, Zn, S | 4, 2, 8 | 14 | In_d, Zn, S |

Safe POTCAR title metadata identifies PAW_PBE Zn, In_d, and S datasets dated 06Sep2000. Within each phase, SHA-256 equality establishes identical POTCAR content across Stage 03, 05-00, 05-01, and 05-02. The actual OUTCAR PAW sequence and ion counts agree. No PAW dataset content is reproduced here or committed.

Each bootstrap CHGCAR is SHA-256-identical to its corresponding Stage 03 CHGCAR. Its structure header is printed to six decimal places: maximum component differences relative to POSCAR are `4.72e-7` Å for the cell and `4.96e-7` in fractional coordinates across the set. These are output-format rounding, not a changed structure; they exceed the reusable launcher's strict `1e-7` structure threshold if that threshold is applied blindly to CHGCAR text. That comparison does not invalidate the completed, converged `ICHARG=1` calculations or their full-precision WAVECAR headers. This audit does not certify that a future CHGCAR-only launch would pass that separate launcher preflight unchanged.

### Regular sampling

The actual KPOINTS files in all four compared stages specify automatic Gamma-centred meshes with zero shift. The signatures match phase by phase. Bootstrap IBZKPT weights sum to the full mesh size.

| Phase | Mesh | Shift | Bootstrap nkpts | Sum of IBZ weights | 05-01 mesh match | 05-02 mesh dimensions match |
|---|---|---|---:|---:|---|---|
| spinel | 4 × 4 × 4 | 0, 0, 0 | 10 | 64 | Yes | Yes; aligned ISYM=2 |
| alpha1 | 12 × 12 × 4 | 0, 0, 0 | 69 | 576 | Yes | Yes; aligned ISYM=2 |
| beta | 12 × 12 × 4 | 0, 0, 0 | 69 | 576 | Yes | Yes; aligned ISYM=2 |
| IIa_prime | 12 × 12 × 4 | 0, 0, 0 | 69 | 576 | Yes | Yes; aligned ISYM=2 |
| IIb | 12 × 12 × 4 | 0, 0, 0 | 57 | 576 | Yes | Yes; aligned ISYM=2 |

## Resolved Bootstrap Settings

All five INCARs have the same scientific settings. Resolved values were checked against OUTCAR and XML; absent explicit controls are distinguished from resolved defaults.

| Control | Actual setting / resolved evidence |
|---|---|
| Functional | `GGA=PE`; XML `LHFCALC=F`, no hybrid calculation |
| Dispersion | `IVDW=12`, also recorded in OUTCAR |
| Basis | `ENCUT=500 eV`; `PREC=Accurate` |
| Convergence | `EDIFF=1e-7 eV`; `NELM=100` |
| Electronic algorithm | `ALGO=Normal`; resolved `IALGO=38` |
| Spin / occupations | `ISPIN=1`; `ISMEAR=-5`; SIGMA absent from INCAR, resolved `0.20 eV` |
| Startup | `ISTART=0`, `ICHARG=1`: new orbitals with retained charge-density initialization |
| Geometry | `IBRION=-1`, `NSW=0`; ISIF absent from INCAR, resolved `2` |
| Projectors / PAW / grid | `LREAL=F`, `LASPH=T`, `ADDGRID=T` |
| Symmetry | `ISYM=2` |
| Restart writing | `LWAVE=T`; `LCHARG=F` |
| NBANDS | Not explicitly set; resolved phase-specific values below |
| Parallelism | NCORE/KPAR/NPAR absent from INCAR; output records eight cores, one k-point group, one core per band, XML `NPAR=8` |

This is a **PBE electronic bootstrap on the existing PBE+D3(BJ) workflow geometry**. IVDW=12 adds the geometry-dependent D3(BJ) correction to energies and associated derivatives; it is not a hybrid-like orbital-Hamiltonian or band-gap correction. [VASP DFT-D3 documentation](https://vasp.at/wiki/DFT-D3).

## Calculation Completion

All five OUTCARs explicitly state that the electronic loop stopped because EDIFF was reached. OSZICAR iteration counts agree with both OUTCAR and the XML scstep counts. Each run has one static calculation, a complete timing footer, and closed XML; the verdict is not inferred from a footer or exit code alone. No scheduler exit code is available.

| Phase | Mesh | Electronic iterations | Converged? | Normal termination | Elapsed seconds | Verdict |
|---|---|---:|---|---|---:|---|
| spinel | 4 × 4 × 4 | 25 / 100 | Yes, EDIFF | Yes | 3967.733 | PASS |
| alpha1 | 12 × 12 × 4 | 15 / 100 | Yes, EDIFF | Yes | 720.660 | PASS |
| beta | 12 × 12 × 4 | 12 / 100 | Yes, EDIFF | Yes | 48.002 | PASS |
| IIa_prime | 12 × 12 × 4 | 18 / 100 | Yes, EDIFF | Yes | 326.940 | PASS |
| IIb | 12 × 12 × 4 | 13 / 100 | Yes, EDIFF | Yes | 200.903 | PASS |

The final Davidson energy changes and band-energy changes are both below the requested `1e-7 eV` threshold. No NELM exhaustion, fatal numerical/MPI error, NaN/Inf token, or evidence of walltime termination was found in OUTCAR, OSZICAR, or XML. Routine PAW reference-energy messages and performance suggestions were not misclassified as failures. Printed forces and stresses are finite; residual force magnitudes are not a new geometry-acceptance test for this frozen-geometry stage.

| Phase | Final OUTCAR TOTEN, eV/cell | Maximum absolute force component, eV/Å | Maximum absolute stress component, kbar |
|---|---:|---:|---:|
| spinel | -233.74241460 | 0.00031694 | 0.07164287 |
| alpha1 | -86.24598338 | 0.00577614 | 0.05665535 |
| beta | -28.75501278 | 0.00816703 | 0.28954956 |
| IIa_prime | -57.45083657 | 0.01043088 | 0.12354651 |
| IIb | -57.50299153 | 0.00618595 | 0.30126065 |

These final totals include the D3(BJ) contribution and agree with XML `e_fr_energy` and rounded OSZICAR F values. The final XML calculation-level `e_0_energy` is a zero-valued field in these files and was **not** used as a physical zero energy; the final electronic scstep energies precede the dispersion contribution. This table records provenance only, with no cross-phase ranking or replacement of Stage 03 energies.

## WAVECAR Validation

For the header-only checks, the existing `scripts/run_vasp.py::wave_header` logic from the local VASP 6.6.1 macOS ARM64 port was reused in memory with its structure/INCAR helpers. Those checks did not invoke a launcher entry point or execution/staging routine. The separately authorized dry-run and private runtime control below did use the public launcher CLI; no tracked port file was modified.

All five WAVECARs have one spin channel and little-endian marker **45200** (complex single-precision coefficient format). Header cutoff is 500 eV, matching OUTCAR and the inputs. Header lattice vectors pass the helper's `1e-7 Å` comparison against the exact POSCAR cell.

| Phase | WAVECAR bytes (MiB) | Header status | Record bytes | nkpts | NBANDS | ENCUT, eV | Static downstream compatibility |
|---|---:|---|---:|---:|---:|---:|---|
| spinel | 745957632 (711.401) | PASS | 244416 | 10 | 304 | 500 | 01: STATICALLY COMPATIBLE; sequential 02 targets aligned |
| alpha1 | 760184128 (724.968) | PASS | 97472 | 69 | 112 | 500 | 01: STATICALLY COMPATIBLE; sequential 02 targets aligned |
| beta | 91679104 (87.432) | PASS | 32384 | 69 | 40 | 500 | 01: STATICALLY COMPATIBLE; sequential 02 targets aligned |
| IIa_prime | 369006000 (351.912) | PASS | 66000 | 69 | 80 | 500 | 01: STATICALLY COMPATIBLE; sequential 02 targets aligned |
| IIb | 300345856 (286.432) | PASS | 65024 | 57 | 80 | 500 | 01: STATICALLY COMPATIBLE; sequential 02 targets aligned |

For every file, the measured length is **exactly** `(2 + nspin × nkpts × (NBANDS + 1)) × record_length`. All k-point metadata records are readable and finite. Their positive plane-wave counts fit the declared coefficient record length; ranges are 30280–30551, 12005–12183, 3974–4047, 8098–8250, and 8007–8127 respectively in table order. The sizes are therefore plausible and no missing/truncated record is indicated.

Header nkpts and NBANDS agree with OUTCAR, XML, EIGENVAL, and IBZKPT as applicable. The maximum fractional k-coordinate residual against IBZKPT is below `4e-15`; against lower-precision XML it is below `3.34e-9`. WAVECAR band-eigenvalue metadata agree with EIGENVAL to below `5e-7 eV` and XML to below `5e-5 eV`, consistent with those text formats. This links each donor to the converged electronic result independently of transfer timestamps.

Only headers and per-k-point band metadata were read, **not wavefunction coefficient arrays**. These checks establish basic integrity and identity, not every coefficient's correctness or arbitrary cross-version compatibility. The separate beta runtime control below supplies the directly exercised restart evidence.

## Downstream HSE Initialization

### Production INCAR policy and evidence

Both HSE branches now use LHFCALC=T, GGA=PE, HFSCREEN=0.2, ALGO=Normal, LFOCKACE=T, HFRCUT=-1, ISYM=2, ISPIN=1, ENCUT=500 eV, ISTART=1, ICHARG=0, IBRION=-1, and NSW=0. The HSE exchange fraction remains the unchanged default; no AEXX override was introduced.

The official [ACE documentation](https://vasp.at/wiki/LFOCKACE) supports Davidson/ALGO=Normal; ACE is not active under Damped/All. The port's VALIDATION_REPORT.md and docs/VALIDATION.md independently record bounded HSE/ACE and reduced-beta KPOINTS_OPT runtime evidence. ALGO=Damped was therefore replaced and the obsolete TIME=0.4 removed in both branches. This selects an already tested algorithm, not a claim of production speedup.

The official [hybrid-band guidance](https://vasp.at/wiki/Band-structure_calculation_using_hybrid_functionals) recommends HFRCUT=-1 for gapped band workflows; the [HFRCUT documentation](https://vasp.at/wiki/HFRCUT) defines its automatic cutoff treatment. Both HSE branches adopt it to keep the finite-mesh Coulomb treatment consistent. This changes the previous default numerical treatment, not the intended HSE06 fraction/screening. No band gap or finite-mesh energy convergence is claimed by setting the tag.

ISYM=2 is a regular-grid restart-consistency choice, not a KPOINTS_OPT requirement. The optional high-symmetry path remains exactly unchanged. Convergence controls, smearing, projectors, PAW options, fixed geometry and output flags were otherwise preserved: 05-01 retains NELMIN=6 and LWAVE=T/LCHARG=T; 05-02 retains NELMIN=8 and LWAVE=F/LCHARG=F. Both retain NELM=150 and EDIFF=1e-7 eV.

### 05-00 → 05-01

All five targets have byte-identical donor POSCAR/POTCAR, matching regular mesh signatures, spin, cutoff and ISYM, and explicit donor-matching NBANDS:

| Phase | Mesh | Expected ISYM=2 regular nkpts | Explicit NBANDS in 01 and 02 | Cross-version evidence |
|---|---|---:|---:|---|
| spinel | 4×4×4 | 10 | 304 | STATICALLY COMPATIBLE FOR INTENDED 05-01 START |
| alpha1 | 12×12×4 | 69 | 112 | STATICALLY COMPATIBLE FOR INTENDED 05-01 START |
| beta | 12×12×4 | 69 | 40 | Static match; separate PBE runtime acceptance PASS |
| IIa_prime | 12×12×4 | 69 | 80 | STATICALLY COMPATIBLE FOR INTENDED 05-01 START |
| IIb | 12×12×4 | 57 | 80 | STATICALLY COMPATIBLE FOR INTENDED 05-01 START |

The expected counts are inherited from the verified donor symmetry sets, not claimed actual future HSE output values. Stored marker 45200, spin and header records are valid in all five files. Counts above NELECT/2 are respectively 56, 19, 9, 18 and 18; NBANDS is accepted for the current band-edge validation, not independently converged high-energy conduction-band science. Explicit counts eliminate reliance on automatic band selection, but the launched parallel layout must still be checked for any padding.

PBE orbitals are a permitted initial guess for HSE; functional equality is not a restart requirement. ICHARG=0 constructs the initial density from the orbitals and permits subsequent self-consistency. [VASP ICHARG documentation](https://vasp.at/wiki/ICHARG).

### 05-01 → 05-02

The actual 05-01 HSE WAVECAR does not yet exist. Its future 05-02 target now matches exact structure and PAW ordering, automatic regular mesh, ISYM=2, explicit NBANDS, ISPIN, cutoff, HSE06/ACE settings, PRECFOCK and HFRCUT. KPOINTS_OPT is the only added path sampling, not a replacement regular mesh. Existing branch-specific convergence and output settings remain intentional.

After 05-01 converges, validate its HSE donor identity and header, then use that donor for 05-02 through the frozen launcher's public CLI. Actual restart acceptance and optional-path convergence must be checked at that later execution. KPOINTS_OPT evaluates path states after regular SCF; the chosen sequential stages provide a provenance boundary and reusable HSE donor, not a claimed runtime advantage. [Official KPOINTS_OPT workflow](https://vasp.at/wiki/KPOINTS_OPT).

### Historical mismatch finding — resolved, not erased

The initial audit found direct reuse of ISYM=2 bootstrap WAVECARs in then-ISYM=0 automatic-mesh band inputs incompatible: stored counts were 10/69/69/69/57, whereas k↔−k-only grid reduction gives 36/292/292/292/292. The latter were analytic counts, not measurements from an unrun VASP job. Inappropriate restart k-point counts can cause fresh-start fallback. [VASP ISTART documentation](https://vasp.at/wiki/ISTART).

This was stale downstream configuration following the representation migration, not donor corruption. The older VASP 5.4.4 explicitly listed weighted-plus-zero-weight method and direct PBE starting route were valid for their implementation. The new production route now takes 05-02 orbitals from converged 05-01 HSE, with matching ISYM=2 regular sampling. No bootstrap rerun was needed; the historical mismatch is no longer a freeze blocker.

### Frozen-launcher beta HSE dry-run

The public scripts/run-vasp.sh CLI was invoked on the actual updated 05-01/beta input, with --restart wavecar and --restart-from pointing to 05-00/beta, --binary std, --ranks 8, --ncore 1, --kpar 1, --mpi-mode synthetic, and --dry-run. The hypothetical output was private/wp1-beta-hse-dry-20260907-202446 under the external port; it remained absent.

**PASS, exit 0, no files created and no MPI/VASP execution.** The plan selected the recorded 6.6.1 std binary, 8 synthetic-MPI ranks, NCORE=1/KPAR=1, ISTART=1/ICHARG=0, the 69-point/40-band donor, and KPOINTS_OPT inactive. There was no compatibility refusal. The CLI does not print every scientific tag: direct inspection of the unchanged input consumed by that plan additionally verified Gamma 12×12×4, ISYM=2, NBANDS=40, HSE06, ALGO=Normal, ACE=T and HFRCUT=-1. These are planned settings, not resolved HSE runtime observations.

NCORE/KPAR/NPAR remain absent from all ten canonical INCARs. The frozen launcher injects the requested parallel tags into execution copies only. No launcher or tracked port file was changed; legacy ignored VASP-5 submission scripts are not the native execution route.

### Private beta cross-version PBE acceptance control

**PASS.** One new private execution was staged through the same public CLI, with --timeout 180 and the same std/8-rank/NCORE=1/KPAR=1/synthetic settings. Its input copied only the unchanged beta POSCAR, POTCAR and KPOINTS. The small PBE INCAR retained GGA=PE, IVDW=12, ENCUT=500, PREC=Accurate, EDIFF=1e-7, ALGO=Normal, ISPIN=1, ISMEAR=-5, fixed geometry, LREAL=F, LASPH=T and ADDGRID=T; it used NELM=40, ISTART=1, ICHARG=0, ISYM=2, NBANDS=40, LWAVE=F and LCHARG=F. There were no HSE tags or KPOINTS_OPT.

| Check | Measured evidence |
|---|---|
| Runtime version | vasp.6.6.1, build 2026-09-07 10:40:42, complex |
| Binary SHA-256 | 4d147f16d515b1c6000b1053cfe5b234a63fb29789bb051346dc9019d59567db |
| Frozen launcher repository | 4d67212dde4148f3aa4085da8aa1867b0bc967ef; clean tracked state |
| Restart acknowledgement | stdout: “the WAVECAR file was read successfully” |
| Effective startup | ISTART=1, ICHARG=0; no fresh-start fallback |
| Regular state dimensions | nkpts=69, NBANDS=40, ENCUT=500 eV, ISPIN=1, ISYM=2, 7 atoms / 62 electrons |
| SCF completion | EDIFF reached in 12 of 40 iterations; normal footer and closed XML |
| Integrity | Finite energy/forces/stress; no NaN/Inf or fatal error; empty stderr |
| Launcher result | child_status=0, reason=completed, scientific_convergence=CONFIRMED (static electronic), restart_observed=WAVECAR ACCEPTED |
| VASP / launcher elapsed | 82.913 s / 83.867406 s; below the original 180 s cap, no extension |
| Reference beta bootstrap TOTEN | -28.75501278 eV/f.u. |
| Restart-control TOTEN | -28.75501263 eV/f.u. |
| Absolute energy difference | 1.5e-7 eV/f.u.; PASS against predeclared 1e-5 eV/f.u. tolerance |

The seven-atom beta cell contains one formula unit. The converged comparison shares the structure, PAW, PBE+D3(BJ), mesh, cutoff and occupation scheme; restart initialization, band explicitness and output controls intentionally differ. Energy agreement was judged only after convergence, not from an assumption of identical startup trajectories.

Local donor SHA-256 was identical before and after the test. Private logs and RUN_METADATA.txt remain in the external port at private/wp1-stage05-restart-validation/20260907-202446/pbe-acceptance; the corresponding input directory retains the small control inputs. The duplicate staged WAVECAR was removed after validation to avoid retaining a large temporary restart copy; the original scientific donor remains untouched. PAW/restart fingerprints are not reproduced in the public report.

This is directly exercised beta VASP-5-to-VASP-6 PBE restart acceptance, not an HSE run or a runtime test of the other four phases. No production HSE convergence, production-mesh ACE speedup, or future HSE-to-band acceptance follows from this control.

## Relationship to Stage 03

Stage 03 remains the authoritative **PBE+D3(BJ) static-energy dataset**. All five Stage 03 INCARs and resolved OUTCARs have LWAVE=F. Stage 05-00 adds the missing PBE orbital files, with LWAVE=T, an explicit bootstrap SYSTEM title, distinct output contents and execution records, and donor eigenvalues matching its own completed result.

The retained CHGCAR lineage is Stage 03, while the newly produced WAVECAR is Stage 05-00. No geometry optimization occurred. No Stage 03 file, energy, or scientific conclusion was replaced or recalculated.

## Phase-Level Acceptance

| Phase | Input lineage | Original bootstrap SCF | WAVECAR | Intended 05-01 start | Decision |
|---|---|---|---|---|---|
| spinel | PASS | PASS, 25 steps | PASS | Statically compatible | PASS |
| alpha1 | PASS | PASS, 15 steps | PASS | Statically compatible | PASS |
| beta | PASS | PASS, 12 steps | PASS | Statically compatible; cross-version PBE control passed | PASS |
| IIa_prime | PASS | PASS, 18 steps | PASS | Statically compatible | PASS |
| IIb | PASS | PASS, 13 steps | PASS | Statically compatible | PASS |

All five decisions accept the donor contract: produce valid PBE starting orbitals for 05-01. No downstream HSE success is presumed. Only beta was runtime-tested across versions; the shared format and matched static contracts support the other four starts without four ceremonial reruns.

## Freeze Decision

**Stage 05-00 freeze decision: PASS.**

All five original donor calculations are accepted and remain unchanged. The sequential production inputs, beta HSE dry-run, and beta bounded cross-version PBE restart control satisfy the requested migration gate. The old 05-02 mismatch is resolved and is not an outstanding bootstrap blocker.

The two authorized local commits are scoped to: (1) the ten changed 05-01/05-02 INCARs, with message calc: align WP1 HSE restart workflow for VASP 6.6.1; (2) only this audit and KPOINTS_SELECTION_JUSTIFICATION.md, with message docs: finalize WP1 HSE restart provenance. Commit identities and exact file scopes are recorded by Git and the completion response, not self-referentially embedded in this document.

No Stage 05-00 calculation commit, empty commit, or output-only substitute is needed. The original 30 untracked bootstrap outputs remain untouched and untracked. No WAVECAR, POTCAR, CHGCAR, private validation output, binary or licensed PAW data enters either commit. No push is authorized or performed.

## Caveats

1. Cross-version restart was directly exercised on **beta only**, under matching PBE conditions. The other four donors are statically compatible, not individually runtime-certified.
2. The first real 05-01 HSE calculation still needs an explicitly authorized production-mesh execution budget, correct donor selection, and observed VASP 6.6.1 restart acceptance/resolved settings. Dense-mesh HSE convergence and practical cost remain unvalidated; the port's earlier dense HSE caps do not establish iteration timing or speedup.
3. Only after 05-01 converges and retains a verified HSE WAVECAR can actual 05-01 → 05-02 restart acceptance and KPOINTS_OPT convergence be tested. This future gate does not prevent freezing the PBE bootstrap.
4. Explicit NBANDS supports current band-edge validation but is not a high-energy conduction-band convergence study. Confirm the actual launch does not pad counts incompatibly under another parallel layout.
5. Original local mtimes reflect a later copied result set, not embedded run dates; no original scheduler/transfer log is available. The retained CHGCAR's six-decimal structure rounding can exceed a separate strict charge-restart preflight tolerance, but it is not a donor geometry change or a failure of the tested WAVECAR route.
