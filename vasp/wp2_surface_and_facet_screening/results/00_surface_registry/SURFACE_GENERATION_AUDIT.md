# WP2 Stage 0 — Surface generation audit

**12 September 2026. Stage 0 structural preparation is complete for human review.**
Eight unrelaxed seeds are retained, including five HOLD models. This is not
scientific approval for production, a reaction shortlist, or a Stage 01 launch.
No thickness or vacuum convergence is claimed.

## Parent structures

The authority is the [WP1 final geometry audit](../../../wp1_polymorph_polytype_benchmark/results/01_geometry_optimisation/FINAL_GEOMETRY_OPT_AUDIT.md),
read with the [WP1 closure summary](../../../wp1_polymorph_polytype_benchmark/results/wp1_summary.md)
and the existing [WP2 roadmap](../../ROADMAP.md). All paths below are repository-relative.

| Phase | Exact accepted source | Identity | Atoms / formula units | Verified network |
| --- | --- | --- | --- | --- |
| beta | `vasp/wp1_polymorph_polytype_benchmark/calculation/01_geometry_optimisation/beta/CONTCAR` | P3m1 (156) | 7 / 1 | One periodic ZnIn2S4 septuple sheet |
| IIb | `vasp/wp1_polymorph_polytype_benchmark/calculation/01_geometry_optimisation/IIb/CONTCAR` | P6_3mc (186) | 14 / 2 | Two periodic ZnIn2S4 septuple sheets |
| spinel | `vasp/wp1_polymorph_polytype_benchmark/calculation/01_geometry_optimisation/spinel/CONTCAR` | Fd-3m (227) | 56 / 8 | One three-dimensional connected framework |

Each source was copied byte-for-byte to `structure/bulk_parents/<phase>/POSCAR`.
Its adjacent metadata records source/destination SHA-256, composition, lattice,
accepted symmetry, coordination and provenance. Parent hashes are pinned in the
generator. No reduction, conventionalization, symmetry refinement, atom
reordering, lattice change or relaxation was applied to the parents.

For each phase, species/order, lattice and fractional coordinates independently
match `vasp/wp1_polymorph_polytype_benchmark/calculation/03_static_scf/<phase>/POSCAR`
within 1e-8. The [Stage 03 audit](../../../wp1_polymorph_polytype_benchmark/results/03_static_scf/STATIC_SCF_ANALYSIS.md)
records this same accepted geometry handoff. Stage 03 was a verification reference,
not the slab source. No CIF, Stage 04 copy, or HSE input supplied a parent.

## Generation protocol

Implementation: [generator](../../script/wp2/generate_surface_candidates.py),
[geometry utilities](../../script/wp2/surface_model_utils.py), and
[focused tests](../../script/wp2/test/test_surface_candidates.py).
The existing `.venv-wp1-plots` environment supplies pymatgen 2026.5.4,
NumPy 2.5.3, SciPy 1.18.1 and spglib 2.7.0. No dependencies were installed.

The maintained [pymatgen SlabGenerator API](https://pymatgen.org/pymatgen.core.html#pymatgen.core.surface.SlabGenerator)
constructs the oriented cell and slabs. Its installed implementation was also
inspected. Raw terminations are the increasing wrapped fractional shifts from
`gen_possible_terminations(ftol=0.01)`, followed by `get_slab(shift=...)`.
The 0.01 A plane tolerance resolves spinel subplanes about 0.04 A apart; the
initially inspected 0.05 A grouping would merge these narrow cuts. Every final
raw cut receives `phase_hkl_tNN` before rejection or deduplication.

Only the seven requested facet families are generated. Settings are
`primitive=False`, `lll_reduce=False`, `center_slab=False`,
`reorient_lattice=False`, `max_normal_search=2`. Neither symmetry-forced atom
removal nor reconstruction/repair is used. A basal seed contains one entire
oriented parent repeat: one beta layer or two IIb layers. Other orientations use
a 12 A minimum **repeat-span** request rounded up to full oriented cells; this
is not a thickness-convergence result. The actual physical thickness is measured
separately and depends on the termination.

The inherited in-plane vectors are rigidly rotated with a parallel to x and
the normal parallel to z. The c vector is replaced by a normal-only vector,
and atomic coordinates receive only that rigid rotation, centering translation
and in-plane periodic wrapping. The external atom-free gap is exactly 20 A:
`vacuum = cell_normal_length - (z_max - z_min)`. Both exterior margins are 10 A.
The POSCAR remains three-dimensionally periodic with a vacuum-separated z image;
the physical bond graph has no bond across that vacuum. No Selective Dynamics
flags are generated.

Each raw row records the oriented-cell integer matrix, repeat count, original
and final lattices, rotation, translation, wrapping integers, parent-site
mapping, normal in the parent frame, and its sign relative to hkl. Some equivalent
facet representations use the negative hkl normal; all final normals point along
+z. This does not preserve the literal parent basis. Inverse transformation tests
recover each corresponding parent site modulo integer lattice translations.

### Coordination and connectivity

Per-site cation–S neighbor distances were inspected to 5 A. The first resolved
gap of at least 0.5 A after four or six neighbors identifies a candidate first
shell. Each phase/species cutoff is the midpoint between the largest included
and smallest excluded distances across its sites. The separation and expected
coordination are checked before use; ambiguity aborts generation.

| Phase | Zn–S included max / excluded min / cutoff (A) | In–S included max / excluded min / cutoff (A) |
| --- | --- | --- |
| beta | 2.53728 / 3.28291 / 2.91009 | 2.70922 / 4.35696 / 3.53309 |
| IIb | 2.54130 / 3.32461 / 2.93296 | 2.71172 / 3.85153 / 3.28162 |
| spinel | 2.36822 / 4.41440 / 3.39131 | 2.61489 / 4.52805 / 3.57147 |

These phase-specific criteria reproduce WP1: every Zn has four S neighbors;
layered In has four/six, and spinel In has six. Sulfur coordination is obtained
from the reciprocal cation–S graph. Lost coordination is referenced to the
actual inherited parent site, not an assumed coordination for all atoms.

Broken Zn–S/In–S estimates count cation-side missing neighbors over the entire
slab, **not bonds per individual face**. An independent periodic-bond crossing
count at each cut equals those deficits for all 67 candidates. Sulfur deficits
represent the other endpoints; they must not be added again to claim more
cleaved bonds. None of these quantities is a surface energy.

Component periodicity is obtained from the rank of bond-graph translation
cycles. Rank-0 atoms/molecular fragments are rejected. Rank-1 edge ribbons and
rank-2 sheets are legitimate inherited networks; IIb's vdW-separated components
are not automatically classified as invalid. Every retained slab has no finite
component or unexpected contact shorter than its parent minimum. A separate
1.8 A floor guards gross overlaps. All retained distances are at least 2.3067 A.

Face fingerprints describe the outer 3.6 A, including ordered species planes
and their coordination. All under-coordinated sites are assigned to the nearer
exterior face for counts. Face-window charge is explicitly partition-dependent.
Whole-slab symmetry must exchange top and bottom before faces are called
equivalent; outermost S composition alone is insufficient.

### Duplicate and polarity semantics

Duplicate screening checks composition, inherited in-plane metric, physical
thickness (1e-4 A), species-resolved layer sequence in either normal direction,
and species/coordination environments. A periodic StructureMatcher then checks
translation/crystal equivalence with `ltol=1e-5`, `stol=0.003`, `angle_tol=0.01`,
no scaling and no primitive reduction. The lowest cut-shift ID is canonical.
Independent matcher grouping agrees for beta (100), IIb (100), and both spinel
families. Tests also cover explicit normal reversal and periodic translation.

Equivalent slab-pairs can have inequivalent faces. Both exposed sides of every
materialized slab receive stable face rows. Equivalent lower faces explicitly
reference their upper face; inequivalent sides keep independent identities.

Formal electrostatics uses Zn2+, In3+, S2−; all 67 current raw cuts are
stoichiometric and nominally neutral. Plane charges and the charge-center dipole
`sum(q_i * (z_i - slab_midpoint))` are recorded. Triage rules are:

- **HIGH_CONCERN:** nonstoichiometric, nonzero nominal total charge, or absolute
  formal dipole/area above 0.1 e/A.
- **WARNING:** otherwise inequivalent faces or absolute dipole/area above 0.001 e/A.
- **LOW_CONCERN:** otherwise equivalent faces with a small formal dipole.

These explicit thresholds are Stage 0 triage choices, not validated physical
polar/nonpolar boundaries. Charged-plane sequence and termination fingerprints
provide supporting context, not independent proof. No electrostatic field,
electronic reconstruction, vacuum-potential slope or dipole-correction behavior
has been calculated. Those require later VASP work.

## Layered basal validation

| Parent | Valid canonical cut | Shift | Verified outer-S gap (A) | Complete layers / atoms | Physical thickness / vacuum / cell height (A) |
| --- | --- | ---: | ---: | --- | --- |
| beta | `beta_001_t07` | 0.984484965391 | 2.72857 | 1 / 7 | 9.41986 / 20 / 29.41986 |
| IIb | `IIb_001_t07` | 0.492362750329 | 2.77917 | 2 / 14 | 21.57260 / 20 / 41.57260 |

These boundaries cross **zero Zn–S and zero In–S first-shell bonds**. Every
retained component is a complete seven-atom ZnIn2S4 sheet with unchanged bulk
coordination. The gap sizes independently reproduce WP1's layer analysis.
Beta's other six cuts split its chemically meaningful layer and are REJECT.
IIb has two equivalent valid gap choices: `t14` duplicates `t07`; its other
12 raw cuts split a layer and are REJECT. No atom was deleted to force symmetry.

For both canonical basal seeds, the outward-to-inward upper fingerprint is
S–Zn–S, whereas the lower fingerprint is S–In–S. Both outermost S atoms have
three cation neighbors, but the underlying cation identities and distances differ.
The faces are **inequivalent**. Full depth/composition/coordination fingerprints
are retained in both registries and per-slab metadata.

Formal dipoles are −2.39158 eA for beta and −4.78870 eA for IIb (about −0.1836
and −0.3667 e/A). Both basal seeds are **HOLD / HIGH_CONCERN**, despite passing
geometric and complete-layer checks. They are retained for inspection, not
silently forced into symmetric or zero-dipole models.

IIb retains both sheets and the accepted stacking sequence. For later deterministic
2/4/6-layer construction, keep the same cut and extend beta by 2/4/6 oriented
repeats or IIb by 1/2/3 full oriented repeats. Site mappings and transformations
record the stacking; no thickness matrix or extra thicker basal seed was generated.

## Edge/prismatic candidates

Beta (100) has three unique neutral, connected 28-atom ribbon seeds from six
raw shifts. All have approximately 12.3152 A physical thickness and 47.115 A²
one-face area. `t01` and `t03` are HOLD for formal dipole densities of approximately
−0.2851 and +0.2851 e/A. `t02` is PASS / WARNING: its dipole is approximately zero,
but its faces are chemically inequivalent. This is a meaningful reactive edge,
not a detached-fragment artifact.

Beta (110) has one unique 56-atom seed from two shifts: `beta_110_t01`,
PASS / LOW_CONCERN, equivalent faces, approximately zero formal dipole,
13.5740 A physical thickness and 81.606 A² area. Its periodic ribbon is connected.
It breaks two Zn–S and six In–S bonds per slab-pair cell, compared with one
Zn–S and four In–S for beta (100) `t02`. These counts do not rank stability.

IIb (100) has two unique 56-atom seeds from six shifts. Both contain two expected
periodic ribbons, not finite fragments. `t01` is HOLD / HIGH_CONCERN, with an
approximately −0.2845 e/A formal dipole density; `t02` is PASS / LOW_CONCERN with
equivalent faces and approximately zero dipole. Their area is 94.566 A².

## Spinel feasibility gate

Both branches are **DEFERRED** in this bounded unreconstructed search. Neither
receives a materialized POSCAR or a reaction-face recommendation.

| Family | Raw / unique | Atoms per raw slab | Canonical seeds | Absolute formal dipole/area (e/A) | Decision |
| --- | --- | ---: | --- | --- | --- |
| (001) | 16 / 2 | 112 | `spinel_001_t01`, `spinel_001_t02` | 0.3767, 1.1302 | DEFERRED |
| (110) | 16 / 2 | 112 | `spinel_110_t01`, `spinel_110_t02` | 0.1884, 0.5651 | DEFERRED |

All are stoichiometric, nominally neutral and connected, with inequivalent
faces and nonzero formal dipoles. Narrow raw cuts (about 0.0401 A for one (001)
class and 0.0567 A for one (110) class) expose different charged-plane sequences;
they remain in the registry rather than disappearing under a broad plane tolerance.

These models require an unresolved electrostatic/compensation decision before
propagation. This is not proof that no viable spinel surface exists. Uniform
lateral replication cannot remove their dipole per area, so no unmodified 2×2
replicas were generated. No chemically justified compensation pattern was
established within this cleavage-only scope. Arbitrary atom removal, changed
electron counts, passivation, chemical-potential optimization and reconstruction
searches were not introduced merely to force a positive gate.

## Candidate inventory

| Phase / facet | Raw | Unique, all statuses | Duplicates | Materialized | PASS | HOLD | REJECT | DEFERRED |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| beta (001) | 7 | 7 | 0 | 1 | 0 | 1 | 6 | 0 |
| beta (100) | 6 | 3 | 3 | 3 | 1 | 5 | 0 | 0 |
| beta (110) | 2 | 1 | 1 | 1 | 1 | 1 | 0 | 0 |
| IIb (001) | 14 | 7 | 7 | 1 | 0 | 2 | 12 | 0 |
| IIb (100) | 6 | 2 | 4 | 2 | 1 | 5 | 0 | 0 |
| spinel (001) | 16 | 2 | 14 | 0 | 0 | 0 | 0 | 16 |
| spinel (110) | 16 | 2 | 14 | 0 | 0 | 0 | 0 | 16 |
| **Total** | **67** | **24** | **43** | **8** | **3** | **14** | **18** | **32** |

Counts of statuses refer to raw rows. Duplicate flags are orthogonal: a valid
duplicate becomes HOLD with `duplicate_of`, while a rejected/deferred duplicate
retains its scientific disposition. Thus 43 duplicates must not be added to the
status total. Of the eight materialized unique seeds, three are PASS and five
are HOLD. There are 16 face rows; explicit lower-to-upper equivalences occur for
beta (110) `t01` and IIb (100) `t02`. The other six retained slab-pairs retain
two inequivalent face identities each.

Deliverables: [surface_registry.csv](surface_registry.csv),
[face_registry.csv](face_registry.csv), three byte-preserved parent POSCARs with
metadata, and eight retained slab POSCARs with metadata. Every rejected,
duplicate or deferred raw structure can be regenerated from the pinned parents,
settings, shifts and generator. Metadata includes generation-code hashes and
software versions; source provenance uses repository-relative paths.

## Pilot recommendation

### Pilot A — beta basal

Recommend **`beta_001_t07`**, one complete septuple layer, seven atoms, area
13.026 A². It is the clean, bond-preserving minimal beta basal seed. Upper/lower
faces are inequivalent and the screening result is HIGH_CONCERN, so retain its
**HOLD** disposition for human electrostatics review before preparing a runtime
pilot. This is an appropriate minimal candidate for subsequent workflow testing;
it is not a bulk-like or thickness-converged slab.

### Pilot B — beta edge

Recommend **`beta_100_t02`**, 28 atoms, area 47.115 A², PASS / WARNING. It is
stoichiometric, neutral, connected, free of detached atoms/finite fragments, and
has approximately zero formal dipole despite inequivalent faces. The upper/lower
under-coordinated counts are Zn 1/0, In 1/1, S 2/2. One Zn–S and four In–S
bonds are broken in the slab-pair cell, providing an actual edge-SCF stress test.
It is smaller than the 56-atom beta (110) alternative and remains chemically
meaningful. No assertion of lower surface energy is involved.

These are recommendations for workflow validation, tractability and representative
electronic difficulty. They are not claims of experimental facet dominance,
lowest surface energy or greatest catalytic activity. Neither pilot was prepared
as a calculation input or launched.

## Validation and repository boundaries

Eight focused unittest cases passed. They cover pinned parent identity and the
Stage 03 handoff; deterministic IDs and byte-identical regeneration; registry
and duplicate/face references; every retained structure's serialization,
orientation, inverse parent mapping and actual vacuum; basal bond-cut rejection;
normal reversal without face-identity collapse; finite-fragment versus periodic
component detection; and no import/help writes. A pre-existing different parent
is refused before any output is replaced. spglib emitted deprecation warnings
about its legacy error handling; no test failed.

From the repository root, reproduce with:

```sh
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv-wp1-plots/bin/python vasp/wp2_surface_and_facet_screening/script/wp2/generate_surface_candidates.py
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv-wp1-plots/bin/python -m unittest discover -s vasp/wp2_surface_and_facet_screening/script/wp2 -p test_surface_candidates.py -v
```

Use `--output-root /private/tmp/<new-directory>` to compare regeneration without
overwriting reviewed outputs. The generator refuses different existing output
bytes and does not rewrite identical files. The audit is a human-review document,
not a calculated property or an automatically updated scientific verdict.

Starting branch: `main`. The pre-existing tree had three modified beta Stage
05-01 inputs (INCAR, KPOINTS, POSCAR), nine untracked beta Stage 05-01 runtime
files, the untracked WP1 summary, and the untracked WP2 scaffold. These WP1 files
were preserved. Only empty obsolete `structure/01_bulk_parents`,
`structure/02_slab_candidates`, and `structure/03_selected_surfaces` placeholders
were removed to implement the requested unnumbered architecture. Existing WP2
roadmap, calculation stages and later results directories were preserved.

No VASP or CMW job was launched; no calculation directory was created or modified;
no INCAR, KPOINTS or POTCAR was generated; no relaxation, convergence calculation,
surface energy, work function, PDOS, adsorption, reconstruction or HSE calculation
was performed. The active WP1 beta 05-02 Scratch directory was not accessed.
No files were staged, committed or pushed. **Stop at Stage 0.**
