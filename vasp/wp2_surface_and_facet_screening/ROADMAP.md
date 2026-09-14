# 5. WP2 — Validated ZnIn2S4 Surfaces for Reaction Chemistry

## Objective

Identify a small, defensible set of ZnIn2S4 surfaces for the subsequent H2O2 / 2e− ORR and CO2RR / HER calculations.

WP2 is a **surface-model validation and screening work package**, not yet an adsorption study.

The central questions are:

1. Which ZnIn2S4 surfaces can be represented by chemically meaningful, numerically converged slab models?
2. How strongly do surface structure and termination modify:
   - energetic accessibility,
   - surface electronic states,
   - vacuum-referenced band-edge behaviour,
   - work function,
   - band-edge localisation,
   - local coordination?
3. Which 2–3 reaction-facing surfaces are sufficiently distinct and computationally tractable to justify propagation into WP3 and WP4?

Bulk band structures alone are insufficient because O2 reduction and CO2 reduction occur at surfaces.

Experimental surface relevance remains provisional until supported by thin-film XRD, preferred-orientation, microscopy, or related structural evidence.

---

## 5.1 Parent-Phase Scope

WP1 defines the starting phase hierarchy.

### Primary layered reference

Use **beta-ZnIn2S4** as the primary layered reference.

### Layered stacking control

Retain **IIb-ZnIn2S4** as a provisional stacking-dependent control because it is energetically close to beta and structurally distinct.

Do not automatically propagate every WP1 polytype into WP2.

### Three-dimensional contrast

Treat **spinel ZnIn2S4** only as a bounded structural contrast.

Its low calculated bulk energy must not by itself be interpreted as evidence that it is the experimentally relevant ambient-pressure thin-film phase.

Spinel surface calculations proceed only if chemically defensible stoichiometric or compensated slab models can be constructed without requiring a large reconstruction or chemical-potential search.

If this gate fails, defer spinel rather than forcing a nominal three-surface shortlist.

---

## 5.2 Initial Candidate Surface Families

Begin with the following candidate facet families:

| Parent | Initial cuts | Purpose |
| --- | --- | --- |
| beta | (001) | layered basal reference |
| beta | (100) | prismatic / broken-bond candidate |
| beta | (110) | structurally distinct prismatic candidate |
| IIb | (001) | matched layered basal control |
| IIb | (100) | stacking-dependent edge control |
| spinel | (001) | three-dimensional low-index contrast |
| spinel | (110) | alternative three-dimensional contrast |

These are candidate **facet families**, not necessarily seven final slab models.

Different termination shifts, exposed faces, and slab thicknesses must retain separate identities until demonstrated to be equivalent.

Do not enforce a target of 6–8 production surfaces if the structural or computational evidence does not justify them.

---

## 5.3 Slab-Construction Rules

All slabs must be derived from accepted WP1 relaxed bulk geometries.

For every generated slab, record:

- parent phase;
- Miller orientation;
- bulk-to-slab transformation;
- in-plane lattice vectors;
- surface-normal direction;
- lateral supercell;
- termination identity;
- upper-face identity;
- lower-face identity;
- composition;
- number of ZnIn2S4 formula units;
- physical slab thickness;
- actual atom-free vacuum separation;
- fixed and relaxed atoms;
- polarity-screen result.

### Layered basal surfaces

For basal slabs:

- cut only through verified interlayer / van der Waals gaps;
- preserve complete chemically meaningful layer units;
- do not silently remove atoms to manufacture symmetric surfaces;
- preserve the actual parent stacking sequence.

### Prismatic / broken-bond surfaces

For edge surfaces:

- enumerate chemically distinct termination shifts;
- reject duplicate structures after explicit structural comparison;
- identify under-coordinated Zn, In, and S sites;
- retain any physically meaningful surface reconstruction during relaxation.

### Inequivalent slab faces

Do not assume the two slab faces are equivalent.

A single asymmetric slab may represent two distinct exposed surfaces.

Surface identity must therefore be defined by:

**phase + facet + termination + exposed face**

rather than by Miller index alone.

### Spinel feasibility gate

For spinel (001)/(110), inspect stoichiometric, neutral and chemically defensible constructions, allowing only bounded lateral reconstruction tests where necessary.

Initial feasibility work may consider up to approximately 2 × 2 lateral cells.

Do not use arbitrary atom deletion, artificial electron-count modification, pseudo-hydrogen termination, or uncontrolled nonstoichiometric compensation merely to obtain a convenient slab.

If meaningful spinel slabs require a substantial surface-thermodynamics or reconstruction campaign, classify spinel as:

**DEFERRED**

for the initial WP2 reaction-surface shortlist.

---

## 5.4 Stage 01 — Runtime and Workflow Pilots

Before launching a large slab matrix, establish the computational cost and stability of the local workflow.

### Pilot A — beta basal slab

Use one complete beta layered unit as a **computational pilot only**, not as a converged bulk-surface model.

Initial local execution:

- VASP 6.6.1;
- PBE+D3(BJ);
- 500 eV ENCUT;
- two MPI ranks;
- one numerical-library thread per rank;
- one job at a time;
- approximately one-hour total pilot budget.

Run:

1. a fixed-geometry SCF test;
2. if the SCF is healthy and budget remains, at most approximately five ionic steps.

Measure:

- wall time per electronic iteration;
- electronic iterations per ionic step;
- memory use;
- SCF behaviour;
- force evolution;
- dipole-correction behaviour where applicable;
- electrostatic-potential output;
- restart / execution provenance.

### Pilot B — layered broken-bond surface

After the basal pilot passes, run one bounded pilot on the most plausible layered edge candidate.

A successful basal calculation does not automatically validate broken-bond surfaces, which may show:

- surface states;
- different SCF behaviour;
- reconstruction;
- spin polarisation.

### Pilot interpretation

Classify each pilot as:

- **PASS** — sufficient runtime and convergence evidence to continue;
- **INCONCLUSIVE** — the bounded runtime was insufficient but no scientific failure was demonstrated;
- **HOLD** — a structural, numerical, electrostatic, or execution problem must be resolved first.

A runtime limit is a resource gate, not a scientific convergence criterion.

---

## 5.5 Computational Resource Gate

WP2 must remain compatible with the subsequent adsorption workload.

The following are project-management targets rather than scientific convergence thresholds:

- electronic iterations should ideally remain on a scale of minutes rather than tens of minutes;
- a sustained electronic-iteration cost above approximately 30 min should trigger review before expanding that slab family;
- routine candidate relaxations should preferably complete on approximately day-scale rather than week-scale;
- single relaxations projected to require more than approximately 48 h require explicit justification before propagation.

Total relaxation cost must be estimated from:

\[
T_\mathrm{relax}
\approx
N_\mathrm{ionic}
\times
\overline{N_\mathrm{SCF/ionic}}
\times
\overline{t_\mathrm{electronic}}
\]

rather than from electronic-iteration time alone.

Do not reduce numerical accuracy merely to meet a runtime target.

Instead reconsider:

- slab dimensions;
- unnecessary candidate models;
- convergence strategy;
- parallel layout;
- whether the calculation changes an unresolved scientific conclusion.

---

## 5.6 Stage 02 — Sequential Numerical Convergence

Convergence testing must be **sequential and adaptive**, not a Cartesian-product sweep.

Do not calculate every possible combination of slab thickness, vacuum, k-point density, cutoff, and constraint scheme.

Use one representative model per relevant slab family and converge approximately in the following order:

1. in-plane k-point density and ENCUT;
2. vacuum separation;
3. slab thickness;
4. relaxed-region / central-layer constraint;
5. one final combined validation using the selected settings.

Only extend the convergence dimension that fails.

---

### 5.6.1 Exchange-Correlation and PAW Baseline

Use as the starting protocol:

- PBE;
- D3(BJ), consistent with WP1;
- Zn / In_d / S PAW-PBE datasets used in WP1;
- ENCUT = 500 eV.

Perform a representative:

500 → 600 eV

check on at least:

- one layered basal surface;
- one broken-bond surface.

Do not propagate 600 eV automatically unless the 500 eV result fails the required convergence criterion.

---

### 5.6.2 k-Point Convergence

Use one k point normal to the slab.

Test approximately comparable in-plane reciprocal-space densities:

- ~0.20 Å⁻¹;
- ~0.16 Å⁻¹;
- ~0.12 Å⁻¹.

Record the actual:

\[
N_1 \times N_2 \times 1
\]

mesh and the reciprocal-spacing convention.

Do not blindly reuse the WP1 bulk k mesh.

Stop once successive refinement satisfies the required energetic and electronic convergence criteria.

---

### 5.6.3 Vacuum Convergence

Test actual atom-free periodic-image separations of approximately:

- 20 Å;
- 25 Å;
- 30 Å.

These refer to physical slab-to-image separation, not total cell height.

Check both:

- surface excess energy;
- vacuum-referenced electronic quantities.

For asymmetric neutral slabs, verify that both vacuum regions are sufficiently field-free.

Do not infer energy convergence solely from a flat electrostatic-potential plateau because periodic dispersion interactions may also depend on vacuum separation.

---

### 5.6.4 Slab-Thickness Convergence

#### Layered basal family

Use complete, termination-compatible layered units.

An initial sequence may be approximately:

- 2 complete layers;
- 4 complete layers;
- 6 complete layers.

These values are candidate convergence points, not mandatory calculations.

Only compare slabs that preserve the same termination / stacking family.

For beta / IIb comparisons, match **physical slab thickness and termination character**, not merely the number of crystallographic unit cells.

#### Prismatic / spinel families

Use approximate physical thicknesses such as:

- ~12 Å;
- ~18 Å;
- ~24 Å;

rounded to chemically valid repeat units.

Check whether a bulk-like interior emerges and whether surface reconstruction changes with thickness.

---

### 5.6.5 Relaxed-Region Convergence

Compare:

1. full ionic relaxation;
2. relaxation with a central bulk-like region fixed.

The purpose is to determine whether:

- the surface geometry is artificially coupled through the slab;
- fixing the interior materially changes the surface structure;
- the selected constraint provides a meaningful bulk-like reference.

Record forces on constrained atoms as diagnostic information rather than interpreting the constrained model as fully force-free.

---

### 5.6.6 Convergence Acceptance

Use the paired surface excess quantity as the common energetic convergence metric where possible.

Accept a numerical setting when successive refinement changes:

\[
|\Delta\Gamma_\mathrm{pair}|
\leq
1\ \mathrm{meV\,Å^{-2}}
\]

and relevant vacuum-referenced levels by approximately:

\[
|\Delta E_\mathrm{vac-ref}|
\leq
0.05\ \mathrm{eV}
\]

without changing the qualitative conclusions regarding:

- reconstruction;
- surface-state localisation;
- spin state;
- termination identity.

If later mechanistic conclusions depend on differences smaller than these tolerances, tighten the corresponding convergence test.

---

## 5.7 Stage 03 — Surface Relaxation and Candidate Screening

### Screening relaxation

Use candidate screening to determine whether a termination remains chemically meaningful after relaxation.

A suitable initial screening target is:

- fixed cell;
- PBE+D3(BJ);
- selected converged numerical settings;
- `ISYM = 0`;
- `EDIFF ≤ 1 × 10−6 eV`;
- maximum relaxed-atom force approximately `< 0.02 eV Å−1`.

Inspect:

- reconstruction;
- broken-bond healing;
- surface rumpling;
- coordination changes;
- layer spacing;
- unexpected atom transfer between surfaces;
- convergence behaviour.

If two independently generated candidates relax to effectively the same structure, merge them only after explicit structural comparison.

### Spin test

Broken-bond surfaces must not automatically be assumed nonmagnetic.

For shortlisted broken-bond surfaces:

- compare non-spin-polarised and physically sensible spin-polarised initialisations;
- inspect local magnetic moments as well as total magnetisation;
- retain the lower-energy physically meaningful solution.

Do not interpret zero total magnetisation as proof of a nonmagnetic state.

### Reconstruction test

For shortlisted broken-bond surfaces, test whether a larger lateral cell allows a lower-energy reconstruction.

Use a physically appropriate:

- 2 × 1;
- 1 × 2;
- or, only where justified, larger reconstruction cell.

Introduce small, documented symmetry-breaking surface displacements before relaxation.

Do not treat an unchanged replicated 1 × 1 structure as proof that no reconstruction exists.

### Final relaxation

Only the final reaction-surface shortlist must be tightened to approximately:

- `EDIFF = 1 × 10−7 eV`;
- maximum relaxed-atom force `< 0.01 eV Å−1`.

Do not apply publication-critical force thresholds automatically to every candidate that will later be discarded.

---

## 5.8 Stage 04 — Matched Static Calculations and Surface Energetics

### Matched bulk references

Recalculate the required parent-phase bulk reference energies using:

- the production VASP 6.6.1 executable;
- the same PAW datasets;
- the same functional and D3 treatment;
- the same ENCUT convention;
- sufficiently comparable Brillouin-zone integration.

Retain the accepted WP1 bulk geometries.

This is a matched static reference, not a new bulk geometry-optimisation campaign.

Use the corresponding parent-phase reference:

- beta slab → beta bulk;
- IIb slab → IIb bulk;
- spinel slab → spinel bulk.

Do not subtract all slabs from the lowest-energy polymorph.

### Static surface calculations

For the accepted relaxed slabs, perform self-consistent fixed-geometry static calculations.

Generate the converged quantities required for:

- total energies;
- charge density;
- layer/site-resolved PDOS;
- electrostatic potential;
- band-edge localisation;
- later selected partial charge-density analysis.

Retain a usable WAVECAR for shortlisted surfaces where band-edge partial charge densities will be required.

Do not retain large restart files for every discarded screening model unless needed for recovery.

---

## 5.9 Surface-Energy Convention

For a stoichiometric slab containing \(n\) ZnIn2S4 formula units and one-face area \(A\), define:

\[
\Gamma_\mathrm{pair}
=
\frac{
E_\mathrm{slab}
-
n e_\mathrm{bulk}^{(p)}
}{A}
\]

where \(e_\mathrm{bulk}^{(p)}\) is the matched bulk energy per formula unit of the corresponding parent phase.

\(\Gamma_\mathrm{pair}\) is the combined excess energy of the two slab faces per one-face area.

### Equivalent faces

Only when the two faces are demonstrably equivalent may the single-face surface energy be reported as:

\[
\gamma
=
\frac{
E_\mathrm{slab}
-
n e_\mathrm{bulk}^{(p)}
}{
2A
}
=
\frac{\Gamma_\mathrm{pair}}{2}
\]

### Inequivalent faces

For asymmetric slabs:

- do not assign half of \(\Gamma_\mathrm{pair}\) to either individual face;
- report the paired surface excess quantity;
- treat any two-face average explicitly as an average rather than an independently determined face energy.

Do not directly rank a single-face \(\gamma\) against a two-face \(\Gamma_\mathrm{pair}\) without converting them to a common convention.

### Nonstoichiometric surfaces

Nonstoichiometric terminations require explicit chemical-potential references.

Keep them outside the initial simple stoichiometric surface-energy ranking unless separately justified.

---

## 5.10 Electrostatics and Vacuum Alignment

For asymmetric neutral slabs:

- orient the cell so the surface normal is unambiguous;
- use the appropriate dipole correction where justified;
- calculate electrostatic potentials suitable for vacuum alignment;
- obtain the Hartree + ionic potential using `LVHAR`.

The dipole correction removes periodic-image artefacts; it does not remove intrinsic physical polarity.

### Two-face vacuum levels

Determine the vacuum level separately on both sides of an asymmetric slab.

For each side verify:

- sufficiently low charge density;
- an approximately field-free plateau;
- negligible residual slope across the chosen sampling interval;
- separation from any dipole-correction discontinuity.

Do not average two physically distinct vacuum levels merely to obtain one number.

### Work function

For each face:

\[
\Phi_{\pm}
=
V_{\mathrm{vac},\pm}
-
E_F
\]

Report the convention used for \(E_F\).

For gapped slabs, the numerical Fermi level inside the band gap is convention-dependent and must not automatically be interpreted as the experimental Fermi level of the thin film.

### Vacuum-referenced band edges

Where a sufficiently bulk-like slab interior exists, determine vacuum-referenced bulk-like band edges using an internally consistent electrostatic-potential alignment.

Do not subtract raw eigenvalues from separate bulk and slab calculations without alignment.

Distinguish explicitly:

- surface-localised states;
- finite-slab states;
- bulk-like interior VBM/CBM.

If the slab remains too thin or shows a persistent internal field / electronic reconstruction, report the limitation rather than forcing a thickness-independent band alignment.

Absolute PBE redox alignment remains provisional; hybrid-functional validation is a later and separately budgeted decision.

---

## 5.11 Stage 05 — Surface Electronic Structure

For shortlisted clean surfaces analyse:

- total DOS;
- element/orbital PDOS;
- layer-resolved PDOS;
- site-resolved contributions from under-coordinated surface atoms;
- surface-state localisation;
- bulk-like interior states;
- electrostatic-potential profile;
- face-specific vacuum levels;
- work function;
- local coordination environments.

### Band-edge charge densities

For selected states, generate band-decomposed / partial charge densities to identify whether the VBM, CBM, or relevant in-gap states are:

- surface-localised;
- face-specific;
- distributed through the slab;
- bulk-like in the interior.

Use these together with PDOS and electrostatic alignment.

Do not use a single DOS plot as evidence of carrier transfer or photocatalytic charge-separation dynamics.

### Slab band structures

Full slab band structures are **not a default WP2 requirement**.

Perform them only if a specific unresolved scientific question requires them.

---

## 5.12 Decision Gate 2 — Reaction-Surface Shortlist

At the end of WP2 select up to approximately three reaction-facing surfaces.

### Surface A — beta basal reference

Use the best-supported beta basal slab.

The exact reaction-facing termination must be explicitly identified.

Until experimental structural evidence becomes available, describe this as the **provisional layered basal reference**, not as the experimentally proven dominant facet.

### Surface B — structurally distinct layered edge

Select the most defensible layered edge / prismatic surface that:

- survives relaxation;
- is numerically converged;
- has chemically meaningful exposed sites;
- remains computationally tractable for repeated adsorption calculations.

Prefer structural and electronic distinctiveness over simply choosing the second-lowest surface energy.

### Surface C — optional contrast

Prefer a defensible spinel surface only if the spinel feasibility gate passes and the model provides a meaningful three-dimensional contrast.

Otherwise consider an IIb surface only if it demonstrates a genuine stacking-dependent difference relative to beta.

If neither option adds sufficient scientific information:

> propagate only two reaction surfaces into WP3/WP4 rather than inventing a third.

The reaction-surface shortlist must be based on separate consideration of:

1. energetic accessibility;
2. structural distinctiveness;
3. electronic character;
4. experimental plausibility;
5. computational tractability.

Clean-surface descriptors do not establish catalytic activity.

---

## 5.13 WP3 / WP4 Handoff

For every shortlisted reaction surface record:

- parent phase;
- Miller orientation;
- termination;
- exact reaction-facing side;
- final slab geometry;
- physical surface area;
- slab thickness;
- vacuum;
- k-point mesh;
- fixed/relaxed atoms;
- spin state;
- clean-surface reference energy;
- work-function / vacuum-alignment convention;
- relevant surface-state information;
- recommended lateral supercell expansion for adsorption;
- remaining adsorption-coverage convergence requirements.

A clean-surface unit cell that is sufficient for WP2 convergence is not automatically large enough for realistic adsorption calculations.

WP3 and WP4 must independently verify lateral adsorbate interactions / coverage where necessary.

---

## Deliverables

### Surface registry

Machine-readable record of:

- parent phase;
- orientation;
- termination;
- exposed face;
- composition;
- slab thickness;
- vacuum;
- constraints;
- numerical settings;
- convergence status;
- runtime provenance.

### Convergence summary

Document:

- k-point convergence;
- ENCUT check;
- vacuum convergence;
- slab-thickness convergence;
- relaxed-region convergence;
- reconstruction/spin checks where required.

### Surface energetics

Provide:

- paired surface excess energies;
- single-face surface energies only where mathematically justified;
- matched bulk references;
- structural ranking.

### Surface electronic properties

Provide:

- layer/site-resolved PDOS;
- surface-state identification;
- electrostatic-potential profiles;
- face-specific vacuum levels;
- work functions;
- vacuum-referenced bulk-like band edges where justified;
- selected band-edge charge localisation.

### Final WP3/WP4 shortlist

Provide a compact justification for:

- Surface A;
- Surface B;
- optional Surface C.

State explicitly which candidate surfaces were:

- rejected;
- merged as duplicates;
- held for additional convergence;
- deferred for computational or physical reasons.

---

## Scope Exclusions

WP2 does **not** include:

- O2 adsorption;
- CO2 adsorption;
- H2O2 intermediates;
- *OOH;
- *COOH;
- *OCHO;
- *H;
- charged reaction slabs;
- sulfur vacancies;
- NiOx interfaces;
- CoOx–Pi interfaces;
- NEB;
- exhaustive surface phase diagrams;
- production HSE06 slab calculations;
- automatic propagation of every WP1 polymorph.

These belong to later work packages or separately approved validation tasks.

---

## Final Principle

WP2 should not maximise the number of surfaces calculated.

Its purpose is to establish a small set of reaction-facing models for which:

**structure → numerical convergence → surface energetics → electrostatics → electronic localisation**

is sufficiently well understood that subsequent reaction calculations are scientifically interpretable.

Every additional termination, slab thickness, reconstruction, or convergence point should resolve a specific remaining uncertainty.

Do not calculate a model merely because a matrix still contains an empty cell.