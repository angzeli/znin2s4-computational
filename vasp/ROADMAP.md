# Publication-Oriented VASP Roadmap for ZnIn2S4 Thin-Film Photocatalysis

## Project Context

This computational project is designed to support the final-year project:

**Precursor-derived ZnIn2S4 thin-film architectures for photocatalytic CO2 reduction**

The experimental project compares four core thin-film architectures:

1. pristine ZnIn2S4
2. NiOx/ZnIn2S4
3. ZnIn2S4/Co-Pi
4. NiOx/ZnIn2S4/Co-Pi

The architecture is spatially divided into three functional regions:

- **ZnIn2S4 layer:** visible-light absorption, charge generation, and reduction chemistry
- **NiOx bottom interface:** buried hole-selective contact intended to extract photogenerated holes away from ZnIn2S4
- **Co-Pi top surface:** reaction-facing hole-management cocatalyst intended to accept or mediate holes while retaining reducing electrons on ZnIn2S4

The experimental programme also proposes photocatalytic O2 reduction to H2O2 as an initial probe of photogenerated-electron utilisation before progressing to photocatalytic CO2 reduction.

The computational programme should therefore not be treated as an isolated DFT exercise. Its purpose is to establish a mechanistic framework connecting:

**crystal structure → exposed surface → electron-transfer chemistry → interface charge management → H2O2 screening → CO2 activation**

The detailed computational workflow below extends beyond the literal scope of the literature review and is proposed as a publication-oriented strategy for making the one-month VASP campaign scientifically coherent.

---

# 1. Central Scientific Question

The main question is:

> How do ZnIn2S4 crystal structure, exposed surface, NiOx bottom-interface engineering, and Co-Pi top-surface engineering collectively control charge localisation and the ability of photogenerated electrons to drive H2O2 formation and CO2 activation?

A second, more specific question is:

> Is photocatalytic H2O2 production a useful mechanistic proxy for the CO2-reduction capability of ZnIn2S4 thin films, or does it mainly report efficient electron delivery to the surface?

A third question concerns the two hole-management components:

> Are NiOx and Co-Pi complementary hole sinks, redundant components, or competing destinations for photogenerated holes?

These questions define the computational scope.

---

# 2. Working Publication Hypothesis

The initial working hypothesis is:

1. Different ZnIn2S4 polymorphs and polytypes exhibit distinct band-edge structures, carrier localisation, and surface chemistry.

2. These structural differences alter:
   - O2 adsorption and activation
   - H2O2-forming 2e− oxygen reduction
   - CO2 activation
   - early CO2-reduction intermediate stabilisation
   - competition with hydrogen evolution

3. NiOx modifies the buried ZnIn2S4 interface and may provide a lower-energy localisation site for photogenerated holes.

4. Co-Pi provides a surface-localised hole sink but may simultaneously alter:
   - ZnIn2S4 surface accessibility
   - surface work function
   - O2 and CO2 adsorption
   - H2O2 stability

5. The full NiOx/ZnIn2S4/Co-Pi architecture is beneficial only if directional hole management improves charge separation without excessively blocking or electronically perturbing ZnIn2S4 reduction sites.

The calculations should be designed so that these hypotheses can also be falsified.

---

# 3. Overall Strategy

The project is divided into three priority tiers.

## P0 — Core publication package

These calculations define the minimum scientifically complete project.

- WP1 — ZnIn2S4 polymorph and polytype benchmark
- WP2 — ZnIn2S4 surface and facet screening
- WP3 — H2O2 / 2e− ORR reaction chemistry
- WP4 — Minimal CO2RR + HER descriptor map
- WP5 — NiOx/ZnIn2S4 buried-interface calculations
- WP6 — Co-Pi/ZnIn2S4 surface-interface calculations

## P1 — High-value mechanistic controls

Perform these after the P0 framework is working.

- WP7 — H2O2 adsorption and decomposition
- WP8 — realistic-film controls:
  - sulfur vacancy
  - partially sulfurised NiOx/ZnIn2S4 interface

## P2 — Optional final integration

Only perform if sufficient computational time remains.

- WP9 — representative NiOx/ZnIn2S4/Co-Pi combined stack

The full trilayer model must not delay the pairwise interface calculations.

---

# 4. WP1 — ZnIn2S4 Crystal-Structure Benchmark

## Objective

Establish which ZnIn2S4 structures are energetically plausible and determine how crystal structure affects the electronic properties relevant to photocatalysis.

## Initial structure set

At minimum evaluate:

- cubic ZnIn2S4
- rhombohedral ZnIn2S4
- representative layered hexagonal ZnIn2S4 structures

Where reliable structures are available, extend the layered family to multiple reported polytypes or stacking arrangements, for example:

- α-type
- β-type
- IIa-type
- IIb-type

The exact crystallographic labels should be verified against the structural source used to generate each model.

## Calculations

For every bulk structure:

### Geometry optimisation

Optimise:

- lattice vectors
- atomic positions
- cell volume

Use one consistent exchange-correlation and dispersion treatment across all structures.

### Relative stability

Calculate:

- total energy
- energy per ZnIn2S4 formula unit
- relative energy against the lowest-energy phase

### Structural descriptors

Record:

- lattice parameters
- cell volume
- density
- characteristic layer spacing
- relevant Zn–S and In–S bond lengths

### Simulated XRD

Generate simulated powder XRD patterns.

Purpose:

- allow direct comparison with experimental thin-film XRD
- help identify which calculated structure is experimentally relevant

### Electronic structure

Calculate:

- band structure
- total DOS
- projected DOS
- VBM orbital character
- CBM orbital character
- band-decomposed charge density

For the most relevant structures:

- HSE06 band gap
- HSE06 band-edge character
- SOC sensitivity
- electron effective mass
- hole effective mass

## Decision Gate 1

At the end of WP1:

### Retain all structures for the bulk comparison figure.

### Shortlist only 2–3 structures for expensive surface chemistry.

Selection should consider:

1. calculated stability
2. experimental XRD compatibility
3. structural distinctiveness
4. electronic-structure differences
5. computational tractability

Do not propagate every candidate structure into the full adsorption matrix.

## Deliverables

- table of optimised structural parameters
- relative-energy diagram
- simulated XRD comparison
- band-gap comparison
- VBM/CBM character comparison
- shortlist for WP2

---

# 5. WP2 — Surface and Facet Screening

## Objective

Identify the surfaces that are most relevant to reduction chemistry and determine whether photocatalytic behaviour is facet-dependent.

Bulk band structures alone are insufficient because O2 reduction and CO2 reduction occur at surfaces.

## Initial slab set

Construct approximately 6–8 candidate surfaces.

### Layered structures

Include:

- basal surface
- representative edge/prismatic surface

If chemically meaningful:

- S-rich termination
- alternative stoichiometric termination

### Cubic or rhombohedral structures

Test 1–2 low-index surfaces.

Prefer:

- stoichiometric slabs
- symmetric slabs
- non-polar or compensated surfaces

Avoid uncritically using strongly polar surfaces.

## Convergence

For each representative slab family test:

- slab thickness
- number of relaxed layers
- vacuum thickness
- k-point density

Target vacuum:

approximately 18–20 Å as an initial value, followed by convergence testing.

Use dipole correction for asymmetric slabs.

## Calculated properties

For each clean surface:

- surface energy
- relaxed geometry
- work function
- surface PDOS
- VBM localisation
- CBM localisation
- electrostatic potential
- surface dipole
- surface reconstruction
- surface Zn/In/S coordination environments

Where absolute band alignment is required:

- align electronic levels to vacuum
- carefully distinguish surface-state energies from bulk-like band edges

## Decision Gate 2

Shortlist approximately three reaction surfaces:

### Surface A
experimentally most plausible dominant surface

### Surface B
structurally distinct surface from the same or another layered phase

### Surface C
cubic or rhombohedral contrast

These surfaces become the common reaction set for WP3 and WP4.

## Deliverables

- surface-energy ranking
- work-function map
- facet-resolved electronic structures
- reaction-surface shortlist

---

# 6. WP3 — H2O2 Formation and 2e− Oxygen Reduction

## Objective

Determine how ZnIn2S4 phase and surface control O2 activation and the thermodynamics of the H2O2-forming oxygen-reduction pathway.

This work directly supports the planned experimental use of H2O2 production as a preliminary photocatalytic screening reaction.

---

## 6.1 O2 Adsorption

For each shortlisted surface test multiple O2 geometries:

- end-on adsorption
- side-on adsorption
- Zn-associated adsorption
- In-associated adsorption
- S-associated adsorption
- bridge configurations where appropriate

For each stable configuration record:

- adsorption energy
- O–O bond length
- Bader charge
- spin density
- magnetic moment
- local structural distortion

The key question is whether adsorbed O2 remains approximately molecular O2 or gains significant superoxide/peroxide-like character.

---

## 6.2 Electron-Added O2 Activation

Perform selected calculations with one additional electron in the slab.

Compare:

### neutral system

surface + O2

against:

### electron-added system

surface + O2 + e−

Analyse:

- electron localisation
- Bader charge transferred to O2
- O–O bond elongation
- spin-density localisation
- energy preference for surface-bound reduced O2

These calculations are intended primarily as a relative localisation probe.

Do not directly interpret uncorrected charged-slab total energies as accurate absolute redox potentials.

Perform supercell-size or vacuum sensitivity checks for the key systems.

---

## 6.3 2e− ORR Intermediates

Calculate the key species:

- O2
- *OOH
- *O
- *OH
- H2O2*

The primary H2O2 route is:

O2 → *OOH → H2O2

The competing O–O-cleavage branch should also be considered:

*OOH → *O-containing products

The purpose is to separate:

### activity

ability to form *OOH

from:

### selectivity

ability to retain the O–O bond and release H2O2 rather than continue toward water or strongly bound oxygen species

---

## 6.4 Free-Energy Treatment

Use a consistent free-energy framework.

Include where feasible:

- electronic energy
- zero-point-energy correction
- entropy correction
- solvation correction

For proton-coupled electron-transfer bookkeeping, a CHE-like treatment may be used as a consistent thermodynamic framework.

However:

> The resulting free-energy profile should not be interpreted as a direct photocatalytic reaction rate or electrocatalytic overpotential.

Where possible, evaluate the sensitivity of the free-energy profile to a physically relevant photogenerated-electron chemical potential derived from:

- experimental band-edge information
- calculated vacuum-aligned CBM
- experimentally estimated flat-band potential

---

## 6.5 H2O2 Desorption

Explicitly calculate:

- adsorbed H2O2 stability
- H2O2 desorption energy

Strong formation but excessively strong H2O2 adsorption may reduce net detectable H2O2 production.

## Deliverables

For each shortlisted surface:

- O2 adsorption map
- electron-transfer analysis
- *OOH free energy
- H2O2 formation profile
- O–O retention/selectivity analysis
- H2O2 desorption energy

Final figure:

**facet-dependent H2O2 activity/selectivity map**

---

# 7. WP4 — Minimal CO2 Reduction and HER Map

## Objective

Test whether surfaces favourable for H2O2 formation are also favourable for the earliest mechanistically important steps in CO2 reduction.

The goal is not to calculate the complete pathway to CH4 or CH3OH.

The goal is to capture:

1. CO2 adsorption
2. CO2 activation
3. first protonation branch
4. early CO/formate selectivity
5. HER competition

---

# 7.1 CO2 Adsorption

For each shortlisted surface test:

- linear physisorbed CO2
- C-bound CO2
- O-bound CO2
- bidentate CO2
- other chemically reasonable configurations

Record:

- adsorption energy
- O–C–O angle
- C–O bond lengths
- Bader charge
- spin density
- surface distortion

---

# 7.2 Electron-Induced CO2 Activation

For representative configurations calculate:

surface + CO2 + e−

Analyse formation of bent CO2-like species.

Key descriptors:

- CO2 bending angle
- charge transferred to CO2
- C–O elongation
- spin localisation
- relative stabilisation of the activated state

This provides a direct counterpart to electron-induced O2 activation.

---

# 7.3 Early CO2RR Branches

Calculate:

- *COOH
- *OCHO
- *CO
- HCOOH-associated state where appropriate

The minimum reaction network is:

CO2 → *COOH → *CO → CO

and

CO2 → *OCHO → HCOOH

The purpose is to compare the preference for:

- carboxyl-type chemistry
- formate-type chemistry

without extending into deep multi-electron reduction.

---

# 7.4 HER Competitor

Calculate:

- *H adsorption

Obtain a consistent ΔG*H descriptor.

This allows each surface to be placed on a simple competition map:

CO2 activation vs H adsorption

---

# 7.5 H2O2–CO2 Correlation Analysis

This is one of the central outputs of the project.

Compare across surfaces:

- ΔG*OOH vs ΔG*COOH
- ΔG*OOH vs ΔG*OCHO
- O2 electron-trapping tendency vs CO2 electron-trapping tendency
- work function vs O2 activation
- work function vs CO2 activation
- CBM position vs activated-CO2 charge
- H2O2 selectivity descriptor vs CO2 branch preference
- ΔG*H vs CO2 activation

## Interpretation A — Strong correlation

If H2O2-active surfaces are also systematically better at CO2 activation:

> H2O2 formation may be a useful preliminary descriptor for reduction-side electron availability in this ZnIn2S4 system.

## Interpretation B — Weak or absent correlation

If the rankings differ:

> H2O2 formation mainly demonstrates efficient electron transfer to O2 and should not be interpreted as direct evidence of strong CO2 activation chemistry.

This negative result would still provide a mechanistically important conclusion.

## Deliverables

- CO2 adsorption map
- CO2 activation map
- *COOH/*OCHO comparison
- HER competition map
- H2O2–CO2 descriptor correlation figure

---

# 8. WP5 — NiOx/ZnIn2S4 Buried Interface

## Objective

Determine whether the NiOx underlayer is electronically capable of acting as a hole-selective contact beneath ZnIn2S4.

The central question is not simply:

"Does NiO form a heterojunction?"

It is:

> Does the realistic NiOx/ZnIn2S4 interface thermodynamically favour hole localisation away from ZnIn2S4 without creating severe recombination-active interface states?

---

# 8.1 Bulk NiO Calibration

NiO must be treated as a correlated magnetic oxide.

Use:

- spin-polarised calculations
- correct antiferromagnetic ordering
- DFT+U or another justified correlated-electron treatment

Calibrate against known qualitative properties:

- band gap
- magnetic moment
- electronic structure

Do not arbitrarily choose U only to force a desired band alignment.

Where practical, perform an HSE06 spot check.

---

# 8.2 Stoichiometric NiO

Construct the clean stoichiometric NiO surface/interface model.

Calculate:

- geometry
- work function
- DOS
- magnetic moments
- band alignment reference

---

# 8.3 Ni-Deficient NiOx

Because experimentally relevant p-type NiOx is commonly associated with nickel vacancies and oxidised Ni environments, construct at least one Ni-deficient model.

Minimum models:

- stoichiometric NiO
- VNi-containing NiOx

Where feasible compare:

- vacancy near interface
- vacancy further from interface

Analyse:

- hole localisation
- Ni oxidation-state changes
- local magnetic moments
- defect states

---

# 8.4 Interface Construction

Perform coincidence-lattice matching between ZnIn2S4 and NiO/NiOx.

Target:

- minimal in-plane strain
- preferably below approximately 3–4% if computationally achievable

Construct at least:

### Interface 1
stoichiometric NiO/ZnIn2S4 — registry A

### Interface 2
stoichiometric NiO/ZnIn2S4 — registry B

### Interface 3
VNi-NiOx/ZnIn2S4 — registry A

### Interface 4
VNi-NiOx/ZnIn2S4 — registry B

If one registry is clearly unstable, reduce the matrix after relaxation.

---

# 8.5 Interface Analysis

For each production interface calculate:

- interface adhesion energy
- interface formation energy where meaningful
- charge-density difference
- Bader charge
- plane-averaged electrostatic potential
- interface dipole
- layer-resolved DOS
- projected DOS
- local magnetic moments
- band offsets

When determining band offsets:

- account for interface potential shifts
- avoid relying only on isolated bulk band energies

---

# 8.6 Hole-Localisation Calculations

This is a key calculation.

Remove one electron from the interface system and test different initial localisation conditions.

Attempt to initialise the hole on:

1. ZnIn2S4
2. NiOx near the interface
3. NiOx away from the interface

Relax each state.

Compare:

- total energies
- spin densities
- Bader charge
- local magnetic moments
- structural distortions

The goal is to determine the thermodynamic preference for hole localisation.

Possible conclusions:

### Outcome A

Hole strongly prefers NiOx.

Interpretation:

NiOx can plausibly act as a hole-selective sink.

### Outcome B

Hole remains on ZnIn2S4.

Interpretation:

The intended hole-extraction mechanism is not strongly supported by the interface energetics.

### Outcome C

Hole localises in a deep interface state.

Interpretation:

NiOx may introduce trapping rather than productive extraction.

---

# 8.7 Important Physical Limitation

The model must explicitly recognise that NiOx is buried.

Even if hole transfer from ZnIn2S4 to NiOx is favourable, static DFT cannot prove that holes are subsequently consumed efficiently.

Therefore distinguish:

- favourable hole localisation
- actual long-range hole transport
- productive oxidation chemistry

The latter two require experimental validation or more advanced dynamics.

## Deliverables

- NiO/NiOx electronic calibration
- interface structures
- band alignment
- charge-density difference
- interface potential profile
- hole-localisation map

Final mechanistic question:

> Is NiOx a hole-selective contact, a passive interface, or a trap-forming layer?

---

# 9. WP6 — Co-Pi/ZnIn2S4 Surface Interface

## Objective

Model how reaction-facing cobalt-phosphate-derived domains influence ZnIn2S4 surface electronic structure and hole localisation.

---

# 9.1 Modelling Principle

Co-Pi is an amorphous and hydrated cobalt oxo/hydroxo-phosphate material.

Therefore:

> Do not claim that a single periodic crystalline structure represents real Co-Pi.

Instead use an ensemble of chemically reasonable:

**CoOx–Pi structural proxies**

The final manuscript should use this terminology explicitly.

---

# 9.2 Structural Proxy Set

Construct approximately 2–3 representative motifs.

Possible motifs include:

- edge-sharing CoO6 dimer
- small CoOx cluster
- phosphate-coordinated CoOx cluster
- hydroxylated/hydrated CoOx–Pi motif

Where computationally feasible vary:

- protonation
- Co oxidation state
- hydration
- local Co coordination

---

# 9.3 Anchoring-Site Comparison

Place the same representative cluster on:

- ZnIn2S4 basal surface
- ZnIn2S4 edge surface

Where useful:

- cubic/rhombohedral contrast surface

Test multiple anchoring geometries.

Calculate:

- binding energy
- Co–surface bonding geometry
- structural distortion
- charge transfer
- work-function change

Main question:

> Does CoOx–Pi preferentially occupy relatively inactive surface regions or block intrinsically active ZnIn2S4 reduction sites?

---

# 9.4 Hole Localisation

Remove one electron from the CoOx–Pi/ZnIn2S4 system.

Test initial hole localisation on:

- ZnIn2S4
- Co centre
- Co–O framework
- nearby surface S

Analyse:

- spin density
- Bader charge
- Co magnetic moments
- projected DOS
- local structural relaxation

Determine whether hole localisation becomes more favourable on the Co-containing domain.

---

# 9.5 Coverage Effect

Calculate at least two qualitative coverage regimes:

### Low coverage

isolated CoOx–Pi domain

### Higher coverage

reduced lateral spacing and/or multiple domains

Compare:

- work function
- exposed ZnIn2S4 area
- surface electronic states
- neighbouring O2 adsorption
- neighbouring CO2 adsorption

This addresses the experimental possibility that excessive Co-Pi loading may:

- shield light
- block reduction sites
- inhibit reactant adsorption
- alter charge-transfer pathways

---

# 9.6 Adjacent Reduction Chemistry

For the most stable CoOx–Pi/ZnIn2S4 structure, test whether O2 and CO2 can still bind to nearby exposed ZnIn2S4.

Compare against bare ZnIn2S4:

- O2 adsorption
- CO2 adsorption
- charge transfer
- key intermediate stabilisation

This distinguishes:

### beneficial electronic modification

from:

### physical site blocking

## Deliverables

- CoOx–Pi structural-proxy ensemble
- anchoring-energy comparison
- basal vs edge preference
- work-function modification
- hole-localisation analysis
- coverage-effect analysis

---

# 10. WP7 — H2O2 Stability and Decomposition

Priority: P1

## Objective

Determine whether Co-Pi may improve charge management while simultaneously decreasing the net accumulation of H2O2.

Experimentally measured H2O2 concentration reflects:

H2O2 formation − H2O2 decomposition

Therefore H2O2 production cannot automatically be interpreted only in terms of formation kinetics.

---

# 10.1 H2O2 Adsorption

Calculate H2O2 adsorption on:

- bare ZnIn2S4
- CoOx–Pi/ZnIn2S4

Compare:

- adsorption energy
- O–O bond length
- charge transfer
- orientation
- hydrogen bonding

---

# 10.2 Decomposition Thermodynamics

Test plausible initial decomposition products.

Examples may include:

- O–O-cleaved configurations
- OH-containing surface species
- dehydrogenated H2O2-derived intermediates

If decomposition is strongly exergonic on Co-containing sites:

proceed to NEB.

---

# 10.3 Optional CI-NEB

Perform one high-value NEB comparison:

H2O2 decomposition on bare ZnIn2S4

versus

H2O2 decomposition at/near CoOx–Pi

Only perform if preliminary thermodynamics justify it.

Use approximately 5–7 images initially.

## Experimental implication

Recommend a dark H2O2 decay control across Co-Pi loadings.

This would separate:

- H2O2 formation
- H2O2 consumption

## Deliverables

- H2O2 adsorption comparison
- decomposition thermodynamics
- optional decomposition barrier

---

# 11. WP8 — Real-Film Chemistry Controls

Priority: P1

These calculations test whether conclusions derived from ideal structures survive more realistic defects and processing chemistry.

---

# 11.1 Surface Sulfur Vacancy

Construct one representative surface sulfur vacancy.

Calculate:

- vacancy formation energy
- defect electronic states
- work-function change
- O2 adsorption
- CO2 adsorption
- *OOH stabilisation
- *COOH and/or *OCHO stabilisation
- *H adsorption

Purpose:

> Determine whether phase/facet trends remain meaningful when a realistic sulfur defect is present.

Sulfur vacancy should remain a control rather than becoming the central project theme.

---

# 11.2 Partial Sulfurisation of NiOx

The experimental fabrication sequence places Zn/In xanthate precursors above a pre-formed NiOx underlayer before thermal conversion to ZnIn2S4.

Therefore the buried oxide interface may experience a sulfur-rich chemical environment during annealing.

Test at least one model representing:

- O → S substitution near the interface

or

- local Ni–S bond formation

Compare against the ideal NiOx/ZnIn2S4 interface.

Calculate:

- sulfurisation reaction energy
- interface stability
- adhesion
- DOS
- interface states
- charge transfer
- hole localisation

Possible outcome:

> The functional buried interface may be better described as a graded Ni–O–S/ZnIn2S4 region rather than an atomically ideal NiOx/ZnIn2S4 junction.

This would be a potentially high-value result if strongly supported.

---

# 12. WP9 — Representative Full NiOx/ZnIn2S4/Co-Pi Stack

Priority: P2

## Objective

Construct one combined model only after the pairwise interfaces are understood.

This calculation is primarily qualitative.

Do not use a thin computational stack to make strong claims about long-range carrier transport in the experimental film.

## Model

Use:

- best-supported NiOx/ZnIn2S4 interface
- best-supported ZnIn2S4/CoOx–Pi surface motif
- minimum ZnIn2S4 thickness that still provides an approximately bulk-like internal region

Calculate:

- layer-resolved DOS
- electrostatic potential
- charge-density difference
- hole localisation
- electron localisation where meaningful

Main question:

> When both hole-management components are present, where does an introduced hole preferentially localise?

Compare approximate hole-state energies or relaxed hole configurations for:

- ZnIn2S4
- NiOx
- CoOx–Pi

Interpret the result as a qualitative competition between hole sinks.

---

# 13. Potential Mechanistic Outcomes

## Scenario 1 — Co-Pi-dominated hole management

Hole localisation:

CoOx–Pi < NiOx < ZnIn2S4

Interpretation:

- surface Co-Pi is the strongest hole sink
- NiOx may provide secondary extraction
- oxidation chemistry remains surface-accessible

---

## Scenario 2 — NiOx-dominated hole management

Hole localisation:

NiOx < CoOx–Pi < ZnIn2S4

Interpretation:

- buried NiOx strongly attracts holes
- the key experimental question becomes whether these holes can be productively consumed or transported

---

## Scenario 3 — Comparable NiOx and Co-Pi hole energies

Interpretation:

Carrier kinetics, film thickness, morphology, and interface quality may determine the actual charge distribution.

This could explain non-additive experimental behaviour.

---

## Scenario 4 — Deep interface trapping

If either modifier introduces strongly localised mid-gap states:

Interpretation:

apparent hole extraction may instead correspond to charge trapping and potentially increased recombination.

---

# 14. One-Month Execution Schedule

# Week 1 — Crystal Structures and Surface Foundation

## Primary goals

Complete:

- bulk structure preparation
- geometry optimisation
- relative phase energies
- simulated XRD
- PBE electronic structures
- initial HSE06 jobs
- candidate surface construction
- slab convergence

## Tasks

### Days 1–2

- organise structure sources
- standardise POSCARs
- standardise pseudopotentials
- perform ENCUT convergence
- perform bulk k-point convergence

### Days 2–4

Run bulk relaxations for:

- cubic
- rhombohedral
- layered/polytype structures

### Days 3–5

Generate:

- energies
- XRD
- DOS
- band structures

Launch HSE06 for shortlisted structures.

### Days 4–7

Construct and relax approximately 6–8 surface models.

Perform:

- vacuum convergence
- slab-thickness convergence
- surface energy
- work-function calculations

## End-of-week decision gate

Select:

- 2–3 important bulk structures
- 3 reaction surfaces

Do not expand the reaction matrix beyond this shortlist.

---

# Week 2 — H2O2 and CO2 Reaction Chemistry

## Primary goals

Complete the reaction-side descriptor framework before starting very large interfaces.

## H2O2 branch

For three shortlisted surfaces:

- O2 adsorption
- electron-added O2
- *OOH
- *O
- *OH
- H2O2 adsorption/desorption

## CO2 branch

For the same surfaces:

- CO2 adsorption
- electron-added CO2
- *COOH
- *OCHO
- *CO
- *H

## Analysis

Generate preliminary:

- H2O2 free-energy diagrams
- CO2 early-pathway diagrams
- HER descriptors
- H2O2–CO2 correlation plots

## End-of-week decision gate

Identify:

- most H2O2-active surface
- most H2O2-selective surface
- most CO2-activating surface
- strongest HER surface
- whether H2O2 and CO2 rankings correlate

This result determines which surface should dominate later interface calculations.

---

# Week 3 — Interface Engineering

Run NiOx and Co-Pi branches in parallel.

## NiOx branch

Complete:

- AFM NiO calibration
- stoichiometric NiO
- VNi-NiOx
- interface matching
- 2–4 interface relaxations
- interface DOS
- charge-density difference
- hole-localisation calculations

## Co-Pi branch

Complete:

- 2–3 CoOx–Pi structural proxies
- basal anchoring
- edge anchoring
- low-coverage model
- higher-coverage model
- hole-localisation calculations
- work-function analysis

## End-of-week decision gate

Determine:

### NiOx

- does it attract holes?
- does VNi improve or worsen the interface?
- are deep traps generated?

### Co-Pi

- where does it preferentially bind?
- does it accept holes?
- does it block reduction sites?

---

# Week 4 — Validation and Publication-Level Completion

Do not indiscriminately start new calculations.

Use Week 4 to strengthen the most important conclusions.

Priority order:

1. convergence and reproducibility
2. HSE06/SOC validation
3. important alternative configurations
4. realistic-film controls
5. NEB
6. full stack

## Recommended Week 4 additions

Select only the most justified calculations from:

- H2O2 decomposition NEB
- surface sulfur vacancy
- partial NiOx sulfurisation
- full NiOx/ZnIn2S4/CoOx–Pi stack

## Final tasks

- rerun failed calculations
- tighten critical geometries
- verify magnetic ground states
- verify adsorption minima
- calculate Bader charges consistently
- generate publication figures
- archive structures
- document all settings
- export machine-readable summary tables

---

# 15. Calculation Standards

The entire project must use a consistent computational protocol.

---

# 15.1 Pseudopotentials

Use one fixed PAW potential set throughout the project.

Record explicitly:

- VASP version
- PAW dataset date/version
- exact POTCAR choices

Do not mix pseudopotential families between comparisons.

---

# 15.2 Exchange-Correlation Functional

Recommended baseline:

- PBE
- one consistent dispersion correction for layered ZnIn2S4

Examples include a D3-type treatment or another validated vdW correction.

The exact method should be selected before production calculations and then kept fixed.

---

# 15.3 Plane-Wave Cutoff

Perform explicit ENCUT convergence.

Production ENCUT should be:

- safely above the maximum ENMAX requirement
- sufficiently converged for energy differences and forces

A value around 500 eV may be a sensible starting point, but the final choice must follow convergence testing.

---

# 15.4 Electronic Convergence

Recommended target:

EDIFF ≤ 1 × 10−6 eV

For difficult magnetic/interface systems:

use robust SCF settings and avoid interpreting unconverged charge states.

---

# 15.5 Ionic Relaxation

Standard production target:

maximum force < 0.02 eV Å−1

For publication-critical structures:

maximum force < 0.01 eV Å−1

Examples:

- lowest-energy interface
- key *OOH state
- key *COOH state
- hole-localised structures
- NEB endpoints

---

# 15.6 k-Point Convergence

Converge bulk k-point density explicitly.

For slabs:

maintain comparable reciprocal-space density in the periodic directions.

Large interface supercells may permit Γ-centred low-density meshes, but convergence must be checked.

---

# 15.7 Spin Polarisation

Mandatory for:

- O2
- reduced O2 states
- NiO
- NiOx
- CoOx–Pi
- charged slabs
- radical intermediates where relevant

Test multiple initial magnetic moments for Co-containing models.

Do not assume the first converged magnetic solution is the ground state.

---

# 15.8 DFT+U

Use only where justified.

Likely required for:

- Ni 3d
- potentially Co 3d

The selected Ueff values must be documented and tested against relevant electronic/magnetic properties.

Perform sensitivity analysis for the most important conclusions where practical.

---

# 15.9 HSE06

HSE06 should be used as a validation layer rather than for every production geometry.

Recommended HSE06 set:

- important ZnIn2S4 phases
- dominant reaction surface
- key NiOx electronic structure
- selected NiOx/ZnIn2S4 interface
- selected CoOx–Pi/ZnIn2S4 model

Approximate target:

8–12 HSE06 validation calculations

depending on computational cost.

---

# 15.10 Spin–Orbit Coupling

Perform SOC spot checks for:

- ZnIn2S4 bulk band edges
- key HSE06 structures if affordable

Purpose:

determine whether SOC materially changes the conclusions.

Do not automatically apply SOC to every adsorption relaxation.

---

# 15.11 Surface Calculations

For slabs:

- converge thickness
- converge vacuum
- use dipole corrections when asymmetric
- minimise adsorbate–adsorbate interactions
- maintain chemically meaningful stoichiometry

For polar/non-stoichiometric surfaces:

do not report naive surface energies without considering the relevant termination and chemical-potential dependence.

---

# 15.12 Charged-Slab Calculations

For electron-added and hole-added slabs:

- use the same supercell for comparisons
- perform supercell-size sensitivity
- check vacuum sensitivity
- analyse localisation rather than only absolute total energy

Charged calculations should primarily support:

- localisation tendencies
- relative stabilisation trends
- qualitative redox-state analysis

---

# 15.13 Solvation

If implicit solvation is available:

perform solvation sensitivity for key adsorbates such as:

- *OOH
- H2O2
- *COOH
- *OCHO

Where practical:

add limited explicit water molecules around strongly hydrogen-bonded intermediates.

The same solvation convention must be used across the H2O2 and CO2 datasets.

---

# 16. Minimum Publishable Calculation Package

The approximate minimum target is:

## Bulk

- approximately 5–6 ZnIn2S4 structural models

## Surfaces

- approximately 6–8 clean surfaces
- approximately 3 shortlisted reaction surfaces

## Reaction chemistry

approximately 24–30 adsorption/intermediate structures covering:

- O2
- *OOH
- *O
- *OH
- H2O2
- CO2
- activated CO2
- *COOH
- *OCHO
- *CO
- *H

## NiOx interfaces

approximately 4 production interface models

## Co-Pi

approximately 4–6 structural/coverage models

## Controls

- 1 sulfur-vacancy model
- 1 partially sulfurised NiOx interface

## NEB

- 1–2 high-value pathways maximum

## Hybrid/SOC validation

approximately 8–12 calculations

Expected total:

approximately 50–70 production calculation chains

Many chains consist of:

geometry optimisation → static calculation → DOS/Bader/potential analysis

rather than completely independent calculations.

---

# 17. Scope Exclusions

The following should NOT be attempted during the one-month campaign unless the entire core roadmap is already complete.

---

## 17.1 Full CO2 → CH4 Mechanism

Do not calculate every PCET step to methane.

Reason:

- too many intermediates
- strong solvent dependence
- protonation-state ambiguity
- many kinetic barriers
- likely to produce a superficially complete but physically weak pathway

---

## 17.2 Full CO2 → CH3OH Mechanism

Same reason.

The early *COOH/*OCHO branching is sufficient for the current mechanistic objective.

---

## 17.3 Exhaustive C2 Chemistry

Do not calculate C–C coupling unless experimental results later justify it.

---

## 17.4 Giant Trilayer Models Before Pairwise Interfaces

Do not begin with:

NiOx/ZnIn2S4/Co-Pi

as one very large model.

First understand:

NiOx/ZnIn2S4

and

ZnIn2S4/Co-Pi

independently.

---

## 17.5 Single DOS as Evidence of Charge Transfer

Do not claim charge separation or photocatalytic transfer pathways from DOS alone.

Use multiple indicators:

- electrostatic potential
- charge-density difference
- Bader charge
- spin density
- hole-localisation calculations
- band alignment

---

## 17.6 S-Scheme Claims from Ground-State DFT Alone

Do not claim an S-scheme mechanism solely from:

- band alignment
- Bader charge
- charge-density difference

Ground-state DFT can support:

- equilibrium charge redistribution
- interface dipoles
- carrier-localisation preference

It cannot by itself prove the actual photoexcited carrier-transfer mechanism.

---

# 18. Decision Rules for Computational Scope

Use explicit stop rules.

## If a bulk phase is > approximately 0.2–0.3 eV/f.u. above the lowest-energy phase

Do not perform the full surface reaction matrix unless experiment suggests its presence.

---

## If two polytypes produce essentially identical reaction-facing surfaces

Retain both in the bulk structural analysis but calculate detailed chemistry for only one representative.

---

## If an interface has extremely high strain

Reject it rather than forcing an unrealistic commensurate cell.

---

## If an adsorbate repeatedly desorbs during optimisation

Treat this as physically meaningful weak adsorption rather than forcing an artificial bound state.

---

## If CoOx–Pi structural proxies give qualitatively different conclusions

Do not select one convenient structure.

Report the model dependence.

---

## If H2O2 and CO2 descriptors clearly decouple

Do not try to force a correlation.

That decoupling becomes a central mechanistic conclusion.

---

# 19. Suggested Data Organisation

Recommended project structure:

VASP-FYP/
│
├── 00_protocol/
│   ├── POTCAR_manifest
│   ├── convergence
│   ├── INCAR_templates
│   └── methodology_notes
│
├── 01_bulk_ZIS/
│   ├── cubic
│   ├── rhombohedral
│   ├── alpha
│   ├── beta
│   ├── IIa
│   └── IIb
│
├── 02_surfaces/
│   ├── phase_surface_A
│   ├── phase_surface_B
│   └── ...
│
├── 03_H2O2_ORR/
│   ├── surface_A
│   ├── surface_B
│   └── surface_C
│
├── 04_CO2RR/
│   ├── surface_A
│   ├── surface_B
│   └── surface_C
│
├── 05_NiOx_interface/
│   ├── bulk_NiO
│   ├── VNi_NiOx
│   ├── interface_01
│   └── ...
│
├── 06_CoPi_interface/
│   ├── proxy_01
│   ├── proxy_02
│   └── coverage
│
├── 07_controls/
│   ├── sulfur_vacancy
│   ├── sulfurised_NiOx
│   └── H2O2_decomposition
│
├── 08_full_stack/
│
├── analysis/
│   ├── bader
│   ├── band_alignment
│   ├── free_energy
│   ├── figures
│   └── tables
│
└── metadata/
    ├── structure_sources
    ├── calculation_registry
    └── failed_jobs

---

# 20. Calculation Registry

Every calculation should be indexed in a machine-readable registry.

Recommended fields:

- calculation_id
- work_package
- material
- phase
- surface
- defect
- adsorbate
- charge
- spin_state
- functional
- Ueff
- ENCUT
- k_mesh
- supercell
- slab_thickness
- vacuum
- convergence_status
- final_energy
- maximum_force
- notes

This prevents a one-month high-throughput campaign from becoming impossible to audit later.

---

# 21. Publication Figure Plan

A publication-quality computational story could be organised around approximately six main figures.

---

## Figure 1 — ZnIn2S4 Structural Landscape

Include:

- crystal structures
- relative energies
- simulated XRD
- band gaps
- VBM/CBM alignment

Message:

> ZnIn2S4 phase/polytype identity materially changes the electronic landscape relevant to photocatalysis.

---

## Figure 2 — Facet-Dependent Reduction Chemistry

Include:

- surface structures
- work functions
- O2 adsorption
- *OOH energetics
- H2O2 selectivity descriptors

Message:

> H2O2 formation is controlled by reaction-facing ZnIn2S4 structure rather than bulk composition alone.

---

## Figure 3 — H2O2 as a CO2RR Proxy

Include:

- CO2 activation geometry
- *COOH
- *OCHO
- *H
- descriptor-correlation plots

Message:

either:

> H2O2 screening captures a transferable reduction-side electronic descriptor.

or:

> H2O2 activity and CO2 activation are mechanistically decoupled.

Both outcomes are scientifically useful.

---

## Figure 4 — NiOx/ZnIn2S4 Buried Interface

Include:

- interface structure
- electrostatic-potential profile
- layer-resolved DOS
- charge-density difference
- hole spin density

Message:

> NiOx changes hole localisation at the buried interface.

---

## Figure 5 — CoOx–Pi/ZnIn2S4 Surface Engineering

Include:

- structural proxies
- basal vs edge anchoring
- work-function changes
- hole localisation
- coverage effects

Message:

> Co-Pi acts as a reaction-facing hole-management domain whose benefit depends on spatial distribution and coverage.

---

## Figure 6 — Unified Architecture Mechanism

Compare:

- pristine ZnIn2S4
- NiOx/ZnIn2S4
- ZnIn2S4/Co-Pi
- NiOx/ZnIn2S4/Co-Pi

Summarise:

- electron-active surface
- preferred hole localisation
- reduction-site accessibility
- H2O2 behaviour
- CO2 activation tendency

Message:

> Bottom-interface and top-surface engineering can be complementary, redundant, or competitive depending on the relative energetics of charge localisation and surface-site accessibility.

---

# 22. Experimental–Computational Cross-Validation

The calculations should generate predictions that can be compared directly with the FYP experiments.

---

## XRD

DFT:

- predicted phase
- simulated XRD

Experiment:

- phase assignment
- crystallinity
- preferred orientation

---

## UV–Vis

DFT:

- phase-dependent band gaps
- electronic transitions qualitatively

Experiment:

- optical band gap
- absorption onset

---

## XPS

DFT:

- charge redistribution
- local coordination changes
- Ni/Co oxidation-state trends qualitatively

Experiment:

- Zn/In/S/Ni/Co chemical environments
- Ni3+/Ni2+ trends
- Co oxidation-state changes
- possible interfacial sulfurisation

---

## PL

DFT:

- indirect evidence through carrier localisation tendencies

Experiment:

- recombination changes

Do not claim quantitative PL intensity prediction from static DFT.

---

## Photocurrent / EIS

DFT:

- interface band alignment
- hole-localisation tendency

Experiment:

- charge extraction
- interfacial transfer resistance

---

## KPFM / Surface Photovoltage

DFT:

- work-function changes
- surface potential
- interface dipoles

Experiment:

- contact-potential differences
- light-induced surface-potential changes

---

## H2O2 Photocatalysis

DFT:

- O2 activation
- *OOH energetics
- H2O2 desorption
- decomposition tendency

Experiment:

- apparent H2O2 accumulation rate

---

## CO2 Photoreduction

DFT:

- CO2 activation
- *COOH/*OCHO preference
- HER competition

Experiment:

- CO
- HCOOH/formate
- CH4
- H2
- other detected products

The DFT should primarily rationalise early-pathway tendencies rather than claim to predict the complete experimental product distribution.

---

# 23. Core Success Criteria

The computational campaign should be considered successful if it answers the following questions robustly.

## Structural

- Which ZnIn2S4 phase/polytype is most relevant?
- Does phase identity materially change the band structure?
- Which facets are plausible reaction surfaces?

## H2O2

- Which surface best activates O2?
- Which surface best stabilises *OOH?
- Which surface is most selective toward retaining the O–O bond?
- Is H2O2 weakly enough bound to desorb?

## CO2

- Which surface best activates CO2?
- Is *COOH or *OCHO preferred?
- How strong is HER competition?
- Does H2O2 activity correlate with CO2 activation?

## NiOx

- Does NiOx attract holes away from ZnIn2S4?
- Does Ni deficiency strengthen hole selectivity?
- Does the interface generate problematic trap states?

## Co-Pi

- Where does CoOx–Pi preferentially bind?
- Does it localise holes?
- Does it block ZnIn2S4 reduction sites?
- Does high loading become detrimental?
- Does it catalyse H2O2 decomposition?

## Combined Architecture

- Are NiOx and Co-Pi complementary?
- Are they competing hole sinks?
- Is one component redundant?
- Can the four experimental architectures be rationalised within one consistent electronic-structure framework?

---

# 24. Priority if Computational Resources Become Limited

If queue time or computational cost becomes severe, preserve the following in this order:

1. ZnIn2S4 bulk phase benchmark
2. surface shortlist
3. H2O2 reaction chemistry
4. minimal CO2RR + HER map
5. NiOx/ZnIn2S4 interface
6. CoOx–Pi/ZnIn2S4 interface
7. HSE06 validation
8. sulfur-vacancy control
9. partial NiOx sulfurisation
10. H2O2 decomposition NEB
11. full trilayer

The full trilayer should be the first major model removed if time becomes limited.

---

# 25. Minimum Scientific Story

If only the P0 package is completed, the resulting study should still support a coherent story:

> ZnIn2S4 photocatalysis is controlled by a hierarchy of structural and interfacial effects. Crystal phase and surface termination define the intrinsic reduction chemistry, including O2 activation, H2O2 formation, and early CO2-reduction pathways. NiOx modifies the buried interface and the thermodynamic preference for hole localisation, whereas Co-Pi introduces a reaction-facing hole-management pathway whose effectiveness depends on its binding geometry and coverage. Comparing these effects establishes whether H2O2 formation can serve as a transferable screening probe for CO2 photoreduction and whether bottom-interface and top-surface engineering act cooperatively in precursor-derived ZnIn2S4 thin films.

---

# 26. Stronger Publication Story if P1 Is Completed

If H2O2 decomposition and realistic-interface controls are also completed:

> The apparent performance of spatially engineered ZnIn2S4 films emerges from competition between reduction-site chemistry, hole localisation, cocatalyst coverage, H2O2 stability, and processing-induced interface chemistry. This framework moves beyond simple band-alignment diagrams and provides an experimentally testable description of how buried NiOx and surface Co-Pi independently and jointly alter photocatalytic function.

---

# 27. Final Principle

The project should not aim to prove that every modification is beneficial.

The strongest possible outcome is a mechanistically resolved answer to:

> What does each architectural component actually do?

A useful computational study may conclude that:

- one ZnIn2S4 phase is electronically superior but chemically less active;
- H2O2 production does not predict CO2 activation;
- NiOx extracts holes but introduces traps;
- Co-Pi improves hole management but decomposes H2O2;
- excessive Co-Pi blocks reduction sites;
- NiOx and Co-Pi compete rather than synergise;
- partial sulfurisation changes the real buried interface.

Any of these outcomes would be more scientifically valuable than forcing the calculations to reproduce the initial working hypothesis.

The target is therefore not a collection of favourable DFT numbers.

The target is a coherent:

**structure → surface → reaction → interface → architecture**

mechanistic model that can be directly tested against the final-year-project experiments.