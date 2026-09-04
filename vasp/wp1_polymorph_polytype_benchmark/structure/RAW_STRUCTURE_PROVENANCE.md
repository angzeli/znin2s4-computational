# ZnIn2S4 Structure Provenance

This document records the provenance, intended role, and required preprocessing of the ZnIn2S4 crystal structures used in WP1:

**Polymorph / Polytype Benchmark**

The purpose is to preserve a clear distinction between:

- experimental source structures;
- reconstructed structures;
- disordered parent structures;
- legacy/control structures;
- final DFT-relaxed production structures.

The raw CIF files must not be assumed to be equivalent to the final structures used for production calculations.

---

# 1. Target Production Structure Set

The intended final WP1 production comparison is:

1. cubic spinel ZnIn2S4
2. ordered α1 / IIIa-derived ZnIn2S4
3. β-ZnIn2S4
4. revised IIa′-ZnIn2S4
5. IIb-ZnIn2S4

The legacy experimental IIa structure is retained only as the starting structure from which IIa′ is expected to form during unconstrained DFT relaxation.

---

# 2. Raw Structure Directory

Raw structures are stored in:

`vasp/wp1_polymorph_polytype_benchmark/structure/raw/`

The files currently relevant to WP1 are:

- spinel_cod_7221365.cif
- IIIa_cod_1525846.cif
- IIa_legacy_cod_1539811.cif
- IIb_cod_1527302.cif
- beta_reconstructed_from_IIb.cif

---

# 3. Cubic Spinel ZnIn2S4

## File

spinel_cod_7221365.cif

## Source

Crystallography Open Database:

COD 7221365

Experimental cubic spinel ZnIn2S4 structure.

The structure originates from the experimental crystallographic literature on spinel ZnIn2S4.

## Expected identity

- phase: cubic spinel
- space group: Fd-3m
- space-group number: 227
- experimental lattice parameter: approximately 10.59 Å

## Role

Direct experimental starting structure for the cubic polymorph.

## Required preprocessing

None beyond ordinary CIF validation and conversion to the VASP structure format.

## WP1 workflow

spinel experimental CIF
→ full geometry optimisation
→ relaxed spinel production structure

## Status

READY AS A DIRECT DFT SEED

---

# 4. α / IIIa ZnIn2S4 Parent Structure

## File

IIIa_cod_1525846.cif

## Source

Crystallography Open Database:

COD 1525846

Experimental IIIa-ZnIn2S4 structure associated with the crystallographic work of Lopez-Rivera et al.

Approximate experimental lattice parameters:

- a ≈ 3.87 Å
- c ≈ 37.07 Å

## Important structural feature

This experimental structure contains occupational disorder on tetrahedral cation sites.

These sites are represented by approximately:

- 50% Zn
- 50% In

Therefore, the experimental CIF is not directly suitable for a deterministic VASP calculation.

VASP requires an explicit atomistic configuration and cannot optimise a fractional Zn/In occupancy into an ordered structure automatically.

## Relation to Lee et al. 2019

Lee et al. treated this occupational disorder by generating symmetry-distinct ordered configurations using the Site Occupancy Disorder (SOD) approach.

They obtained three representative configurations:

- α1
- α2
- α3

The α1 configuration was the lowest-energy structure.

Reported optB86b result for α1:

- space group: R3m
- a ≈ 3.89 Å
- c ≈ 36.86 Å

α2 and α3 were substantially higher in energy.

## Role

The current CIF is therefore:

**the experimental disordered α / IIIa parent structure**

and not the final α production structure.

## Required preprocessing

Resolve the Zn/In occupational disorder before VASP.

Preferred workflow:

IIIa disordered experimental CIF
→ generate ordered Zn/In configurations
→ identify/reconstruct α1
→ full geometry optimisation
→ relaxed α1 production structure

It is not necessary to carry α2 and α3 into the later surface/reaction calculations once α1 has been established as the representative ordered structure.

## Status

VALID RAW EXPERIMENTAL SOURCE

NOT YET A DIRECT PRODUCTION DFT SEED

---

# 5. β-ZnIn2S4

## Parent structure

IIb_cod_1527302.cif

## Reconstructed file

beta_reconstructed_from_IIb.cif

## Reconstruction provenance

The β structure was reconstructed from experimental IIb-ZnIn2S4.

Lee et al. describe the layered structures in terms of an:

S–Zn–S–In–S–In–S

septuple layer.

β contains one septuple layer per unit cell.

IIb contains two such septuple layers per unit cell, with the two layers related by a 60° rotation about the c axis.

The reconstruction therefore consisted of:

experimental IIb
→ identify one complete septuple layer
→ use that layer as the periodic repeat
→ construct β

The two independent septuple layers present in IIb produce crystallographically equivalent β structures up to the expected 60° rotational relationship.

## Reconstructed structure

- composition: ZnIn2S4
- atom count: 7
- space group: P3m1
- space-group number: 156
- a = b = 3.850000 Å
- c = 12.340000 Å
- α = β = 90°
- γ = 120°

A very small symmetry standardisation was applied:

- coordinates near 1/3 and 2/3 were snapped to exact thirds
- maximum displacement ≈ 0.00022 Å

This displacement is crystallographically negligible.

## Reference DFT target

Lee et al. reported for β using optB86b:

- a ≈ 3.89 Å
- c ≈ 12.28 Å

These values can be used as a geometry-relaxation validation target.

## Role

Experimental-geometry-derived reconstructed β starting structure.

It is not an independently experimentally refined β CIF.

## WP1 workflow

experimental IIb
→ reconstructed β seed
→ full geometry optimisation
→ relaxed β production structure

## Status

READY AS A DIRECT DFT SEED

---

# 6. Legacy IIa-ZnIn2S4

## File

IIa_legacy_cod_1539811.cif

## Source

Crystallography Open Database:

COD 1539811

Experimental historical IIa structure.

Approximate experimental structure:

- space group: P-3m1
- space-group number: 164
- a ≈ 3.85 Å
- c ≈ 24.68 Å

## Important scientific status

This structure is intentionally retained even though Lee et al. 2019 argued that it is not the stable ground-state IIa structure.

Their calculations showed that when the legacy IIa structure is allowed to undergo a fully unconstrained geometry optimisation, the atoms spontaneously reorganise into a different structure referred to as:

IIa′

The restructuring occurs through Zn–S bond rearrangement.

Lee et al. reported that this behaviour occurred with multiple exchange-correlation treatments.

The legacy IIa structure was also reported to be mechanically/dynamically unstable compared with IIa′.

## Critical workflow decision

Do NOT manually construct IIa′ before the initial WP1 geometry optimisation.

Instead, use the experimental legacy IIa structure directly as the initial DFT structure.

The geometry optimisation itself should test whether the Lee et al. restructuring can be reproduced.

## WP1 workflow

legacy experimental IIa
→ fully unconstrained geometry optimisation
→ spontaneous structural rearrangement
→ IIa′
→ relaxed IIa′ production structure

## Reference DFT target

Lee et al. reported for IIa′ using optB86b:

- space group: P-3m1
- a ≈ 3.90 Å
- c ≈ 24.76 Å

For PBE+D3BJ they reported approximately:

- a ≈ 3.89 Å
- c ≈ 24.56 Å

These values provide useful validation targets.

## Important calculation requirement

The geometry optimisation must allow both:

- ionic relaxation
- lattice relaxation

No selective-dynamics constraint should prevent the relevant Zn–S rearrangement.

The calculation should not artificially preserve the legacy geometry.

## Role

Historical experimental starting structure and reconstruction seed for IIa′.

## Status

INTENTIONAL LEGACY / TRANSFORMATION SEED

NOT THE FINAL PRODUCTION STRUCTURE

---

# 7. IIb-ZnIn2S4

## File

IIb_cod_1527302.cif

## Source

Crystallography Open Database:

COD 1527302

Experimental IIb-ZnIn2S4 structure.

The structure originates from the crystallographic study of the two-pack IIb polytype.

## Expected identity

- space group: P63mc
- space-group number: 186
- a ≈ 3.85 Å
- c ≈ 24.68 Å
- two ZnIn2S4 formula units per cell
- 14 atoms after complete symmetry expansion
- no occupational disorder

## Structural interpretation

IIb contains two ZnIn2S4 septuple layers per crystallographic repeat.

The two layers are related by an approximately 60° rotation about the c axis.

## Role

Direct experimental starting structure for IIb.

It is also the crystallographic parent used to reconstruct the β seed.

## WP1 workflow

IIb experimental CIF
→ full geometry optimisation
→ relaxed IIb production structure

## Reference DFT target

Lee et al. reported with optB86b:

- a ≈ 3.89 Å
- c ≈ 24.58 Å

## Status

READY AS A DIRECT DFT SEED

---

# 8. Summary of Structure Status

| Structure | Current raw source | Directly VASP-ready? | Required action before final production structure |
|---|---|---:|---|
| spinel | spinel_cod_7221365.cif | Yes | Full geometry optimisation |
| α1 | IIIa_cod_1525846.cif | No | Resolve Zn/In disorder, generate α1, then optimise |
| β | beta_reconstructed_from_IIb.cif | Yes | Full geometry optimisation |
| IIa′ | IIa_legacy_cod_1539811.cif | Yes as transformation seed | Full unconstrained optimisation should generate IIa′ |
| IIb | IIb_cod_1527302.cif | Yes | Full geometry optimisation |

---

# 9. Final WP1 Structure-Generation Logic

Only one structure requires explicit preprocessing before ordinary VASP optimisation:

**α1**

because the experimental IIIa parent contains fractional Zn/In occupancies.

The other structures can enter the geometry-optimisation stage directly.

The expected workflow is therefore:

spinel
→ relax
→ spinel

IIIa disordered parent
→ resolve cation ordering
→ α1
→ relax
→ α1

β reconstructed from IIb
→ relax
→ β

legacy IIa
→ fully unconstrained relax
→ IIa′

IIb experimental
→ relax
→ IIb

---

# 10. Provenance Principle

Never describe all five input structures simply as "experimental CIFs".

Their provenance is different:

- spinel: experimental
- α1: ordered model derived from a disordered experimental parent
- β: reconstructed from the experimental IIb structure
- IIa′: DFT-relaxed structure generated from the experimental legacy IIa structure
- IIb: experimental

These distinctions should be preserved in all later calculation records, figures, methods sections, and publication data.
