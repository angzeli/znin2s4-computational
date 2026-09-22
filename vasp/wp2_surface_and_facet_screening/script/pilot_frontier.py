"""Electron-count-based frontier diagnostics for the audited nonmagnetic pilots.

No file I/O, gap certification, spin folding or occupation rescaling occurs here.
The reader must establish a completed final solution and its geometry identity.
"""

import numpy as np


FRONTIER_RESOLUTION_EV = 1e-4  # Reporting/triage choice, not physical accuracy.
NELECT_INTEGER_TOLERANCE = 1e-6
WEIGHT_SUM_TOLERANCE = 1e-6
ARRAY_PRECISION_TOLERANCE = 5.1e-5  # XML prints energies/occupations to 4 decimals.
KPOINT_PRECISION_TOLERANCE = 1e-7
OCCUPATION_TOLERANCE = 1e-3
OCCUPATION_CONVENTION = 'VASP_6.6.1_ISPIN1_FULL_1'


def electronic_diagnostic(eigenvalues, *, nelect=None, ispin=None, nbands=None,
                          kpoints=None, weights=None, noncollinear=None, soc=None,
                          occupation_convention=None):
    """Use bands N and N+1 at every k point, with N=NELECT/2 (one-based).

    Arrays are {one spin channel: (nk, nb, [energy_eV, occupation])}, with
    increasing band order, fractional reciprocal coordinates (nk, 3), and
    normalized nonnegative weights (nk,). Occupations are unchanged, full=1.
    Electron-count tolerance is nb*1e-6 electrons: a conservative bound for
    six-decimal EIGENVAL occupations after normalizing printed k weights.
    Unsupported or inconsistent evidence returns unavailable values and a reason.
    """
    result = dict(frontier_band_separation_eV=None,
                  minimum_direct_frontier_separation_eV=None,
                  frontier_valence_band_index=None, frontier_conduction_band_index=None,
                  frontier_interpretation='UNSUPPORTED_OR_UNRESOLVED',
                  validation_status='UNSUPPORTED_OR_UNRESOLVED',
                  frontier_resolution_eV=FRONTIER_RESOLUTION_EV,
                  occupation_convention=occupation_convention,
                  partial_occupation_present=None, partial_occupation_count=None)

    def unresolved(reason):
        return dict(result, reason=reason)

    if ispin != 1 or noncollinear is not False or soc is not False:
        return unresolved('Only explicitly non-SOC, collinear ISPIN=1 is supported')
    if occupation_convention != OCCUPATION_CONVENTION:
        return unresolved('Unsupported or unverified occupation convention; no rescaling guessed')
    if (not isinstance(nelect, (int, float, np.number)) or isinstance(nelect, (bool, np.bool_))
            or not np.isfinite(nelect) or nelect <= 0):
        return unresolved('Completed-output NELECT must be finite and positive')
    integer_nelect = int(round(nelect))
    if abs(nelect-integer_nelect) > NELECT_INTEGER_TOLERANCE or integer_nelect % 2:
        return unresolved('NELECT is not an even integer within 1e-6 electrons')
    n = integer_nelect//2
    if not isinstance(eigenvalues, dict) or len(eigenvalues) != 1:
        return unresolved('Exactly one nonmagnetic eigenvalue channel is required')
    try:
        values = np.asarray(next(iter(eigenvalues.values())), dtype=float)
        coordinates = np.asarray(kpoints, dtype=float)
        k_weights = np.asarray(weights, dtype=float)
    except (TypeError, ValueError):
        return unresolved('Eigenvalues, coordinates or weights are not numeric arrays')
    if (values.ndim != 3 or values.shape[-1] != 2 or not values.shape[0]
            or not np.isfinite(values).all()):
        return unresolved('Finite (nk, nb, 2) eigenvalue/occupation arrays are required')
    nk, nb, _ = values.shape
    if (not isinstance(nbands, (int, np.integer)) or isinstance(nbands, (bool, np.bool_))
            or nbands != nb or nb < n+1):
        return unresolved('Completed-output NBANDS must match the array and include N+1')
    if (coordinates.shape != (nk, 3) or k_weights.shape != (nk,)
            or not np.isfinite(coordinates).all() or not np.isfinite(k_weights).all()
            or np.any(k_weights < 0) or abs(k_weights.sum()-1) > WEIGHT_SUM_TOLERANCE):
        return unresolved('Finite matching coordinates and normalized k weights are required')
    if nk > 1:
        differences = coordinates[:, None, :] - coordinates[None, :, :]
        differences -= np.rint(differences)
        duplicate = np.max(np.abs(differences), axis=2) < KPOINT_PRECISION_TOLERANCE
        np.fill_diagonal(duplicate, False)
        if duplicate.any():
            return unresolved('Duplicate periodic k coordinates; ordering/weights unresolved')
    energies, occupations = values[:, :, 0], values[:, :, 1]
    if np.any(np.diff(energies, axis=1) < -1e-6):
        return unresolved('Bands are not in increasing VASP energy order')
    if np.any(occupations < -1e-6) or np.any(occupations > 1+1e-6):
        return unresolved('Occupations violate the audited full=1 convention')
    normalized_weights = k_weights/k_weights.sum()
    electrons = float(2*np.sum(normalized_weights[:, None]*occupations))
    count_tolerance = nb*1e-6
    result.update(nelect_check=electrons, electron_count_tolerance=count_tolerance,
                  original_weight_sum=float(k_weights.sum()),
                  occupation_min=float(occupations.min()), occupation_max=float(occupations.max()))
    if abs(electrons-nelect) > count_tolerance:
        return unresolved('Weighted occupations do not reproduce NELECT with one factor of two')
    valence, conduction = energies[:, n-1], energies[:, n]
    separation = float(conduction.min()-valence.max())
    partial = (occupations > OCCUPATION_TOLERANCE) & (occupations < 1-OCCUPATION_TOLERANCE)
    interpretation = ('SAMPLED_BAND_OVERLAP' if separation < -FRONTIER_RESOLUTION_EV
                      else 'FINITE_SAMPLED_FRONTIER_SEPARATION' if separation > FRONTIER_RESOLUTION_EV
                      else 'NEAR_TOUCHING_UNRESOLVED')
    return dict(result, frontier_band_separation_eV=separation,
                minimum_direct_frontier_separation_eV=float(np.min(conduction-valence)),
                frontier_valence_band_index=n, frontier_conduction_band_index=n+1,
                frontier_interpretation=interpretation, validation_status='VALIDATED',
                partial_occupation_present=bool(partial.any()), partial_occupation_count=int(partial.sum()),
                reason='Sampled frontier manifold only; no production gap or experimental metallicity claim')


def crosscheck_eigenval_xml(values, xml_values, kpoints, xml_kpoints, weights, xml_weights, nelect):
    """Check the full unmodified arrays in their original k/band ordering.

    Coordinates tolerate printed precision and periodic integer offsets, but
    arrays are never sorted or paired by nearest energy. XML electron-count
    tolerance is nb*1e-4 electrons, reflecting four-decimal occupations.
    """
    arrays = [np.asarray(x, dtype=float) for x in (values, xml_values, kpoints, xml_kpoints, weights, xml_weights)]
    a, b, k, xml_k, w, xml_w = arrays
    if (any(not np.isfinite(x).all() for x in arrays) or a.shape != b.shape
            or a.ndim != 3 or a.shape[-1] != 2 or k.shape != xml_k.shape
            or k.shape != (a.shape[0], 3) or w.shape != xml_w.shape or w.shape != (a.shape[0],)):
        raise ValueError('EIGENVAL/XML shapes or finite-data evidence disagree')
    difference = k-xml_k
    difference -= np.rint(difference)
    if np.max(np.abs(difference)) > KPOINT_PRECISION_TOLERANCE:
        raise ValueError('EIGENVAL/XML k-coordinate ordering disagrees')
    if (np.any(w < 0) or np.any(xml_w < 0) or abs(w.sum()-1) > WEIGHT_SUM_TOLERANCE
            or abs(xml_w.sum()-1) > WEIGHT_SUM_TOLERANCE
            or np.max(np.abs(w/w.sum()-xml_w/xml_w.sum())) > KPOINT_PRECISION_TOLERANCE):
        raise ValueError('EIGENVAL/XML k weights disagree')
    maximum = float(np.max(np.abs(a-b)))
    if maximum > ARRAY_PRECISION_TOLERANCE:
        raise ValueError('EIGENVAL/XML energies or occupations disagree beyond printed precision')
    xml_electrons = float(2*np.sum((xml_w/xml_w.sum())[:, None]*b[:, :, 1]))
    tolerance = a.shape[1]*1e-4
    if abs(xml_electrons-nelect) > tolerance:
        raise ValueError('XML full=1 occupations do not reproduce NELECT')
    return dict(eigenval_xml_max_difference=maximum, xml_nelect_check=xml_electrons,
                xml_electron_count_tolerance=tolerance, eigenval_xml_status='VALIDATED')
