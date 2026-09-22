"""Pure slab electrostatics; Stage 01 window thresholds and semantics are preserved.

Potential is LVHAR in eV; density is electrons/A^3. Candidate plateaus require
human review and never release a workflow gate. No file or execution side effects.
"""

import numpy as np


def fit_window(z, potential, density=None):
    """Fit a candidate window without assuming it is a field-free plateau."""
    if len(z)<3 or not np.isfinite(potential).all() or np.ptp(z)<=0:
        return dict(status='AMBIGUOUS', reason='Insufficient finite samples')
    slope, intercept=np.polyfit(z,potential,1)
    residual=float(np.max(np.abs(potential-(slope*z+intercept))))
    width=float(np.ptp(z))
    flat=width>=2.0 and abs(slope)<=0.005 and residual<=0.02 and np.ptp(potential)<=0.03
    maximum_density=float(np.max(np.abs(density))) if density is not None else None
    if maximum_density is not None and maximum_density>1e-5:
        status='CHARGE_IN_VACUUM_WARNING'
    else:
        status='CANDIDATE_PLATEAU' if flat else 'NONFLAT_OR_TOO_SHORT'
    return dict(status=status,z_start_A=float(z[0]),z_end_A=float(z[-1]),width_A=width,
                points=len(z),mean_potential_eV=float(np.mean(potential)),slope_eV_per_A=float(slope),
                detrended_max_residual_eV=residual,potential_range_eV=float(np.ptp(potential)),
                mean_density_e_per_A3=float(np.mean(density)) if density is not None else None,
                maximum_abs_density_e_per_A3=maximum_density,
                density_status='AVAILABLE' if density is not None else 'UNAVAILABLE')


def vacuum_windows(z, potential, slab_min, slab_max, height, density=None):
    """Inspect both geometric vacuum regions, splitting visible potential jumps.

    Exclude 2 A next to each exposed face and 1 A at the periodic boundary.
    A discontinuity is a neighboring potential jump over 0.1 eV and over the
    median jump plus ten median absolute deviations. Retain every segment;
    competing candidate plateaus produce an ambiguous result.
    """
    result={}
    for side,low,high in [('lower',1.0,slab_min-2.0),('upper',slab_max+2.0,height-1.0)]:
        indices=np.flatnonzero((z>=low)&(z<=high))
        if len(indices)<3:
            result[side]=dict(status='AMBIGUOUS',reason='Insufficient geometric vacuum')
            continue
        full=fit_window(z[indices],potential[indices],None if density is None else density[indices])
        differences=np.abs(np.diff(potential[indices]))
        median=np.median(differences)
        threshold=max(0.1,float(median+10*np.median(np.abs(differences-median))))
        splits=np.flatnonzero(differences>threshold)+1
        segments=np.split(indices,splits)
        windows=[fit_window(z[part],potential[part],None if density is None else density[part])
                 for part in segments if len(part)>=3]
        candidates=[w for w in windows if w['status']=='CANDIDATE_PLATEAU']
        if len(candidates)==1:
            status='CANDIDATE_PLATEAU'
            selected=candidates[0]
        elif len(candidates)>1:
            status='AMBIGUOUS_MULTIPLE_WINDOWS'
            selected=None
        else:
            status='NO_CONFIRMED_PLATEAU'
            selected=None
        result[side]=dict(status=status,geometric_window=full,segments=windows,selected_window=selected,
                          discontinuity_count=len(splits),jump_threshold_eV=threshold)
    lower=result['lower'].get('selected_window')
    upper=result['upper'].get('selected_window')
    result['upper_minus_lower_potential_eV']=(upper['mean_potential_eV']-lower['mean_potential_eV']
                                              if upper and lower else None)
    result['interpretation']='Candidate windows only; human review required. No forced plateau or vacuum alignment.'
    return result


def same_geometry(first, second):
    """Require equivalent site order, cell and periodic coordinates for output pairing."""
    if first.species!=second.species or not np.allclose(first.lattice.matrix,second.lattice.matrix,atol=1e-6,rtol=0):
        return False
    difference=first.frac_coords-second.frac_coords
    return bool(np.allclose(difference-np.rint(difference),0,atol=1e-6,rtol=0))


def vacuum_window_sensitivity(z, potential, slab_min, slab_max, height, density=None):
    """Retain both faces at fixed 2/4/5/6 A setbacks; never select a best window.

    The boundary exclusion stays 1 A and the original automatic assessment is
    unchanged. A jump uses the same threshold as vacuum_windows; an interval
    containing one is ambiguous and is never linearly fitted across the jump.
    Missing density cannot certify a depleted-vacuum window.
    """
    z, potential = np.asarray(z), np.asarray(potential)
    density = None if density is None else np.asarray(density)
    valid = (z.ndim == 1 and len(z) >= 3 and potential.shape == z.shape
             and np.isfinite(z).all() and np.isfinite(potential).all()
             and np.all(np.diff(z) > 0) and np.isfinite([slab_min, slab_max, height]).all()
             and 0 <= slab_min < slab_max <= height)
    density_valid = density is not None and density.shape == z.shape and np.isfinite(density).all()
    automatic = vacuum_windows(z, potential, slab_min, slab_max, height,
                               density if density_valid else None) if valid else {}
    records = []
    for setback in (2, 4, 5, 6):
        for side, low, high in [('lower', 1., slab_min-setback),
                                ('upper', slab_max+setback, height-1.)]:
            record = dict(side=side, setback_A=setback, boundary_exclusion_A=1.,
                          requested_z_start_A=low, requested_z_end_A=high,
                          automatic_status=automatic.get(side, {}).get('status', 'UNAVAILABLE'),
                          status='INVALID_WINDOW', reason=None,
                          width_criterion=None, slope_criterion=None, residual_criterion=None,
                          range_criterion=None, density_criterion=None, discontinuity_count=None)
            records.append(record)
            if not valid or low >= high or low < 0 or high > height:
                record['reason'] = 'Invalid profile, geometry or requested physical interval; window not moved'
                continue
            indices = np.flatnonzero((z >= low) & (z <= high))
            record.update(points=len(indices), z_start_A=float(z[indices[0]]) if len(indices) else None,
                          z_end_A=float(z[indices[-1]]) if len(indices) else None,
                          width_A=float(np.ptp(z[indices])) if len(indices) else None)
            if len(indices) < 3 or np.ptp(z[indices]) < 2.:
                record.update(status='INSUFFICIENT_WINDOW', width_criterion=False,
                              reason='Fewer than three samples or less than 2 A width; window not moved')
                continue
            differences = np.abs(np.diff(potential[indices]))
            median = np.median(differences)
            threshold = max(.1, float(median+10*np.median(np.abs(differences-median))))
            jumps = int(np.count_nonzero(differences > threshold))
            record.update(discontinuity_count=jumps, jump_threshold_eV=threshold, width_criterion=True)
            if jumps:
                record.update(status='AMBIGUOUS_DISCONTINUITY',
                              reason='Potential jump inside the requested interval; no cross-jump fit or selection')
                continue
            record.update(fit_window(z[indices], potential[indices],
                                     density[indices] if density_valid else None))
            record.update(width_criterion=record['width_A'] >= 2.,
                          slope_criterion=abs(record['slope_eV_per_A']) <= .005,
                          residual_criterion=record['detrended_max_residual_eV'] <= .02,
                          range_criterion=record['potential_range_eV'] <= .03,
                          density_criterion=record['maximum_abs_density_e_per_A3'] <= 1e-5 if density_valid else None)
            if not density_valid:
                record.update(status='INSUFFICIENT_DENSITY', reason='Matching finite charge density unavailable')
            elif record['status'] != 'CANDIDATE_PLATEAU':
                failed = [name for name in ('width', 'slope', 'residual', 'range', 'density')
                          if record[name+'_criterion'] is False]
                record['reason'] = 'Failed screening criteria: '+', '.join(failed)
    return records



def planar_profile(locpot, *, charge=None, nelect=None):
    """Analyze already parsed LVHAR/CHGCAR objects without any filesystem access.

    Caller must establish LVHAR output type and common self-consistent source.
    VASP's raw CHGCAR grid mean is NELECT; dividing its planar mean by cell
    volume gives electrons/A^3. Mismatched density is unavailable, not zero.
    This stricter reusable adapter does not change the Stage 01 file reader.
    """
    structure = locpot.structure
    matrix = structure.lattice.matrix
    if (not np.allclose(matrix[:2, 2], 0, atol=1e-10)
            or not np.allclose(matrix[2, :2], 0, atol=1e-10)):
        raise ValueError('A slab-normal third lattice vector is required')
    potential = np.asarray(locpot.get_average_along_axis(2), dtype=float)
    if potential.ndim != 1 or len(potential) < 3 or not np.isfinite(potential).all():
        raise ValueError('Finite plane-averaged potential samples are required')
    height = float(matrix[2, 2])
    z = np.arange(len(potential)) * height / len(potential)
    density = None; density_status = 'UNAVAILABLE'; integral = None
    if charge is not None:
        try:
            if not same_geometry(charge.structure, structure):
                raise ValueError('CHGCAR geometry differs from LOCPOT')
            if charge.dim != locpot.dim:
                raise ValueError('Density/potential grids differ; no interpolation')
            if nelect is None or not np.isfinite(nelect) or nelect <= 0:
                raise ValueError('NELECT validation is required')
            integral = float(np.mean(charge.data['total']))
            if not np.isfinite(integral) or abs(integral-nelect) > max(1e-3, 1e-5*nelect):
                raise ValueError('CHGCAR integral does not reproduce NELECT')
            values = np.asarray(charge.get_average_along_axis(2)) / structure.volume
            if values.shape != z.shape or not np.isfinite(values).all():
                raise ValueError('Finite matching density/potential grids required; no interpolation')
            density = values; density_status = 'VALIDATED_NELECT_INTEGRAL'
        except ValueError as exc:
            density_status = 'UNAVAILABLE: ' + str(exc)
    atom_z = structure.cart_coords[:, 2]
    return dict(z_A=z, lvhar_potential_eV=potential, density_e_per_A3=density,
                charge_density_status=density_status, charge_integral_electrons=integral,
                vacuum=vacuum_windows(z, potential, float(atom_z.min()), float(atom_z.max()), height, density))
