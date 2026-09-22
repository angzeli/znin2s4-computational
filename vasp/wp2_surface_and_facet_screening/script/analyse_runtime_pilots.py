"""Read later WP2 outputs into bounded SCF/electrostatic diagnostics.

This tool never runs VASP, changes CMW state, releases a relaxation, or edits the
preparation registry. Potential slopes are diagnostics, not automatic approval.
No surface energy, work function, DOS/PDOS, or production band gap is inferred.
"""

import argparse
import csv
import re
from pathlib import Path

import numpy as np
from pymatgen.io.vasp import Chgcar, Incar, Locpot, Poscar, Vasprun
from pymatgen.io.vasp.outputs import Eigenval

from prepare_runtime_pilots import JOBS, WP2, SOURCE, digest
from electrostatics import (fit_window, vacuum_windows, same_geometry, planar_profile,
                            vacuum_window_sensitivity)
from pilot_frontier import electronic_diagnostic, crosscheck_eigenval_xml, OCCUPATION_CONVENTION


def result_states(*, normal, electronic_converged, timeout, fatal, nsw, ibrion,
                  ediffg, forces, ionic_marker):
    """Separate termination and the executed force target; static forces are diagnostic.

    Forces are (ionic step, atom, Cartesian component), in eV/A. A negative
    EDIFFG and the explicit ionic marker are both required for force convergence.
    No final or minimum-force geometry is substituted for another step.
    """
    result = dict(program_termination='NORMAL' if normal else 'UNCONFIRMED',
                  electronic_convergence_status='CONVERGED_ALL_STEPS' if electronic_converged else 'UNRESOLVED',
                  ionic_convergence_status='UNRESOLVED', force_converged=None,
                  geometry_acceptance_status='UNRESOLVED', force_target_eV_per_A=None)
    if fatal or timeout:
        result['program_termination'] = 'NUMERICAL_FAILURE' if fatal else 'TIMEOUT_EVIDENCE'
    if nsw == 0 and ibrion == -1:
        result.update(ionic_convergence_status='NOT_APPLICABLE',
                      geometry_acceptance_status='STATIC_CONTROL_ONLY')
    else:
        vectors = np.asarray(forces, dtype=float)
        if (vectors.ndim != 3 or vectors.shape[-1] != 3 or not vectors.size
                or not np.isfinite(vectors).all() or ediffg is None
                or not np.isfinite(ediffg) or ediffg >= 0):
            return result
        maxima = np.linalg.norm(vectors, axis=2).max(axis=1)
        target = abs(ediffg)
        converged = bool(normal and electronic_converged and not timeout and not fatal
                         and ionic_marker and maxima[-1] < target)
        result.update(force_converged=converged, force_target_eV_per_A=target,
                      minimum_max_force_eV_per_A=float(maxima.min()),
                      minimum_force_step=int(np.argmin(maxima))+1,
                      ionic_convergence_status='FORCE_CONVERGED' if converged else
                      'NSW_LIMIT_FORCE_UNCONVERGED' if len(maxima) >= nsw and not ionic_marker else 'FORCE_UNCONVERGED',
                      geometry_acceptance_status='FORCE_CONVERGED_PILOT' if converged
                      else 'CANDIDATE_CONTINUATION_REVIEW_REQUIRED')
    if not normal or not electronic_converged or timeout or fatal:
        result['geometry_acceptance_status'] = 'UNRESOLVED'
    return result


def final_frontier(directory, run, evidence):
    """Read EIGENVAL only after final XML/geometry/SCF identity is established."""
    unavailable = electronic_diagnostic({})
    try:
        if not (evidence['normal_termination'] and evidence['electronic_converged']
                and evidence.get('geometry_identity_verified') is True
                and not evidence['timeout_evidence'] and not evidence['fatal_evidence']):
            raise ValueError('Completed final electronic solution at the intended geometry is unverified')
        params = run.parameters
        if run.vasp_version != '6.6.1':
            raise ValueError('Occupation convention is audited only for these VASP 6.6.1 outputs')
        eigenval = Eigenval(directory/'EIGENVAL')
        if len(eigenval.eigenvalues) != 1 or len(run.eigenvalues) != 1:
            raise ValueError('Only one ISPIN=1 channel is supported')
        values = next(iter(eigenval.eigenvalues.values()))
        xml = next(iter(run.eigenvalues.values()))
        # pymatgen ignores printed band indices: verify them before relying on N.
        lines = [line.split() for line in (directory/'EIGENVAL').read_text().splitlines()[6:] if line.strip()]
        nk, nb, _ = values.shape
        if len(lines) != nk*(nb+1) or any(
                int(lines[k*(nb+1)+band][0]) != band for k in range(nk) for band in range(1, nb+1)):
            raise ValueError('EIGENVAL printed band indices or block count disagree')
        if eigenval.ispin != params.get('ISPIN') or eigenval.nelect != params['NELECT']:
            raise ValueError('EIGENVAL/XML electron count or spin metadata disagree')
        result = electronic_diagnostic(eigenval.eigenvalues, nelect=params['NELECT'],
            ispin=params.get('ISPIN'), nbands=params.get('NBANDS'),
            kpoints=eigenval.kpoints, weights=eigenval.kpoints_weights,
            noncollinear=params.get('LNONCOLLINEAR'), soc=params.get('LSORBIT'),
            occupation_convention=OCCUPATION_CONVENTION)
        if result['validation_status'] != 'VALIDATED':
            return result
        result.update(crosscheck_eigenval_xml(values, xml, eigenval.kpoints, run.actual_kpoints,
                      eigenval.kpoints_weights, run.actual_kpoints_weights, params['NELECT']))
        result.update(source_eigenval=str(directory/'EIGENVAL'), source_xml=str(directory/'vasprun.xml'),
                      occupation_convention_evidence='Unscaled EIGENVAL/XML arrays agree; both weighted electron counts validated')
        return result
    except (OSError, ValueError, KeyError, IndexError, TypeError) as exc:
        return dict(unavailable, reason=str(exc))


def output_evidence(directory):
    """Separate normal termination, electronic convergence, timeout, and ionic progress."""
    directory=Path(directory)
    if not directory.exists():
        return dict(scientific_status='PREPARED',reason='Canonical runtime does not exist'),None
    outcar=(directory/'OUTCAR').read_text(errors='replace') if (directory/'OUTCAR').exists() else ''
    normal='General timing and accounting informations for this job:' in outcar
    marker='aborting loop because EDIFF is reached' in outcar
    timing=re.findall(r'Elapsed time \(sec\):\s*([-+\d.Ee]+)',outcar)
    loop_times=[float(x) for x in re.findall(r'LOOP:\s+cpu time\s+\S+:?\s+real time\s+([\d.Ee+-]+)',outcar)]
    ionic_times=[float(x) for x in re.findall(r'LOOP\+:\s+cpu time\s+\S+:?\s+real time\s+([\d.Ee+-]+)',outcar)]
    metadata=(directory/'RUN_METADATA.txt').read_text(errors='replace') if (directory/'RUN_METADATA.txt').exists() else ''
    timeout=bool(re.search(r'(?:reason|termination)\s*[:=]\s*timeout|(?:exit_code|child_status|status)\s*[:=]\s*124',metadata,re.I))
    fatal=bool(re.search(r'BRMIX: very serious problems|Error EDDDAV|ZBRENT: fatal error|VERY BAD NEWS',outcar))
    result=dict(normal_termination=normal,ediff_marker=marker,timeout_evidence=timeout,fatal_evidence=fatal,
                walltime_s=float(timing[-1]) if timing else None,
                mean_electronic_step_walltime_s=float(np.mean(loop_times)) if loop_times else None,
                ionic_step_walltimes_s=ionic_times,scientific_status='INCONCLUSIVE',
                reason='Incomplete or unconfirmed scientific evidence; no automatic retry/release',
                start_time=None,end_time=None)
    # Managed hard deadlines may terminate the launcher before metadata finalization.
    # Missing metadata never becomes a fabricated timeout or a numerical failure.
    for key in ('start_utc','end_utc'):
        match=re.findall(r'\b'+key+r'\s*[:=]\s*(\S+)',metadata)
        result['start_time' if key=='start_utc' else 'end_time']=match[-1] if match else None
    try:
        run=Vasprun(directory/'vasprun.xml',parse_dos=False,parse_eigen=True,
                    parse_projected_eigen=False,parse_potcar_file=False,exception_on_bad_xml=True)
    except Exception as exc:
        result['xml_status']='UNAVAILABLE_OR_INCOMPLETE: '+str(exc)
        if fatal:
            result.update(scientific_status='FAILED_NUMERICAL',reason='Explicit numerical-failure evidence; XML incomplete')
        return result,None
    steps=run.ionic_steps
    iterations=[len(s['electronic_steps']) for s in steps]
    # Map each EDIFF marker to its ionic step instead of inferring all-step success
    # from the last SCF loop or from one marker somewhere in OUTCAR.
    ediff_steps = []
    ionic_index = None
    for match in re.finditer(r'Iteration\s+(\d+)\(\s*\d+\)|aborting loop because EDIFF is reached', outcar):
        if match.group(1):
            ionic_index = int(match.group(1))
        else:
            ediff_steps.append(ionic_index)
    converged=bool(run.converged_electronic and iterations
                   and all(2<=n<int(run.parameters.get('NELM',60)) for n in iterations)
                   and ediff_steps == list(range(1,len(steps)+1)))
    forces=[float(np.max(np.linalg.norm(s['forces'],axis=1))) for s in steps if 'forces' in s]
    displacements=[]
    previous=run.initial_structure
    for step in steps:
        structure=step.get('structure')
        if structure is None:
            displacements.append(None)
            continue
        # Minimum-image displacement is a bounded trend; it does not classify a reconstruction.
        delta=structure.frac_coords-previous.frac_coords
        delta-=np.rint(delta)
        displacements.append(float(np.max(np.linalg.norm(delta@structure.lattice.matrix,axis=1))))
        previous=structure
    result.update(xml_status='COMPLETE',electronic_converged=converged,
                  electronic_iterations=iterations[-1] if iterations else None,
                  electronic_iterations_per_ionic_step=iterations,ionic_steps=len(steps),
                  final_energy_eV=float(run.final_energy),nelect=float(run.parameters['NELECT']),
                  max_force_eV_per_A_by_step=forces,first_five_max_forces_eV_per_A=forces[:5],
                  max_displacement_A_by_ionic_step=displacements,ediff_step_indices=ediff_steps)
    try:
        initial = Poscar.from_file(directory/'POSCAR').structure
        final = Poscar.from_file(directory/'CONTCAR').structure
        result['geometry_identity_verified'] = (same_geometry(initial,run.initial_structure)
            and same_geometry(final,run.final_structure) and same_geometry(final,steps[-1]['structure']))
        executed_nelect = re.findall(r'NELECT\s*=\s*([-+\d.Ee]+)',outcar)
        if not executed_nelect or abs(float(executed_nelect[-1])-result['nelect'])>1e-6:
            raise ValueError('OUTCAR/XML NELECT mismatch or missing executed electron count')
        incar = Incar.from_file(directory/'INCAR')
        # These launcher parallel tags are absent from the audited XML INCAR;
        # confirm their recorded OUTCAR values rather than guessing defaults.
        for key,value in incar.items():
            if key in ('NCORE','KPAR') and key not in run.incar:
                found = re.findall(r'\b'+key+r'\s*=\s*(\d+)',outcar)
                if not found or int(found[-1]) != value:
                    raise ValueError(key+' disagrees with the executed OUTCAR')
            elif run.incar.get(key) != value:
                raise ValueError('Archived INCAR differs from executed XML: '+key)
        result['input_settings_verified'] = True
    except (OSError,ValueError,KeyError,IndexError) as exc:
        result.update(geometry_identity_verified=False,input_settings_verified=False,provenance_reason=str(exc))
    result.update(result_states(normal=normal,electronic_converged=converged,timeout=timeout,fatal=fatal,
                  nsw=run.parameters.get('NSW'),ibrion=run.parameters.get('IBRION'),
                  ediffg=run.parameters.get('EDIFFG'),forces=[s.get('forces') for s in steps],
                  ionic_marker='reached required accuracy - stopping structural energy minimisation' in outcar))
    result['electronic_diagnostic'] = final_frontier(directory,run,result)
    ionic_ok=result['ionic_convergence_status'] in ('NOT_APPLICABLE','FORCE_CONVERGED')
    if fatal:
        result.update(scientific_status='FAILED_NUMERICAL',reason='Explicit numerical-failure marker')
    elif timeout:
        result.update(scientific_status='INCONCLUSIVE',reason='Safety timeout; resource review, not structural invalidity')
    elif normal and converged and ionic_ok and result['geometry_identity_verified'] and run.vasp_version=='6.6.1':
        result.update(scientific_status='PASS',reason='Numerical completion criteria met; electrostatics and relaxation release still require human review')
    elif result['ionic_convergence_status']=='NSW_LIMIT_FORCE_UNCONVERGED':
        result.update(scientific_status='NSW_LIMIT_FORCE_UNCONVERGED',reason='Normal NSW-limited endpoint; force convergence and structural review remain pending')
    return result,run


def analyse_job(directory, spec=None):
    """Read final potential/density profiles only when their geometry and output type agree."""
    directory=Path(directory)
    if directory.exists() and spec is not None:
        try:
            if digest((directory/'POSCAR').read_bytes())!=spec['source_structure_hash']:
                raise ValueError('Runtime input geometry is not the selected Stage 0 seed')
            if 'incar' in spec and dict(Incar.from_file(directory/'INCAR'))!=dict(Incar(spec['incar'])):
                raise ValueError('Runtime scientific settings differ from the prepared pilot')
        except (OSError,ValueError) as exc:
            return dict(scientific_status='HOLD',reason='Input provenance mismatch: '+str(exc)),None
    result,run=output_evidence(directory)
    if run is None:
        return result,None
    if not (result.get('normal_termination') and result.get('electronic_converged')
            and result.get('geometry_identity_verified')):
        result['potential_status']='UNAVAILABLE: completed final solution/geometry unverified'
        return result,None
    incar=Incar.from_file(directory/'INCAR')
    if not incar.get('LVHAR',False):
        result['potential_status']='UNAVAILABLE: LVHAR is not enabled'
        return result,None
    try:
        locpot=Locpot.from_file(directory/'LOCPOT')
        if not same_geometry(locpot.structure,run.final_structure):
            raise ValueError('LOCPOT geometry differs from final XML')
        charge=None
        try:
            charge=Chgcar.from_file(directory/'CHGCAR')
        except Exception as exc:
            result['charge_reader_reason']=str(exc)
        data=planar_profile(locpot,charge=charge,nelect=result['nelect'])
        z,potential,density=data['z_A'],data['lvhar_potential_eV'],data['density_e_per_A3']
        height=float(locpot.structure.lattice.c)
        result.update(charge_density_status=data['charge_density_status'],
                      charge_integral_electrons=data['charge_integral_electrons'],vacuum=data['vacuum'])
        atom_z=locpot.structure.cart_coords[:,2]
        result['vacuum_window_checks']=vacuum_window_sensitivity(z,potential,float(atom_z.min()),float(atom_z.max()),height,density)
        result['potential_status']='AVAILABLE_LVHAR'
        profile=[dict(z_A=float(x),lvhar_potential_eV=float(v),
                      mean_electron_density_e_per_A3=float(density[i]) if density is not None else None)
                 for i,(x,v) in enumerate(zip(z,potential))]
        return result,profile
    except Exception as exc:
        result['potential_status']='UNAVAILABLE_OR_AMBIGUOUS: '+str(exc)
        return result,None


def pilot_summary(short, directory, result):
    """Flatten observed axes and reconcile only the five already-audited pilot roles."""
    surface = 'beta_001_t07' if short.startswith('A') else 'beta_100_t02'
    roles = dict(A0='ACCEPTED_UNCORRECTED_STATIC_CONTROL', A1='ACCEPTED_CORRECTED_BASAL_STATIC_CONTROL',
                 B0='ACCEPTED_UNRELAXED_EDGE_STATIC_CONTROL', A2='FORCE_CONVERGED_SINGLE_LAYER_SEED',
                 B1='CANDIDATE_CONTINUATION_GEOMETRY_REQUIRING_REVIEW')
    complete = (result.get('normal_termination') and result.get('electronic_converged')
                and result.get('geometry_identity_verified') and not result.get('timeout_evidence')
                and not result.get('fatal_evidence'))
    role = roles[short] if complete else 'UNRESOLVED'
    if short=='A2' and result.get('force_converged') is not True:
        role='UNRESOLVED'
    if short=='B1' and result.get('ionic_convergence_status')!='NSW_LIMIT_FORCE_UNCONVERGED':
        role='REQUIRES_NEW_REVIEW'
    directory=Path(directory).resolve()
    source=str(directory.relative_to(WP2)) if directory.is_relative_to(WP2) else str(directory)
    row=dict(pilot=short,surface_id=surface,lower_face_id=surface+'_lower',upper_face_id=surface+'_upper',
             source_directory=source,source_eigenval=source+'/EIGENVAL',source_xml=source+'/vasprun.xml',
             source_structure_sha256=SOURCE[short[0]][1],accepted_role=role,
             face_identity_note='Inherited face labels; no equivalence or endpoint coordination claim',
             production_surface_gap_status='NOT_ESTABLISHED')
    for key in ('normal_termination','program_termination','electronic_converged','electronic_convergence_status',
                'ionic_convergence_status','force_converged','force_target_eV_per_A','geometry_acceptance_status',
                'geometry_identity_verified','input_settings_verified','timeout_evidence','ionic_steps','walltime_s',
                'nelect','final_energy_eV','minimum_max_force_eV_per_A','minimum_force_step',
                'charge_density_status','charge_integral_electrons'):
        row[key]=result.get(key)
    maxima=result.get('max_force_eV_per_A_by_step',[])
    row['final_max_force_eV_per_A']=maxima[-1] if maxima else None
    row['total_scf_iterations']=sum(result.get('electronic_iterations_per_ionic_step',[]))
    diagnostic=result.get('electronic_diagnostic',electronic_diagnostic({}))
    row.update({k:v for k,v in diagnostic.items() if k not in ('source_eigenval','source_xml')})
    checks=result.get('vacuum_window_checks',[])
    distant=[r for r in checks if r['setback_A']==6]
    row['automatic_lower_status']=result.get('vacuum',{}).get('lower',{}).get('status')
    row['automatic_upper_status']=result.get('vacuum',{}).get('upper',{}).get('status')
    passing=len(distant)==2 and all(r['status']=='CANDIDATE_PLATEAU' for r in distant)
    row['distant_window_status']='CANDIDATE_PLATEAUS_BOTH_FACES' if passing else 'NO_CONFIRMED_FIELD_FREE_REFERENCE'
    row['distant_upper_minus_lower_potential_eV']=(next(r['mean_potential_eV'] for r in distant if r['side']=='upper')
        -next(r['mean_potential_eV'] for r in distant if r['side']=='lower')) if passing else None
    windows=[dict(pilot=short,surface_id=surface,face_id=surface+'_'+r['side'],**r) for r in checks]
    return row,windows


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runtime-root',type=Path,default=WP2/'calculation/01_runtime_pilot',
                        help='Archived five-pilot root; no queue/controller is contacted')
    parser.add_argument('--output-directory',type=Path,required=True,help='New analysis directory; never overwrites prepared registry or audit')
    args=parser.parse_args()
    results,windows=[],[]
    for short,_,subdir,_,_ in JOBS:
        directory=args.runtime_root/subdir
        result,_=analyse_job(directory,{'source_structure_hash':SOURCE[short[0]][1]})
        row,checks=pilot_summary(short,directory,result)
        results.append(row);windows.extend(checks)
    args.output_directory.mkdir(parents=True,exist_ok=False)
    for name,rows in [('runtime_pilot_summary.csv',results),('vacuum_window_checks.csv',windows)]:
        fields=list(dict.fromkeys(key for row in rows for key in row))
        with (args.output_directory/name).open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=fields)
            writer.writeheader();writer.writerows(rows)


if __name__=='__main__':
    main()
