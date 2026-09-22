"""Prepare five WP2 pilots without admitting jobs, starting CMW, or running VASP.

Execution uses the local calculation folders, including their Git-ignored POTCAR,
following the WP1 convention. POTCAR is never returned as a publishable artifact.
Optional --stage-private retains the earlier private-copy preparation path.
The frozen launcher owns exclusive canonical runtime creation.
"""

import argparse
import csv
import hashlib
import io
import json
import math
import re
from pathlib import Path

import numpy as np
from pymatgen.io.vasp import Incar, Poscar

WP2 = Path(__file__).resolve().parents[1]
REPO = WP2.parents[1]
PORT = Path('/Users/liangze/Desktop/squiddy tools/vasp-6.6.1-macos-arm64-port')
CMW = Path('/Users/liangze/Desktop/squiddy tools/computational-modelling-workflow')
LEGACY_STATE = Path('/Users/liangze/Library/Application Support/cmw/wp1-IIb-05-01')
STATE = Path('/Users/liangze/Library/Application Support/cmw/wp2-stage01')
SCRATCH = Path('/Volumes/Scratch/wp2_surface_facet_screening_runs/01_runtime_pilot')
PRIVATE = PORT / 'private/wp2-stage01-inputs'
CALC = Path('calculation/01_runtime_pilot')
RESULTS = Path('results/01_runtime_pilot')
POTCAR_DONOR = REPO / 'vasp/wp1_polymorph_polytype_benchmark/calculation/01_geometry_optimisation/beta/POTCAR'
SOURCE = {
    'A': ('beta/001/beta_001_t07', '4904ae4438034ce96765d8bc49ffa974e66e063fcda4b345a784cff3ed386cf2'),
    'B': ('beta/100/beta_100_t02', 'a9093aeae6649f41fcf8c675d227b58b03419ccc8e20124dbf772034a5ff5f47'),
}
JOBS = (
    ('A0', 'beta001_kpoints-static-no-dipole', '01_beta_basal/00_static_no_dipole', False, False),
    ('A1', 'beta001_kpoints-static-dipole', '01_beta_basal/01_static_dipole', True, False),
    ('B0', 'beta100-static-dipole', '02_beta_edge/00_static_dipole', True, False),
    ('A2', 'beta001_kpoints-relax', '01_beta_basal/02_relaxation', True, True),
    ('B1', 'beta100-relax', '02_beta_edge/01_relaxation', True, True),
)
BASELINE = dict(GGA='PE', IVDW=12, ENCUT=500, PREC='Accurate', LASPH=True,
                LREAL=False, ISYM=0, ISPIN=1, NCORE=4, KPAR=1, EDIFF=1e-6,
                NELM=150, NELMIN=2, ALGO='Normal', ISMEAR=0, SIGMA=0.05,
                ISTART=0, ICHARG=2, LCHARG=True, LWAVE=False, LVHAR=True)
BLOCK_REASON = ('Prepared specification; actual queue membership is recorded in runtime_pilot_registry.csv. '
                'Admission requires verified retirement of the preserved WP1 state; relaxations require human release.')


def digest(data):
    return hashlib.sha256(data).hexdigest()


def json_bytes(data):
    return (json.dumps(data, indent=2, sort_keys=True, allow_nan=False) + '\n').encode()


def mesh_record(structure):
    """Use ceil(|b_i|/0.20) with reciprocal vectors including 2*pi, and N3=1."""
    lengths = np.linalg.norm(structure.lattice.reciprocal_lattice.matrix, axis=1)
    mesh = [int(math.ceil(lengths[0] / 0.20)), int(math.ceil(lengths[1] / 0.20)), 1]
    spacings = [float(lengths[i] / mesh[i]) for i in (0, 1)]
    return dict(mesh=mesh, reciprocal_lengths_inv_A=lengths.tolist(),
                actual_inplane_spacings_inv_A=spacings, maximum_inplane_spacing_inv_A=max(spacings),
                convention='b includes 2*pi; N_i=ceil(|b_i|/0.20) for i=1,2; Gamma centered; N3=1',
                converged=False)


def source_structure(family):
    """Read pinned Stage 0 bytes and verify the centered, unconstrained slab."""
    directory, expected = SOURCE[family]
    path = WP2 / 'structure/slab_candidates' / directory / 'POSCAR'
    data = path.read_bytes()
    metadata = json.loads(path.with_name('metadata.json').read_text())
    if digest(data) != expected or metadata['structure_sha256'] != expected:
        raise ValueError(f'{family}: Stage 0 identity changed')
    poscar = Poscar.from_str(data.decode())
    structure = poscar.structure
    z = structure.cart_coords[:, 2]
    if (poscar.selective_dynamics is not None or metadata['diagnostics']['validation_errors']
            or not np.allclose(structure.lattice.matrix[:2, 2], 0, atol=1e-10)
            or not np.allclose(structure.lattice.matrix[2, :2], 0, atol=1e-10)
            or abs((z.min() + z.max()) / structure.lattice.c - 1) > 1e-10):
        raise ValueError(f'{family}: centered Stage 0 seed contract failed')
    # Repeated POSCAR species blocks are intentional; never reorder the seed.
    blocks = data.decode().splitlines()[5].split()
    return path, data, structure, blocks


def build_plan():
    """Produce public inputs and explicit future CMW argv, with no queue mutation."""
    files, specs, rows = {}, [], []
    for short, suffix, subdir, dipole, relaxation in JOBS:
        source, poscar, structure, blocks = source_structure(short[0])
        mesh = mesh_record(structure)
        incar = dict(BASELINE, SYSTEM=f'WP2-P01-{short}-{suffix}',
                     IBRION=2 if relaxation else -1, NSW=100 if relaxation else 0, ISIF=2,
                     LDIPOL=dipole)
        if dipole:
            incar.update(IDIPOL=3, DIPOL=[0.5, 0.5, 0.5])
        if relaxation:
            incar['EDIFFG'] = -0.02
        input_dir = WP2 / CALC / subdir
        private_dir = PRIVATE / subdir
        runtime = SCRATCH / subdir
        timeout = 172800 if relaxation else 21600
        argv = [str(PORT / 'scripts/run-vasp.sh'), '--input', str(input_dir), '--output', str(runtime),
                '--binary', 'std', '--ranks', '8', '--ncore', '4', '--kpar', '1',
                '--mpi-mode', 'synthetic', '--restart', 'none', '--timeout', str(timeout),
                '--stop-before', '0', '--managed-foreground']
        environment = dict(PYTHONPATH=str(CMW / 'src'), CMW_MANAGED_PYTHON=str(CMW / '.venv/bin/python'),
                           CMW_JOBS_OWN_SESSION='1', CMW_JOBS_STATE=str(STATE),
                           OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1')
        name = f'WP2-P01-{short}-{suffix}'
        cmw_add = ['jobs', '--state', str(STATE), 'add', '--name', name, '--engine', 'VASP',
                   '--cwd', str(input_dir), '--cpus', '8', '--mpi-ranks', '8', '--threads-per-rank', '1',
                   '--on-failure', 'pause', '--role', 'primary', '--write-scope', str(SCRATCH), '--json']
        if relaxation:
            cmw_add += ['--hold']
        for key, value in environment.items():
            cmw_add += ['--env', f'{key}={value}']
        cmw_add += ['--', *argv]
        spec = dict(pilot_job=short, name=name, surface_id=source.parent.name,
                    source_structure=str(source.relative_to(REPO)), source_structure_hash=digest(poscar),
                    input_directory=str(input_dir.relative_to(REPO)), private_input_directory=str(private_dir),
                    runtime_root=str(runtime), launcher_argv=argv, cmw_state=str(STATE), cmw_add_arguments=cmw_add,
                    environment=environment, scheduling_role='primary', allow_auxiliary=False, write_scope=str(SCRATCH),
                    ready_or_hold='HOLD' if relaxation else 'READY',
                    requires_explicit_human_release=relaxation, predecessor={'A0':'verified retirement of WP1 queue','A1':'A0','B0':'A1','A2':'human review A0/A1','B1':'human review B0; after A2 if both released'}[short],
                    mpi_ranks=8, ncore=4, kpar=1, timeout_seconds=timeout, stop_before_seconds=0,
                    kpoints=mesh, incar=incar, poscar_species_blocks=blocks,
                    potcar_dataset_order=[{'Zn':'Zn','In':'In_d','S':'S'}[e] for e in blocks],
                    potcar_family='PAW_PBE; WP1 datasets dated 06Sep2000; local Git-ignored POTCAR',
                    dipole_center_convention='Direct lattice coordinates; geometric z midpoint verified 0.5; only z is physically relevant for IDIPOL=3',
                    dispatch_status='PENDING_HUMAN_RELEASE' if relaxation else 'PREPARED',
                    admission_blocker=BLOCK_REASON)
        inputs = {'POSCAR':poscar, 'INCAR':str(Incar(incar)).encode(),
                  'KPOINTS':(f'{name}; pilot spacing <=0.20 inverse A; not converged\n0\nGamma\n'
                             + ' '.join(map(str,mesh['mesh'])) + '\n0 0 0\n').encode()}
        spec['input_sha256'] = {key:digest(value) for key,value in inputs.items()}
        for key,value in inputs.items():
            files[CALC / subdir / key] = value
        files[CALC / subdir / 'job.json'] = json_bytes(spec)
        validate_inputs(spec, inputs)
        specs.append(spec)
        rows.append(dict(pilot_job=name, surface_id=spec['surface_id'], task='relaxation' if relaxation else 'static SCF',
                         dipole_correction=dipole, input_directory=spec['input_directory'], runtime_root=str(runtime),
                         cmw_job_id='', cmw_state=str(STATE), dispatch_status=spec['dispatch_status'],
                         scientific_status='HOLD' if relaxation else 'PREPARED', mpi_ranks=8, ncore=4, kpar=1,
                         k_mesh='x'.join(map(str,mesh['mesh'])), encut=500, ediff=1e-6,
                         ediffg=-0.02 if relaxation else '', nsw=100 if relaxation else 0, timeout=timeout,
                         source_structure_hash=digest(poscar), start_time='',end_time='',walltime='',
                         electronic_iterations='',ionic_steps='',notes=BLOCK_REASON))
    stream = io.StringIO(newline='')
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
    writer.writeheader()
    writer.writerows(rows)
    files[RESULTS / 'runtime_pilot_registry.csv'] = stream.getvalue().encode()
    return files, specs


def validate_inputs(spec, inputs):
    """Fail closed on seed/order, control-pair, sampling and electronic-setting drift."""
    short = spec['pilot_job']
    _, source, _, blocks = source_structure(short[0])
    if inputs['POSCAR'] != source or digest(source) != spec['source_structure_hash']:
        raise ValueError('POSCAR differs from reviewed Stage 0 bytes')
    if spec['potcar_dataset_order'] != [{'S':'S','In':'In_d','Zn':'Zn'}[e] for e in blocks]:
        raise ValueError('POTCAR block order mismatch')
    actual = dict(Incar.from_str(inputs['INCAR'].decode()))
    if actual != dict(Incar(spec['incar'])):
        raise ValueError('INCAR differs from reviewed specification')
    forbidden = {'LHFCALC','AEXX','HFSCREEN','METAGGA','LDAU','LDAUU','LSORBIT','LNONCOLLINEAR','NELECT','NPAR'}
    if forbidden.intersection(actual) or any(actual.get(k) != v for k,v in dict(Incar(BASELINE)).items()):
        raise ValueError('Unexpected scientific or parallel settings')
    mesh = list(map(int, inputs['KPOINTS'].decode().splitlines()[3].split()))
    if mesh != spec['kpoints']['mesh'] or mesh[2] != 1:
        raise ValueError('Unexpected slab sampling')
    if short=='A0' and (actual['LDIPOL'] or 'IDIPOL' in actual or 'DIPOL' in actual):
        raise ValueError('A0 is not an uncorrected control')
    if short!='A0' and (not actual['LDIPOL'] or actual['IDIPOL']!=3 or actual['DIPOL']!=[.5,.5,.5]):
        raise ValueError('Incorrect slab dipole correction')
    relaxation = short in ('A2','B1')
    if actual['NSW']!=(100 if relaxation else 0) or actual['IBRION']!=(2 if relaxation else -1) or actual['ISIF']!=2:
        raise ValueError('Incorrect ionic policy')
    if relaxation and (actual['EDIFFG'] != -.02 or '--hold' not in spec['cmw_add_arguments']):
        raise ValueError('Relaxation force target or human hold missing')


def write_unchanged_or_new(path, data):
    """Never silently replace an existing prepared artifact."""
    if path.exists():
        if path.read_bytes() != data:
            raise FileExistsError(f'Existing prepared artifact differs: {path}')
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('xb') as stream:
        stream.write(data)


def stage_private(specs):
    """Assemble regular private input files in the POSCAR's exact species-block order.

    Dataset bytes are never returned, printed, or written into the scientific
    repository. Complete dataset blocks retain their original bytes.
    """
    payload = POTCAR_DONOR.read_bytes()
    pieces = re.findall(rb'.*?End of Dataset[^\n]*\n', payload, flags=re.S)
    if b''.join(pieces).strip() != payload.strip() or len(pieces)!=3:
        raise ValueError('Cannot resolve the three complete WP1 PAW datasets')
    datasets = {}
    for piece in pieces:
        match = re.search(rb'TITEL\s*=\s*PAW_PBE\s+(\S+)\s+06Sep2000', piece)
        if not match:
            raise ValueError('Unexpected PAW family/date')
        datasets[match.group(1).decode()] = piece
    if set(datasets) != {'S','In_d','Zn'}:
        raise ValueError('Unexpected PAW variants')
    for spec in specs:
        source = REPO / spec['input_directory']
        inputs = {n:(source/n).read_bytes() for n in ('INCAR','POSCAR','KPOINTS')}
        validate_inputs(spec, inputs)
        private = Path(spec['private_input_directory'])
        if private.is_symlink() or not private.resolve().is_relative_to((PORT/'private').resolve()):
            raise ValueError('Private staging escaped the ignored port tree')
        for name,data in inputs.items():
            write_unchanged_or_new(private/name,data)
        write_unchanged_or_new(private/'POTCAR',b''.join(datasets[name] for name in spec['potcar_dataset_order']))
        # The launcher must create the final leaf exclusively at execution time.
        runtime = Path(spec['runtime_root'])
        if runtime.exists() or runtime.is_symlink():
            raise FileExistsError(f'Canonical runtime already exists: {runtime}')
        runtime.parent.mkdir(parents=True,exist_ok=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage-private',action='store_true',help='Stage licensed data only in ignored private inputs; no execution')
    args=parser.parse_args()
    files,specs=build_plan()
    for path,data in files.items():
        target=WP2/path
        if target.exists() and target.read_bytes()!=data:
            raise FileExistsError(f'Existing prepared artifact differs: {target}')
    for path,data in files.items():
        write_unchanged_or_new(WP2/path,data)
    if args.stage_private:
        stage_private(specs)
    for spec in specs:
        print(spec['pilot_job'],spec['ready_or_hold'],spec['kpoints']['mesh'],spec['runtime_root'])
    print('Five specifications prepared; no CMW mutation, controller start, or VASP execution.')


if __name__=='__main__':
    main()
