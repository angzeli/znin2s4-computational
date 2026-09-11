"""Execute one of the two plotting notebooks in a fresh current-environment kernel.

Usage: python execute_wp1_plot_notebook.py NOTEBOOK [--cwd DIRECTORY]
Only the named notebook is updated. Its cells own the explicit figure rerun policy.
No kernel installation or package/environment mutation is performed.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import sys
from tempfile import TemporaryDirectory


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('notebook', type=Path)
    parser.add_argument('--cwd', type=Path, help='Execution directory; defaults to the notebook directory')
    args = parser.parse_args()
    path = args.notebook.resolve()
    if path.name not in ('04_01_DOS_PDOS_PLOTTING.ipynb', '04_02_BAND_STRUCTURE_PLOTTING.ipynb'):
        parser.error('Select one of the two WP1 plotting notebooks')
    for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
        os.environ[name] = '1'
    os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
    os.environ['MPLBACKEND'] = 'Agg'
    import nbformat
    from nbclient import NotebookClient
    from jupyter_client import KernelManager
    notebook = nbformat.read(path, as_version=4)
    nbformat.validate(notebook)
    with TemporaryDirectory(prefix='wp1-notebook-') as runtime:
        # Local IPC needs no listening TCP port; all connection state is temporary.
        os.environ['IPYTHONDIR'] = str(Path(runtime) / 'ipython')
        manager = KernelManager(kernel_name='python3', transport='ipc',
                                connection_file=str(Path(runtime) / 'kernel.json'))
        manager.kernel_spec.argv = [sys.executable, '-m', 'ipykernel_launcher', '-f', '{connection_file}']
        client = NotebookClient(notebook, km=manager, timeout=600,
                                resources={'metadata': {'path': str((args.cwd or path.parent).resolve())}})
        try:
            client.execute()
        finally:
            if manager.has_kernel:
                manager.shutdown_kernel(now=True)
            else:
                manager.cleanup_resources()
    # Execution timestamps add noise but no scientific provenance; retain cell counts/results.
    for cell in notebook.cells:
        cell.metadata.pop('execution', None)
    nbformat.validate(notebook)
    nbformat.write(notebook, path)
    print(f'Executed without errors: {path}')


if __name__ == '__main__':
    main()
