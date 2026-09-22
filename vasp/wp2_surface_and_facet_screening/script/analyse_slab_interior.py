"""Reproduce the endpoint interior audit; never write to calculation directories.

Run with the repository's .venv-wp1-plots Python. --check-only validates and
prints measurements without exporting profiles or figures. All potentials are
LVHAR Hartree-plus-ionic energies in eV; positions and window lengths are in Å.
Exports default to the numerical-convergence results directory; figures are
optional and require --figures. Use --output-directory for regeneration checks.
"""
import argparse
import csv
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True
import numpy as np
from pymatgen.core import Structure
from pymatgen.io.vasp.outputs import Chgcar, Locpot

HERE = Path(__file__).resolve().parent
WP2 = HERE.parent
RESULTS = WP2 / 'results/02_numerical_convergence'
CALC = WP2 / 'calculation/02_numerical_convergence/03_slab_thickness/01_relaxation'
CASES = {2: ('beta001_2L', 'J33.1'), 4: ('beta001_4L', 'J34.1'),
         6: ('beta001_6L/retry_nelm300_r1', 'J36.1')}
# Original construction repeat, measured from equivalent sites in the inputs.
REPEAT = 12.14843313327493
U = np.linspace(-REPEAT / 4, REPEAT / 4, 2001)


class PeriodicLinear:
    """Exact integrals of the periodic piecewise-linear native-grid profile."""
    def __init__(self, values, height):
        self.height = height
        self.x = np.linspace(0., height, len(values) + 1)
        self.y = np.r_[values, values[0]]
        self.dx = height / len(values)
        self.cumulative = np.r_[0., np.cumsum((self.y[:-1] + self.y[1:]) * self.dx / 2)]

    def primitive(self, z):
        z = np.asarray(z)
        periods = np.floor(z / self.height)
        local = z - periods * self.height
        i = np.minimum((local / self.dx).astype(int), len(self.y) - 2)
        t = local - self.x[i]
        return (periods * self.cumulative[-1] + self.cumulative[i]
                + self.y[i] * t + (self.y[i+1] - self.y[i]) * t*t / (2*self.dx))

    def mean(self, lo, hi):
        return (self.primitive(hi) - self.primitive(lo)) / (np.asarray(hi) - lo)

    def box(self, z, width=REPEAT):
        return self.mean(np.asarray(z) - width/2, np.asarray(z) + width/2)

    def direct_mean(self, lo, hi):
        """Independent trapezoidal quadrature at all enclosed knots (no cumulative sums)."""
        assert 0 <= lo < hi <= self.height
        knots = np.r_[lo, self.x[(self.x > lo) & (self.x < hi)], hi]
        return np.trapezoid(np.interp(knots, self.x, self.y), knots) / (hi - lo)


def same_geometry(a, b):
    assert a.species == b.species
    # LOCPOT/CHGCAR headers print six decimals for both cell and fractions.
    np.testing.assert_allclose(a.lattice.matrix, b.lattice.matrix, atol=5.1e-7, rtol=0)
    delta = a.frac_coords - b.frac_coords
    np.testing.assert_allclose(delta - np.rint(delta), 0, atol=5.1e-7, rtol=0)


def read_case(n, relative, job):
    path = CALC / relative
    out = (path / 'OUTCAR').read_text()
    assert 'reached required accuracy' in out and 'General timing and accounting' in out
    assert 'campaign_cmw_display_id: ' + job in (path / 'RUN_METADATA.txt').read_text()
    last = ET.parse(path / 'vasprun.xml').getroot().findall('calculation')[-1]
    parameters = ET.parse(path / 'vasprun.xml').getroot().find('parameters')
    assert parameters.find(".//i[@name='LVHAR']").text.strip() == 'T'
    final = Structure.from_file(path / 'CONTCAR')
    basis = [[float(t) for t in v.text.split()] for v in last.findall("structure/crystal/varray[@name='basis']/v")]
    frac = [[float(t) for t in v.text.split()] for v in last.findall("structure/varray[@name='positions']/v")]
    same_geometry(final, Structure(basis, final.species, frac))
    loc = Locpot.from_file(path / 'LOCPOT')
    chg = Chgcar.from_file(path / 'CHGCAR')
    same_geometry(final, loc.structure)
    same_geometry(final, chg.structure)
    assert loc.dim == chg.dim and len(final) == 7*n
    height = final.lattice.matrix[2, 2]
    np.testing.assert_allclose(final.lattice.matrix[:2, 2], 0., atol=1e-10)
    np.testing.assert_allclose(final.lattice.matrix[2, :2], 0., atol=1e-10)
    raw = np.asarray(loc.get_average_along_axis(2))
    np.testing.assert_allclose(raw, np.mean(loc.data['total'], axis=(0, 1)), atol=1e-10)
    density = np.asarray(chg.get_average_along_axis(2)) / final.volume
    charge = float(np.mean(chg.data['total']))
    assert abs(charge - 62*n) < 1e-3
    assert np.isfinite(raw).all() and np.isfinite(density).all()
    z = np.arange(len(raw)) * height / len(raw)
    az = final.cart_coords[:, 2]
    assert 0 < az.min() < az.max() < height
    bounds = np.array([[az[7*k:7*k+7].min(), az[7*k:7*k+7].max()] for k in range(n)])
    for k in range(n):
        assert [str(s) for s in final.species[7*k:7*k+7]] == ['S']*4 + ['In']*2 + ['Zn']
    assert np.all(bounds[1:, 0] > bounds[:-1, 1])
    gap_mid = (bounds[:-1, 1] + bounds[1:, 0]) / 2
    register = float(gap_mid[n//2-1])
    slab_center = float((az.min() + az.max()) / 2)
    model = PeriodicLinear(raw, height)
    charge_model = PeriodicLinear(density, height)
    macro = model.box(z)
    # The common half-repeat comparison plus its averaging footprint is wholly
    # inside even the 2L atomic envelope; no vacuum is averaged into the audit.
    assert register + U[0] - REPEAT/2 > az.min()
    assert register + U[-1] + REPEAT/2 < az.max()
    central = model.box(register + U)
    independent = np.array([model.direct_mean(t - REPEAT/2, t + REPEAT/2) for t in register + U])
    np.testing.assert_allclose(central, independent, atol=2e-11, rtol=0)
    np.testing.assert_allclose(model.box(z + height), macro, atol=2e-11, rtol=0)
    slope = float(np.polyfit(U, central, 1)[0])
    np.testing.assert_allclose(slope, np.dot(U, independent-independent.mean()) / np.dot(U, U), atol=1e-12)
    mean = float(np.trapezoid(central, U) / np.ptp(U))
    np.testing.assert_allclose(mean, np.trapezoid(independent, U)/np.ptp(U), atol=2e-11)
    refs = {}
    for side, lo, hi in [('lower', 1., az.min()-6.), ('upper', az.max()+6., height-1.)]:
        mask = (z >= lo) & (z <= hi)
        assert np.ptp(z[mask]) >= 2 and np.ptp(raw[mask]) < .03
        assert np.max(np.abs(density[mask])) < 1e-5
        refs[side] = float(raw[mask].mean())
    lo, hi = bounds[n//2-1, 1], bounds[n//2, 0]
    inter_lo, inter_hi = lo + .25*(hi-lo), hi - .25*(hi-lo)
    gap_potential = float(model.mean(inter_lo, inter_hi))
    gap_density = float(charge_model.mean(inter_lo, inter_hi))
    np.testing.assert_allclose(gap_potential, model.direct_mean(inter_lo, inter_hi), atol=2e-11)
    np.testing.assert_allclose(gap_density, charge_model.direct_mean(inter_lo, inter_hi), atol=2e-11)
    return dict(n=n, path=relative, height=height, z=z, raw=raw, macro=macro, density=density,
                bounds=bounds, register=register, slab_center=slab_center, refs=refs,
                central=central, central_mean=mean, slope=slope, span=float(np.ptp(central)),
                gap_mid=gap_mid, gap_macro=model.box(gap_mid), gap_potential=gap_potential,
                gap_density=gap_density, model=model, integral=charge,
                valid=(z >= az.min()+REPEAT/2) & (z <= az.max()-REPEAT/2))


def describe(cases):
    print('L =', REPEAT, 'Å; central comparison width =', np.ptp(U), 'Å')
    for n, d in cases.items():
        print('\nSLAB', n, 'source', d['path'])
        print('slab_center / central_gap_center / layer_envelopes (Å):', d['slab_center'], d['register'], d['bounds'].tolist())
        print('central mean lower/upper (eV):', *[d['central_mean']-d['refs'][s] for s in ['lower','upper']])
        print('central slope (eV/Å) / range (eV):', d['slope'], d['span'])
        print('gap-centre macro lower-aligned (eV):', (d['gap_macro']-d['refs']['lower']).tolist())
        print('central gap inner-half mean lower/upper (eV) / density (e/Å³):',
              *[d['gap_potential']-d['refs'][s] for s in ['lower','upper']], d['gap_density'])
        if n == 6:
            interior_gap = d['gap_mid'][1:-1]
            independent_gap = np.array([d['model'].direct_mean(t-REPEAT/2, t+REPEAT/2)
                                        for t in interior_gap])
            np.testing.assert_allclose(d['gap_macro'][1:-1], independent_gap, atol=2e-11)
            print('inner three gap-centre span (Å) / potential range (eV):',
                  np.ptp(interior_gap), np.ptp(independent_gap))
    for a, b in [(2, 4), (4, 6)]:
        print('\nCOMPARISON', a, 'to', b)
        for side in ['lower', 'upper']:
            delta = (cases[b]['central']-cases[b]['refs'][side])-(cases[a]['central']-cases[a]['refs'][side])
            print(side, 'central mean delta / max|delta| / RMS delta (eV):',
                  float(np.trapezoid(delta,U)/np.ptp(U)), float(np.abs(delta).max()),
                  float(np.sqrt(np.trapezoid(delta**2,U)/np.ptp(U))))
        print('density relative change:', cases[b]['gap_density']/cases[a]['gap_density']-1)


def export(cases, *, output_directory=RESULTS, figures=False):
    output_directory = Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    with (output_directory/'03_SLAB_INTERIOR_PROFILES.csv').open('w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['layers','z_A','z_from_central_gap_A','planar_LVHAR_eV','macro_LVHAR_eV',
                         'planar_density_e_per_A3','lower_vacuum_eV','upper_vacuum_eV','macro_inside_slab'])
        for n,d in cases.items():
            for i,z in enumerate(d['z']):
                writer.writerow([n, f'{z:.12g}', f'{z-d["register"]:.12g}', f'{d["raw"][i]:.12g}',
                                 f'{d["macro"][i]:.12g}', f'{d["density"][i]:.12g}',
                                 f'{d["refs"]["lower"]:.12g}', f'{d["refs"]["upper"]:.12g}',int(d['valid'][i])])
    if not figures:
        return
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'figure.facecolor':'white','axes.facecolor':'white','savefig.facecolor':'white',
        'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','Liberation Sans','DejaVu Sans'],
        'mathtext.fontset':'stixsans','text.color':'black','axes.edgecolor':'black','axes.labelcolor':'black',
        'axes.labelsize':13,'axes.labelweight':'bold','axes.titlesize':14,'axes.titleweight':'bold',
        'axes.linewidth':1.8,'axes.grid':False,'xtick.labelsize':10,'ytick.labelsize':10,
        'xtick.direction':'in','ytick.direction':'in','xtick.major.width':1.8,'ytick.major.width':1.8,
        'xtick.major.size':4,'ytick.major.size':4,'xtick.top':False,'ytick.right':False,
        'lines.linewidth':2.,'legend.fontsize':9,'legend.frameon':True,'legend.framealpha':1.,
        'legend.facecolor':'white','legend.edgecolor':'black','savefig.dpi':600})
    # Thickness is ordered: light-to-dark shades of the same blue, plus line styles.
    colours={2:'#8DBBD5',4:'#397FAB',6:'#004C73'}
    styles={2:':',4:'--',6:'-'}
    fig=plt.figure(figsize=(8,6),layout='constrained');grid=fig.add_gridspec(2,6)
    axes=[]
    for k,(n,d) in enumerate(cases.items()):
        ax=fig.add_subplot(grid[0,2*k:2*k+2]);axes.append(ax)
        ax.plot(d['z']-d['register'],d['raw']-d['refs']['lower'],color=colours[n],label='Planar')
        ax.plot(d['z'][d['valid']]-d['register'],d['macro'][d['valid']]-d['refs']['lower'],color='black',label='12.148 Å mean')
        for lo,hi in d['bounds']:
            ax.plot([lo-d['register'],hi-d['register']],[-.6,-.6],color='#4D4D4D',lw=3)
        ax.set_xlabel('z − gap centre (Å)');ax.set_title(f'{chr(97+k)}  {n}L',loc='left')
        ax.set_ylabel('V − lower vacuum (eV)' if k==0 else '')
        if k==0:ax.legend(loc='lower right',prop={'size':8,'weight':'bold'})
    for k,side in enumerate(['lower','upper']):
        ax=fig.add_subplot(grid[1,3*k:3*k+3]);axes.append(ax)
        for n,d in cases.items():
            mask=d['valid'];ax.plot(d['z'][mask]-d['register'],d['macro'][mask]-d['refs'][side],
                color=colours[n],ls=styles[n],label=f'{n}L')
        ax.axvspan(U[0],U[-1],color='#EEEEEE',zorder=0)
        ax.set_xlabel('z − gap centre (Å)');ax.set_ylabel('Averaged potential (eV)')
        ax.set_title(f'{chr(100+k)}  {side.capitalize()}-vacuum alignment',loc='left')
        ax.legend(loc='upper left',prop={'size':9,'weight':'bold'})
    for ax in axes:
        ax.grid(False,which='both');ax.minorticks_off()
        for spine in ax.spines.values():spine.set_visible(True);spine.set_linewidth(1.8)
        for label in ax.get_xticklabels()+ax.get_yticklabels():label.set_fontweight('bold')
    for extension in ['png','pdf']:
        fig.savefig(output_directory/f'SLAB_INTERIOR_ELECTROSTATICS.{extension}',bbox_inches='tight',facecolor='white',transparent=False)
    plt.close(fig)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check-only',action='store_true')
    parser.add_argument('--output-directory',type=Path,default=RESULTS)
    parser.add_argument('--figures',action='store_true',help='Also export the optional PNG and PDF figures')
    args=parser.parse_args()
    initial = Structure.from_file(CALC / CASES[2][0] / 'POSCAR')
    np.testing.assert_allclose(initial.cart_coords[7:, 2] - initial.cart_coords[:-7, 2],
                               REPEAT, atol=1e-10, rtol=0)
    cases={n:read_case(n,*spec) for n,spec in CASES.items()}
    describe(cases)
    if not args.check_only:export(cases,output_directory=args.output_directory,figures=args.figures)
