"""Make blazar detectabiliity plot."""
from __future__ import print_function
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np

table = Table.read('data/results_cut1.00TeV.fits')
print('Number of sources:', len(table))
table = table[table['Sigma'] >= 5.]
print('Number of detected sources:', len(table))

fig, ax = plt.subplots(figsize=(6, 4))

ax.axhline(y=5, c='green', lw=3)
ax.text(x=0.022, y=460, s='TECHNICAL PLOT', bbox=dict(color='white'),
        fontsize=15, fontweight='bold', color='black', alpha=0.3)

# Draw line to connect the two significance estimates for each sources
ax.plot([table['Redshift'], table['Redshift']],
        [table['Sigma'], table['SigmaCut']],
        c='b', lw=1, alpha=0.4)


def plot_sources(label, source_classes, sigma_col, **opts):
    """Helper function to select and plot sources of a given class"""
    mask = np.array([row['CLASS'] in source_classes for row in table])
    t = table[mask]
    ax.errorbar(
        x=t['Redshift'], y=t[sigma_col], yerr=t[sigma_col + 'Err'],
        label=label, markersize=6, alpha=0.8, **opts
    )


opts = dict(source_classes={'bll', 'BLL'}, c='blue', fmt='o')
plot_sources(label='BL Lac', sigma_col='Sigma', **opts)
plot_sources(label='BL Lac (with cutoff)', sigma_col='SigmaCut', fillstyle='none', **opts)

opts = dict(source_classes={'fsrq', 'FSRQ'}, c='red', fmt='s')
plot_sources(label='FSRQ', sigma_col='Sigma', **opts)
plot_sources(label='FSRQ (with cutoff)', sigma_col='SigmaCut', fillstyle='none', **opts)

ax.set_yscale('log')
ax.set_xscale('log')

ax.set_xlim(0.02, 2)
ax.set_ylim(2, 700)

ax.set_xlabel('Redshift')
ax.set_ylabel('Significance (20 h observation time)')

ax.grid(which='both')
ax.legend(loc='upper right')
fig.tight_layout()

filename = 'plots/cut_off_agn_pop.png'
print('Writing', filename)
fig.savefig(filename, format='png', dpi=300)

filename = 'plots/cut_off_agn_pop.pdf'
print('Writing', filename)
fig.savefig(filename, format='pdf')

plt.show()
