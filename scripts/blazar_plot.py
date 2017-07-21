"""
Small script to plot results
"""
import os
import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.table import Table
import astropy.units as u

import matplotlib.pyplot as plt
from matplotlib import collections  as mc

# Read results
t_orig = Table.read('./data/results_cut1.00TeV.fits')

condition = np.logical_or.reduce((t_orig['CLASS'] == 'bll',
                                  t_orig['CLASS'] == 'BLL'))
index = np.where(condition)
t_orig_bll = t_orig[index]

condition = np.logical_or.reduce((t_orig['CLASS'] == 'fsrq',
                                  t_orig['CLASS'] == 'FSRQ'))
index = np.where(condition)
t_orig_fsrq = t_orig[index]


# Select detected sources
t = t_orig[t_orig['Sigma'] >= 5.]

# Select bll and fsrqs
condition = np.logical_or.reduce((t['CLASS'] == 'bll',
                                  t['CLASS'] == 'BLL'))
index = np.where(condition)
t_bll = t[index]

condition = np.logical_or.reduce((t['CLASS'] == 'fsrq',
                                  t['CLASS'] == 'FSRQ'))
index = np.where(condition)
t_fsrq = t[index]

# Start plot
plt.figure()
ax = plt.gca()

# Significance VS redshift (with and without cut-off)
ax.errorbar(t_bll['Redshift'], t_bll['Sigma'],
            xerr=0, yerr=t_bll['SigmaErr'],
            label='BL Lacs', c='blue', fmt='o')

ax.errorbar(t_fsrq['Redshift'], t_fsrq['Sigma'],
            xerr=0, yerr=t_fsrq['SigmaErr'],
            label='FSRQs', c='red', fmt='^')

ax.errorbar(t_bll['Redshift'], t_bll['SigmaCut'],
            xerr=0, yerr=t_bll['SigmaCutErr'],
            label='BL Lacs with cut-off', c='blue', fmt='o', fillstyle='none')

ax.errorbar(t_fsrq['Redshift'], t_fsrq['SigmaCut'],
            xerr=0, yerr=t_fsrq['SigmaCutErr'],
            label='FSRQs with cut-off', c='red', fmt='^', fillstyle='none')

# Draw segment line to join same sources
lines_list = list()
for row in t:
    x1 = row['Redshift']
    y1 = row['Sigma']    
    x2 = row['Redshift']
    y2 = row['SigmaCut']    
    lines_list.append([(x1, y1), (x2, y2)])
lc = mc.LineCollection(lines_list, colors='b', linewidths=1, alpha=0.2)
ax.add_collection(lc)

ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlim(0.02, 2.)

ax.set_xlabel('Redshift')
ax.set_ylabel('Significance (20 h livetime)')

ax.plot([0.02,2.],[5.,5.], c='green', lw=1, ls='--')

ax.text(0.35, 0.93, 'TECHNICAL PLOT', style='italic', fontweight='bold',
        bbox={'facecolor':'white', 'alpha':1, 'pad':3},
        verticalalignment='bottom', horizontalalignment='right',
        transform=ax.transAxes, fontsize=15, color='black')

plt.grid(True, which='both')
ax.legend(loc='best', numpoints=1)
plt.tight_layout()
plt.show()

out_file = './plots/'
if not os.path.exists(out_file):
    os.makedirs(out_file)

out_file += 'cut_off_agn_pop.pdf'

plt.savefig(out_file, format='pdf')
