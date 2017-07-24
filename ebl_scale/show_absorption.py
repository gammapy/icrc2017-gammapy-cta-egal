import numpy as np

import matplotlib.pyplot as plt

import astropy.units as u

from gammapy.spectrum.models import Absorption
from gammapy.utils.energy import EnergyBounds

import os

def cut_off(energy):
    """Value of the attenuation due to the cut-off"""
    return np.exp(-energy / (1 * u.TeV / (1 + z)))

# EBL model
absorption = Absorption.read_builtin('dominguez')

# Intervals for the attenuation computation
zmin = 0.001
zmax = 0.5
zstep= 0.001

emin = 0.1 * u.TeV
emax = 30 * u.TeV
venergy = EnergyBounds.equal_log_spacing(emin=emin,
                                         emax=emax,
                                         nbins=1000)

# Initialise outputs
n = int(round((zmax - zmin) / zstep, 0))
z_ebl = np.zeros(n)
e_ebl = np.zeros(n) * u.TeV

z_cutoff = np.zeros(n)
e_cutoff = np.zeros(n)  * u.TeV

# Initialise attenuation (here, 63 %)
sigma1 = np.exp(-1)

z = zmin
count = 0
while (z <= zmax):

    # Absorption from EBL
    tau_ebl = absorption.evaluate(venergy, parameter=z)
    idx = np.abs(tau_ebl - sigma1).argmin()
    energy_ebl = venergy[idx]
    e_ebl[count] = energy_ebl
    z_ebl[count] = z

    # Absorption from cut-off
    tau_cutoff = cut_off(venergy)
    idx = np.abs(tau_cutoff - sigma1).argmin()
    energy_cutoff = venergy[idx]
    e_cutoff[count] = energy_cutoff
    z_cutoff[count] = z

    z = z + zstep
    count = count + 1

# Plot results
plt.figure()
ax = plt.gca()
ax.plot(z_ebl, e_ebl.value, label=r'EBL - 63% absorption', color='green')
ax.plot(z_cutoff, e_cutoff.value, label=r'Cut-off - 63% attenuation', color='red', ls='--')

ax.set_yscale("log", nonposx='clip')
ax.set_ylim([emin.value, emax.value])
ax.set_xlim([zmin, zmax])
ax.set_ylabel('Energy [{}]'.format(venergy.unit))
ax.set_xlabel('Redshift')

plt.grid(which='both')
plt.legend(loc='best')
plt.tight_layout()
plt.show()


# Save figure in pdf
out_file = './plots/'
if not os.path.exists(out_file):
    os.makedirs(out_file)

out_file += 'attenuation.pdf'

plt.savefig(out_file, format='pdf')
