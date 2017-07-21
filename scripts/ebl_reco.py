import os
import sys

import numpy as np

import astropy.units as u
from gammapy.spectrum.models import AbsorbedSpectralModel, PowerLaw, Absorption
from gammapy.scripts import CTAPerf

# Utilties functions
from utils import simulate

# Output directory
outdir = './data/'

# CTA performance (point-like responses from prod2)
irf_dir = '$GAMMAPY_EXTRA/datasets/cta/perf_prod2/point_like_non_smoothed/'
irf_file = 'South_5h.fits.gz'
perf = CTAPerf.read(irf_dir + irf_file)

# EBL model
absorption = Absorption.read_builtin('dominguez')

# Simulate 2 ~random sources, PWL + EBL absorption
# the "simu" object belongs to the gamma.spectrum.SpectrumObservation class

# 1st source
pwl = PowerLaw(amplitude=1.e-12 * u.Unit('1 / (cm2 s TeV)'),
               index=2.5,
               reference=1 * u.TeV)
pwl_absorbed = AbsorbedSpectralModel(
    spectral_model=pwl,
    absorption=absorption,
    parameter=0.2  # redshift
)
simu = simulate(
    model=pwl_absorbed,
    perf=perf,
    livetime=100 * u.h,
    emin=0.03 * u.TeV,
    emax=10. * u.TeV,
    obs_id=0)
simu.write(outdir=outdir, use_sherpa=True)

# 2nd source (change redshift and obs_id)
pwl_absorbed = AbsorbedSpectralModel(
    spectral_model=pwl,
    absorption=absorption,
    parameter=0.3  # redshift
)
simu = simulate(
    model=pwl_absorbed,
    perf=perf,
    livetime=100 * u.h,
    emin=0.03 * u.TeV,
    emax=10. * u.TeV,
    obs_id=1)
simu.write(outdir=outdir, use_sherpa=True)

# Sherpa imports
# http://cxc.cfa.harvard.edu/contrib/sherpa/
import sherpa.astro.datastack as sh
from utils import SherpaAbsorption

# Add sherpa model
sh.add_model(SherpaAbsorption)

# Load data, set model and initialise parameters
obs_id = 0
sh.load_data(id=obs_id, filename=outdir + '/pha_obs' + str(obs_id) + '.fits')
sherpa_model = 'powlaw1d.pwl_' + str(obs_id) + '*sherpaabsorption.ebl_' + str(obs_id)
sh.set_source(obs_id, sherpa_model)
pwl_0.gamma = 2.
pwl_0.ref = (1 * u.TeV).to('keV')  # sherpa works in keV
pwl_0.ampl = (1.e-12 * u.Unit('1 / (cm2 s TeV)')).to('1 / (cm2 s keV)')   # sherpa works in keV
ebl_0.redshift = 0.2
ebl_0.alpha = 1.  # EBL scale

obs_id = 1
sh.load_data(id=obs_id, filename=outdir + '/pha_obs' + str(obs_id) + '.fits')
sherpa_model = 'powlaw1d.pwl_' + str(obs_id) + '*sherpaabsorption.ebl_' + str(obs_id)
sh.set_source(obs_id, sherpa_model)
pwl_1.gamma = 2.
pwl_1.ref = (1 * u.TeV).to('keV')  # sherpa works in keV
pwl_1.ampl = (1.e-12 * u.Unit('1 / (cm2 s TeV)')).to('1 / (cm2 s keV)')   # sherpa works in keV
ebl_1.redshift = 0.3
ebl_1.alpha = 1.  # EBL scale

# Free alpha and link the parameters for the two models
ebl_0.alpha.thaw()
ebl_1.alpha.link = ebl_0.alpha

# Print
print(sh.get_model(0))
print(sh.get_model(1))

# Set threshold by hand, normally it should be already
# defined when calling the write method ==> TO BE CHECK !!!!
emin = 0.03 * u.TeV
emax = 10 * u.TeV
sh.notice((emin.to('keV')).value,(emax.to('keV')).value)

# Fit with Cash statistics (full forward folding)
# http://cxc.harvard.edu/sherpa/statistics/#cash
sh.set_stat('WStat')
sh.fit()

# Compute confidence intervals
sh.conf()

# Access to fitted parameters
res = sh.get_conf_results()
pvals1 = zip(res.parvals, res.parmins, res.parmaxes)
pvals2 = [(v, l, h) for (v,l,h) in pvals1]
dres = dict(zip(res.parnames, pvals2))

alpha = dres['ebl_0.alpha'][0]
alpha_errmin = dres['ebl_0.alpha'][1]
alpha_errmax = dres['ebl_0.alpha'][2]

print('EBL scale={}, Err-={}, Err+={}'.format(
    alpha,
    alpha+alpha_errmin,
    alpha+alpha_errmax)
)

# Clean memory
sh.clean()
