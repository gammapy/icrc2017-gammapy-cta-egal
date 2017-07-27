from __future__ import print_function
import json
import numpy as np
import astropy.units as u
from gammapy.spectrum.models import (PowerLaw, PowerLaw2,
                                     AbsorbedSpectralModel, Absorption)
from gammapy.scripts import CTAPerf
from gammapy.scripts.cta_utils import (CTAObservationSimulation, Target,
                                       ObservationParameters)

import gammapy
print(gammapy.__version__)
print(gammapy)

# Define start and stop times for the simulation
t_start = np.logspace(np.log10(20), np.log10(172800), num=20) * u.s
t_stop = t_start  # could choose a different time binning here

# GRB spectral template and other parameters
grb = dict(
    name='GRB 080916C',
    redshift=3.,
    index_template=2.0,
    int_flux=500e-5 * u.Unit('1 / (cm2 s)'),
    int_emin=0.1 * u.GeV,
    int_emax=10 * u.GeV,
    index_decay=1.7,
    tp=6.5 * u.s
)

# Add starting observation time, 20 s after peak
grb['t0'] = 20 * u.s + grb['tp']

# Add differential flux at 1 GeV for the simulation
grb['eref'] = 1 * u.GeV
grb['amplitude_template'] = PowerLaw2.evaluate(
    energy=grb['eref'],
    index=grb['index_template'],
    amplitude=grb['int_flux'],
    emin=grb['int_emin'],
    emax=grb['int_emax']
)
# Add differential flux at t0
grb['amplitude_template_t0'] = grb['amplitude_template'] * np.power(
    grb['t0'] / grb['tp'], -grb['index_decay'])

# EBL absorption
absorption = Absorption.read_builtin('dominguez')

# Perf
irf_dir = '$GAMMAPY_EXTRA/datasets/cta/perf_prod2/point_like_non_smoothed/'
irf_file = irf_dir + 'North_0.5h.fits.gz'
perf = CTAPerf.read(irf_file)


def compute_sigma(start, stop):
    """Compute significance for a single start and stop time."""
    delta_t = stop - start  # livetime

    if delta_t <= 0:  # only interested in physical intervals
        return 0

    # Obs params
    livetime = delta_t
    emin = 0.03 * u.TeV
    emax = 1. * u.TeV
    alpha = 0.2 * u.Unit('')
    obs_param = ObservationParameters(alpha=alpha,
                                      livetime=livetime,
                                      emin=emin,
                                      emax=emax)

    # Computation of the amplitude at the average time interval
    # average time ** decay = 1/(1-decay) * int t**-decay dt / int dt * t0**decay
    index_decay = grb['index_decay']
    amplitude_template_t0 = grb['amplitude_template_t0']
    t0 = grb['t0']

    time_from_t0 = np.power(stop, 1. - index_decay) - np.power(start, 1. - index_decay)
    time_from_t0 /= (stop - start)
    time_from_t0 /= 1 - index_decay
    amplitude_t = amplitude_template_t0 * time_from_t0 * np.power(t0, index_decay)

    # Model
    pwl = PowerLaw(index=grb['index_template'],
                   amplitude=amplitude_t,
                   reference=grb['eref'])

    # EBL absorbed model
    spectral_model = AbsorbedSpectralModel(spectral_model=pwl,
                                           absorption=absorption,
                                           parameter=grb['redshift'])
    # Target
    target = Target(name='GRB 080916C', model=spectral_model)

    # Simulation
    simu = CTAObservationSimulation.simulate_obs(
        perf=perf, target=target, obs_param=obs_param, random_state=0,
    )
    stats = simu.total_stats_safe_range
    return stats.sigma


# Compute significance for all start and stop times
sigma = np.zeros((len(t_start), len(t_stop)))
for i, start in enumerate(t_start):
    for j, stop in enumerate(t_stop):
        sigma[i][j] = compute_sigma(start, stop)

# Store all results in a Python dictionary, then a JSON file
results = {}
grb_info = dict(name=grb['name'], redshift=grb['redshift'])
results['grb_info'] = grb_info
results['sigma'] = sigma.tolist()
results['t_start'] = t_start.value.tolist()
results['t_stop'] = t_stop.value.tolist()

filename = 'grb_data.json'
print('Writing', filename)
json.dump(results, open(filename, 'w'), indent=4)
