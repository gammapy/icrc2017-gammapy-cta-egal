import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from gammapy.spectrum.models import (PowerLaw, PowerLaw2,
                                     AbsorbedSpectralModel, Absorption)
from gammapy.scripts import CTAPerf
from gammapy.scripts.cta_utils import (CTAObservationSimulation, Target,
                                       ObservationParameters)

# Time bins definition
n_bins = 20
tmin = 20
tmax = 172800

xmin, xmax = np.log10([tmin, tmax])
t_start = np.logspace(xmin, xmax, n_bins) * u.s

xmin, xmax = np.log10([tmin, tmax])
t_stop = np.logspace(xmin, xmax, n_bins) * u.s

# GRB template
grb_template = dict(name='GRB 080916C',
                    redshift=3.,
                    index_template=2.0,
                    int_flux=500e-5 * u.Unit('1 / (cm2 s)'),
                    int_emin=0.1 * u.GeV,
                    int_emax=10 * u.GeV,
                    index_decay=1.7,
                    tp=6.5 * u.s)

# Add starting observation time, 20 s after peak
grb_template['t0'] = 20 * u.s + grb_template['tp']

# Add differential flux at 1 GeV for the simulation
grb_template['eref'] = 1 * u.GeV
grb_template['amplitude_template'] = PowerLaw2.evaluate(
    energy=grb_template['eref'],
    index=grb_template['index_template'],
    amplitude=grb_template['int_flux'],
    emin=grb_template['int_emin'],
    emax=grb_template['int_emax']
)
# Add differential flux at t0
grb_template['amplitude_template_t0'] = grb_template['amplitude_template'] * np.power(
    grb_template['t0'] / grb_template['tp'], -grb_template['index_decay'])

# EBL absorption
absorption = Absorption.read_builtin('dominguez')

# Perf
irf_dir = '$GAMMAPY_EXTRA/datasets/cta/perf_prod2/point_like_non_smoothed/'
irf_file = irf_dir + 'North_0.5h.fits.gz'
perf = CTAPerf.read(irf_file)

# Significance results
sigma = np.zeros((len(t_start), len(t_stop)))

# Loop 
for i, start in enumerate(t_start):

    for j, stop in enumerate(t_stop):

        delta_t = stop - start  # livetime

        if delta_t <= 0:  # only interested in physical intervals
            sigma[i][j] = 0.
            continue

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
        index_decay = grb_template['index_decay']
        amplitude_template_t0 = grb_template['amplitude_template_t0']
        t0 = grb_template['t0']

        time_from_t0 = np.power(stop, 1. - index_decay) - np.power(start, 1. - index_decay)
        time_from_t0 /= (stop - start)
        time_from_t0 /= 1 - index_decay
        amplitude_t = amplitude_template_t0 * time_from_t0 * np.power(t0, index_decay)

        # Model
        pwl = PowerLaw(index=grb_template['index_template'],
                       amplitude=amplitude_t,
                       reference=grb_template['eref'])

        # EBL absorbed model 
        spectral_model = AbsorbedSpectralModel(spectral_model=pwl,
                                               absorption=absorption,
                                               parameter=grb_template['redshift'])
        # Target
        target = Target(name='GRB 080916C', model=spectral_model)

        # Simulation
        simu = CTAObservationSimulation.simulate_obs(perf=perf,
                                                     target=target,
                                                     obs_param=obs_param)
        stats = simu.total_stats_safe_range
        sigma[i][j] = stats.sigma

# Get ride of non extreme significance values (inf or negative for esthetic)
mask = np.isfinite(sigma)
sigma[np.where(mask == False)] = 0.

mask = sigma < 0.
sigma[np.where(mask == True)] = 0.

# Change the orientation of the axes
sigma.transpose()

# Plot the significance for each time interval
plt.figure()
x, y = np.meshgrid(t_start.value, t_stop.value)
plt.pcolormesh(x, y, sigma, cmap='hot')
plt.ylim(t_start.value[-1], t_start.value[0])

plt.xlabel('tstop (s)')
plt.ylabel('tstart (s)')
plt.xscale('log')
plt.yscale('log')
plt.colorbar(label='Detection significance')
plt.title('{} (z={})'.format(grb_template['name'], grb_template['redshift']))
plt.tight_layout()

filename = 'plots/grb_twindow.png'
print('Writing {}'.format(filename))
plt.savefig(filename, format='png', dpi=300)
