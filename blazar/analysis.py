"""
Example script to analyse 3FHL data. Does:
 - source selection
 - source simulation with and without cut-off in Fermi spectra
 - save results (significance and so on)
"""
from __future__ import print_function
import numpy as np
import astropy.units as u
from astropy.table import Table, Column, vstack
from gammapy.spectrum.models import AbsorbedSpectralModel, Absorption, SpectralModel
from gammapy.utils.modeling import Parameter, ParameterList
from gammapy.scripts import CTAPerf
from gammapy.scripts.cta_utils import (CTAObservationSimulation, Target,
                                       ObservationParameters)
from gammapy.catalog import SourceCatalog3FHL


class AbsorbedSpectralModelExpoCutOff(SpectralModel):
    """
    Class to handle any spectra with a cut-off
    """

    def __init__(self, spectral_model, cut_off):
        self.spectral_model = spectral_model
        self.cut_off = cut_off

        param_list = []
        for param in spectral_model.parameters.parameters:
            param_list.append(param)

        # Add parameter to the list
        par = Parameter('cut_off', cut_off,
                        parmin=10 * u.MeV, parmax=100 * u.TeV,
                        frozen=True)
        param_list.append(par)

        self.parameters = ParameterList(param_list)

    def evaluate(self, energy, **kwargs):
        """Evaluate the model at a given energy."""
        flux = self.spectral_model(energy=energy)
        absorption = np.exp(- energy / (self.cut_off.to(energy.unit)))
        return flux * absorption


def average_simu(target, perf, params, n=20):
    """Compute average significance for `n` trials."""
    sigmas = np.empty(n)
    for i in range(n):
        simu = CTAObservationSimulation.simulate_obs(
            perf=perf, target=target, obs_param=params, random_state=i,
        )
        sigmas[i] = simu.total_stats_safe_range.sigma

    return sigmas.mean(), sigmas.std()


fermi = SourceCatalog3FHL()
table_3fhl = fermi.table
print('Fermi sources:', len(table_3fhl))

# Select blazars
condition = np.logical_or.reduce((table_3fhl['CLASS'] == 'bll    ',
                                  table_3fhl['CLASS'] == 'BLL    ',
                                  table_3fhl['CLASS'] == 'fsrq   ',
                                  table_3fhl['CLASS'] == 'FSRQ   '))
index = np.where(condition)
table = table_3fhl[index]
print('Fermi blazars:', len(table))

# Select sources with redshift
condition = np.isfinite(table['Redshift']) == True
index = np.where(condition)
table = table[index]
print('Fermi blazars with z:', len(table))

# Select sources according to zenith angle
south_site_lat = -25.
north_site_lat = 18.
zen_lim = 30.

# South
south_lat_max = south_site_lat + zen_lim
south_lat_min = south_site_lat - zen_lim
condition = np.logical_and(table['DEJ2000'] >= south_lat_min,
                           table['DEJ2000'] <= south_lat_max)
index = np.where(condition)
table_south = table[index]
print('South blazars: {} (bll={}, fsrq={})'.format(
    len(table_south),
    len(table_south[table_south['CLASS'] == 'bll    ']),
    len(table_south[table_south['CLASS'] == 'fsrq   '])))

# North
lat_max = north_site_lat + zen_lim
lat_min = north_site_lat - zen_lim
condition = np.logical_and(table['DEJ2000'] >= lat_min,
                           table['DEJ2000'] <= lat_max)
index = np.where(condition)
table_north = table[index]
print('North blazars: {} (bll={}, fsrq={})'.format(
    len(table_north),
    len(table_north[table_north['CLASS'] == 'bll    ']),
    len(table_north[table_north['CLASS'] == 'fsrq   '])))

# Concatenate table
table = vstack([table_south, table_north])
print('Selected blazars: {} (bll={}, fsrq={})'.format(
    len(table),
    len(table[table['CLASS'] == 'bll    ']),
    len(table[table['CLASS'] == 'fsrq   '])))

# Replace table, dirty but efficient
fermi.table = table

# Performance
irf_dir = '$GAMMAPY_EXTRA/datasets/cta/perf_prod2/point_like_non_smoothed/'
south_irf_file = 'South_5h.fits.gz'
north_irf_file = 'North_5h.fits.gz'

cta_perf_south = CTAPerf.read(irf_dir + south_irf_file)
cta_perf_north = CTAPerf.read(irf_dir + north_irf_file)

# EBL
absorption = Absorption.read_builtin('dominguez')

# Cut-off
cut_off = 1. * u.TeV

# Variables to save
n_src = len(table)
src_sigma = np.zeros(n_src, dtype=float)
src_sigma_err = np.zeros(n_src, dtype=float)
src_sigma_cut = np.zeros(n_src, dtype=float)
src_sigma_cut_err = np.zeros(n_src, dtype=float)
src_class = np.zeros(n_src, dtype=np.dtype('a4'))
src_name = np.zeros(n_src, dtype=np.dtype('a18'))
src_ana = np.zeros(n_src, dtype=np.dtype('a1'))
src_z = np.zeros(n_src, dtype=float)

# Loop on sources
for idx, src in enumerate(fermi):

    name = src.data['Source_Name']
    redshift = src.data['Redshift']
    fermi_model = src.spectral_model
    dec = src.data['DEJ2000'].value
    is_south = False
    ana = 'N'
    src_type = (src.data['CLASS']).strip()
    if south_lat_min <= dec < south_lat_max:
        is_south = True
        ana = 'S'
    print('Processing {} {} ({}/{})'.format(src_type, name, idx + 1, len(fermi.table)))

    # Observation parameters
    emin = 0.05 * u.TeV
    emax = 100. * u.TeV
    alpha = 0.2 * u.Unit('')
    livetime = 20. * u.h

    obs_param = ObservationParameters(
        alpha=alpha, livetime=livetime,
        emin=emin, emax=emax,
    )

    # Performance
    cta_perf = cta_perf_south if is_south else cta_perf_north

    # Absorbed spectral model
    abs_model = AbsorbedSpectralModel(spectral_model=fermi_model,
                                      absorption=absorption,
                                      parameter=redshift)

    # Absorbed spectral model with cut-off
    cut_off_abs_model = AbsorbedSpectralModelExpoCutOff(
        spectral_model=abs_model,
        cut_off=cut_off / (1 + redshift),
    )

    # Simulations
    target = Target(name=name, model=abs_model)
    av_sigma, rms_sigma = average_simu(target, cta_perf, obs_param, 20)

    target_cut = Target(name=name, model=cut_off_abs_model)
    av_sigma_cut, rms_sigma_cut = average_simu(target_cut, cta_perf, obs_param, 20)

    print('==> sigma={:.2f}, sigma_err={:.2f}, sigma_cut={:.2f}, sigma_cut_err={}'.format(
        av_sigma,
        rms_sigma,
        av_sigma_cut,
        rms_sigma_cut)
    )

    src_sigma[idx] = av_sigma
    src_sigma_err[idx] = rms_sigma
    src_sigma_cut[idx] = av_sigma_cut
    src_sigma_cut_err[idx] = rms_sigma_cut
    src_class[idx] = src_type.strip()
    src_name[idx] = name
    src_ana[idx] = ana
    src_z[idx] = redshift

t = Table()
t['Source_Name'] = Column(src_name, unit='', description='')
t['CLASS'] = Column(src_class, unit='', description='')
t['Redshift'] = Column(src_z, unit='', description='')
t['Site'] = Column(src_ana, unit='', description='')
t['Sigma'] = Column(src_sigma, unit='', description='')
t['SigmaErr'] = Column(src_sigma_err, unit='', description='')
t['SigmaCut'] = Column(src_sigma_cut, unit='', description='')
t['SigmaCutErr'] = Column(src_sigma_cut_err, unit='', description='')

# Save results
filename = 'results_cut{:.2f}{}.fits'.format(cut_off.value, cut_off.unit)
t.write(filename, format='fits', overwrite=True)
