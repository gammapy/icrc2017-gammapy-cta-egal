import numpy as np
import astropy.units as u

from gammapy.scripts.cta_utils import (CTAObservationSimulation,
                                       ObservationParameters, Target)
from gammapy.spectrum.models import SpectralModel, Absorption
from gammapy.utils.modeling import Parameter, ParameterList

def get_differential_flux_from_integral_flux(emin, emax,
                                             eref,
                                             gamma,
                                             int_flux):
    """
    Compute differential flux from integrated flux
    Parameters
    ----------
    emin : `~astropy.units.Quantity`
        Minimal energy for integrated flux
    emax : `~astropy.units.Quantity`
        Maximal energy for integrated flux
    eref : `~astropy.units.Quantity`
        Reference energy for the differential flux
    gamma : `~astropy.units.Quantity`
        Spectral index
    int_flux : `~astropy.units.Quantity`
        Integrated flux beteween emin and emax
    """
    diff_flux = int_flux
    diff_flux *= np.power(eref, -gamma)
    diff_flux /= (np.power(emin, 1.-gamma) - np.power(emax, 1.-gamma))
    diff_flux *= (gamma-1.)
    return diff_flux


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
                        parmin=10 *  u.MeV, parmax=100 * u.TeV,
                        frozen=True)
        param_list.append(par)

        self.parameters = ParameterList(param_list)

    def evaluate(self, energy, **kwargs):
        """Evaluate the model at a given energy."""
        flux = self.spectral_model(energy=energy)
        absorption = np.exp(- energy/(self.cut_off.to(energy.unit)))
        return flux * absorption

def average_simu(target, perf, params, n=20):
    """
    Function to compute average significance for N trials
    """
    
    sigma_list = []
    for i in range(n):
        simu = CTAObservationSimulation.simulate_obs(perf=perf,
                                                     target=target,
                                                     obs_param=params)

        stats = simu.total_stats_safe_range
        sigma = stats.sigma
        sigma_list.append(sigma)

    av_sigma = 0
    for sigma in sigma_list:
        av_sigma += sigma
    av_sigma = av_sigma / float(len(sigma_list))

    sigma_rms = 0
    for sigma in sigma_list:
        sigma_rms += (sigma-av_sigma)**2
    sigma_rms = sigma_rms / float(len(sigma_list))

    return av_sigma, np.sqrt(sigma_rms)


def simulate(model, perf, livetime, emin, emax, obs_id):
    """
    Simulate a source for a set of parameters
    """
    # observation parameters
    alpha = 0.2 * u.Unit('')
    livetime = livetime
    obs_param = ObservationParameters(alpha=alpha,
                                      livetime=livetime,
                                      emin=emin,
                                      emax=emax)

    # target (not really useful)
    target = Target(name=str(obs_id), model=model)

    # simulation
    simu = CTAObservationSimulation.simulate_obs(perf=perf,
                                                 target=target,
                                                 obs_param=obs_param,
                                                 obs_id=obs_id)
    return simu


from sherpa.models import ArithmeticModel
from sherpa.models import Parameter as SherpaParameter
from gammapy.utils.scripts import make_path
from astropy.table import Table
from gammapy.utils.nddata import NDDataArray, BinnedDataAxis
from gammapy.utils.energy import EnergyBounds

class SherpaAbsorption(ArithmeticModel):
    """
    Sherpa model for EBL.

    Read table in gammapy-extra and create an interpolator for exp[-alpha * tau(E,Z)]
    http://pysherpa.blogspot.kr/2010/06/user-defined-sherpa-model-types-using.html
    """

    def __init__(self,
                 name='EBL',
                 filename='dominguez'):

        # Create EBL data array the same way it is done in Gammapy
        absorption = Absorption.read_builtin(filename)

        self.ebl = absorption.data

        # define redshift aprameter
        par_min = absorption.data.axes[0].lo[0]
        par_max = absorption.data.axes[0].hi[-1]
        self.redshift = SherpaParameter(name, 'z', par_min.value,
                                        frozen=True,
                                        min=par_min.value,
                                        max=par_max.value)

        # define parameter alpha
        self.alpha = SherpaParameter(name, 'alpha', 1., frozen=True, min=0, max=3)

        ArithmeticModel.__init__(self, name, (self.redshift, self.alpha))

    def calc(self, p, xlo, xhi=None, *args, **kwargs):
        """
        Params
        `p`   list of ordered parameter values.
        `x`   ndarray of domain values, bin midpoints or lower
              bin edges.
        `xhi` ndarray of upper bin edges.

        returns ndarray of calculated function values.
        """

        z = p[0]
        alpha = p[1]

        energies = xlo * u.keV   # Cause we're using sherpa

        if xhi is not None:
            energies += xhi * u.keV   # Cause we're using sherpa
            energies *= 0.5

        return np.power(self.ebl.evaluate(parameter=z, energy=energies), alpha)
