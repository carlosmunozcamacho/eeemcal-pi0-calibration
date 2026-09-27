# SPDX-License-Identifier: LGPL-3.0-or-later
# Adapted from calibration/pi0/src/peak_diagnostics.py and its earlier pilot.find_peak parent.
"""Small robustness extension of pilot.find_peak's binned Poisson model.

Gaussian bin integrals plus an exponential or positive linear background.
Parameter scales are normalized for the increased counts; the peak and width
remain floating. Significance is local, with no scan/trials correction.
"""
import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.special import ndtr


def background_integral(lo, hi, slope, form):
    lo, hi = np.asarray(lo), np.asarray(hi)
    if form == 'linear':
        return (hi-lo) + slope*((hi-.14)**2-(lo-.14)**2)/.2
    # Stable even around slope=0, where subtracting two exponentials divided
    # by slope corrupts numerical optimizer gradients for a steep continuum.
    x = slope*(hi-lo)/.2
    with np.errstate(divide='ignore', invalid='ignore'):
        sinhc = np.where(np.abs(x) < 1e-4, 1+x*x/6+x**4/120, np.sinh(x)/x)
    return (hi-lo)*np.exp(slope*((hi+lo)/2-.14)/.1)*sinhc


def window_yields(fit, lo, hi):
    mu, sigma = fit['peak_GeV'], fit['sigma_GeV']
    signal = fit['signal_yield']*(ndtr((hi-mu)/sigma)-ndtr((lo-mu)/sigma))
    background = fit['background_yield_in_fit_range'] * float(background_integral(lo, hi, fit['background_slope'], fit['background_form']) /
                                                               background_integral(*fit.get('fit_range_GeV', [.04,.24]), fit['background_slope'], fit['background_form']))
    return dict(window_GeV=[float(lo), float(hi)], signal=float(signal), background=float(background),
                purity=float(signal/(signal+background)) if signal+background > 0 else 0.)


def fit_peak(mass, bins=50, background_form='exponential', fixed_shape=None, fit_range=(.04,.24), criteria=None):
    criteria = criteria or dict(entries=60, signal=15, significance=4., mean_bounds=[.076,.194], sigma_bounds=[.0041,.0279])
    if background_form not in ['exponential', 'linear']:
        raise ValueError('Unsupported background form')
    edges = np.linspace(*fit_range, bins+1)
    y, _ = np.histogram(mass, edges)
    n = max(float(y.sum()), 1.)
    slope_bounds = (-.98, .98) if background_form == 'linear' else (-4., 4.)
    def background(q):
        integrals = background_integral(edges[:-1], edges[1:], q[1], background_form)
        return n*np.exp(q[0])*integrals/integrals.sum()
    def expected(q):
        return background(q) + n*q[2]*np.diff(ndtr((edges-.1*q[3])/(.01*q[4])))
    def loss(mu):
        mu = np.maximum(mu, 1e-12)
        return float(np.sum(mu-y*np.log(mu)))
    # With normalized background shape the null total has the exact MLE sum(y).
    # Only its slope needs a bounded scalar minimization; do not rely on a
    # potentially stalled multidimensional null fit for significance.
    bfit = minimize_scalar(lambda slope: loss(background([0., slope])*float(y.sum())/n),
                           bounds=slope_bounds, method='bounded', options={'xatol': 1e-10})
    bounds = [(-16., 2.), slope_bounds, (0., 2.), (.75, 1.95), (.4, 2.8)]
    if fixed_shape:
        bounds[3] = (fixed_shape[0]/.1, fixed_shape[0]/.1)
        bounds[4] = (fixed_shape[1]/.01, fixed_shape[1]/.01)
    candidates = []
    for center in ([fixed_shape[0]] if fixed_shape else [.09, .12, .15, .18]):
        q0 = [np.log(.8), float(bfit.x), .2, center/.1, fixed_shape[1]/.01 if fixed_shape else 1.2]
        candidates.append(minimize(lambda q: loss(expected(q)), q0, bounds=bounds,
                                   method='L-BFGS-B', options={'maxiter': 600, 'ftol': 1e-12, 'gtol': 1e-7}))
    best = min(candidates, key=lambda r: r.fun)
    q = best.x
    significance = float(np.sqrt(max(0., 2*(bfit.fun-best.fun))))
    mu, sigma, signal = float(.1*q[3]), float(.01*q[4]), float(n*q[2])
    interior = criteria['mean_bounds'][0] < mu < criteria['mean_bounds'][1] and criteria['sigma_bounds'][0] < sigma < criteria['sigma_bounds'][1]
    passed = bool(best.success and y.sum() >= criteria['entries'] and signal >= criteria['signal'] and significance >= criteria['significance'] and interior)
    result = dict(fit_range_GeV=list(fit_range), identification_criteria=criteria, identified=passed, fit_converged=bool(best.success), entries_in_fit_range=int(y.sum()),
                  signal_yield=signal, peak_GeV=mu, sigma_GeV=sigma,
                  local_likelihood_significance=significance, background_form=background_form,
                  background_slope=float(q[1]), background_yield_in_fit_range=float(n*np.exp(q[0])),
                  bin_width_GeV=float(edges[1]-edges[0]), hist_edges=edges.tolist(), hist_counts=y.tolist(),
                  fit_expected=expected(q).tolist(), background_expected=background(q).tolist(),
                  fixed_shape=fixed_shape, uncertainties=None,
                  null_fit_converged=bool(bfit.success), null_background_slope=float(bfit.x),
                  null_nll=float(bfit.fun), signal_nll=float(best.fun),
                  criteria='Original >=60 entries, >=15 fitted signal, local significance >=4 and interior mean/width; no trials correction')
    # Local expected-information curvature, not a gain pseudoinverse error.
    # Report only when all five parameters have an invertible information matrix.
    if interior and signal > 0 and not fixed_shape:
        steps = [1e-4]*5
        derivatives = []
        for k, step in enumerate(steps):
            qp, qm = q.copy(), q.copy(); qp[k] += step; qm[k] -= step
            derivatives.append((expected(qp)-expected(qm))/(2*step))
        jac = np.array(derivatives).T
        information = jac.T @ (jac/np.maximum(expected(q), 1e-12)[:, None])
        try:
            if np.linalg.cond(information) < 1e10:
                covariance = np.linalg.inv(information)
                if np.all(np.diag(covariance) > 0):
                    errors = np.sqrt(np.diag(covariance))
                    result['uncertainties'] = dict(peak_GeV=float(.1*errors[3]), sigma_GeV=float(.01*errors[4]),
                                                   signal_yield=float(n*errors[2]),
                                                   method='local expected-information standard errors; no model/scan systematic')
        except np.linalg.LinAlgError:
            pass
    result['windows'] = {str(width)+'sigma': window_yields(result, mu-width*sigma, mu+width*sigma) for width in [1., 1.5, 2.]}
    result['windows']['fraction030'] = window_yields(result, mu*.7, mu*1.3)
    return result


def robustness(mass):
    return {name: fit_peak(mass, bins=bins, background_form=form)
            for name, bins, form in [('5MeV_exponential', 40, 'exponential'),
                                     ('4MeV_exponential', 50, 'exponential'),
                                     ('8MeV_exponential', 25, 'exponential'),
                                     ('5MeV_linear', 40, 'linear')]}
