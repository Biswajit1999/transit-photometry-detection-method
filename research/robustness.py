"""Synthetic transit-survey diagnostics; no observational detections."""
import numpy as np


def ar1_noise(n, sigma, rho, seed=0):
    """Stationary Gaussian AR(1) noise with marginal standard deviation sigma."""
    if n < 1 or sigma < 0 or not np.isfinite(sigma) or not -1 < rho < 1:
        raise ValueError('invalid AR(1) parameters')
    rng = np.random.default_rng(seed)
    x = np.empty(n)
    x[0] = rng.normal(scale=sigma)
    for i in range(1, n):
        x[i] = rho*x[i-1] + rng.normal(scale=sigma*np.sqrt(1-rho*rho))
    return x


def depth_gls(flux, in_transit, sigma, rho=0):
    """Fit baseline and box depth with exact stationary AR(1) precision.

    Assumes regularly spaced samples; use separate segments for gaps.
    """
    y = np.asarray(flux, float)
    mask = np.asarray(in_transit, bool)
    if y.ndim != 1 or mask.shape != y.shape or len(y) < 3 or not mask.any() or mask.all():
        raise ValueError('need in- and out-of-transit samples')
    if not np.isfinite(y).all() or sigma <= 0 or not -1 < rho < 1:
        raise ValueError('invalid data or covariance')
    X = np.column_stack([np.ones(len(y)), -mask.astype(float)])
    W = np.empty_like(X)
    z = np.empty_like(y)
    W[0], z[0] = X[0]/sigma, y[0]/sigma
    scale = sigma*np.sqrt(1-rho*rho)
    W[1:] = (X[1:]-rho*X[:-1])/scale
    z[1:] = (y[1:]-rho*y[:-1])/scale
    info = W.T @ W
    covariance = np.linalg.inv(info)
    beta = covariance @ (W.T @ z)
    return float(beta[1]), float(np.sqrt(covariance[1,1]))


def window_recovery_probability(period, duration, baseline, cadence, gaps=(), epochs=500):
    """Fraction of uniformly sampled transit epochs with >=3 observed events.

    Geometric window completeness only: no noise or detection threshold.
    """
    if min(period,duration,baseline,cadence) <= 0 or duration >= period or epochs < 1:
        raise ValueError('invalid survey parameters')
    t = np.arange(0,baseline,cadence)
    keep = np.ones(len(t), bool)
    for a,b in gaps:
        if not 0 <= a < b <= baseline:
            raise ValueError('invalid gap')
        keep &= ~((t>=a)&(t<b))
    t=t[keep]
    fractions=[]
    for epoch in np.linspace(0,period,epochs,endpoint=False):
        centers=np.arange(epoch,baseline+period,period)
        counts=sum(np.any(np.abs(t-c) <= duration/2) for c in centers if c < baseline)
        fractions.append(counts>=3)
    return float(np.mean(fractions))
