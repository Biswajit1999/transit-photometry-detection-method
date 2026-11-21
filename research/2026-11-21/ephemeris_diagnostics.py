"""Conditional linear ephemeris diagnostics. Input times and uncertainties share units."""
import numpy as np

def weighted_linear_ephemeris(cycle, midpoint, sigma):
    """Fit T(E)=T0+P*E with independent Gaussian timing uncertainties."""
    e,t,s=(np.asarray(x,dtype=float) for x in (cycle,midpoint,sigma))
    if e.ndim!=1 or t.shape!=e.shape or s.shape!=e.shape or e.size<3 or not np.all(np.isfinite(e+t+s)) or np.any(s<=0) or np.unique(e).size<2:
        raise ValueError("invalid timing measurements")
    center=np.average(e,weights=1/s**2)
    x=e-center
    A=np.column_stack((np.ones(e.size),x))/s[:,None]
    beta,_,rank,_=np.linalg.lstsq(A,t/s,rcond=None)
    if rank!=2: raise ValueError("rank-deficient ephemeris")
    cov=np.linalg.inv(A.T@A)
    t0=beta[0]-beta[1]*center
    transform=np.array([[1,-center],[0,1]])
    return np.array([t0,beta[1]]),transform@cov@transform.T,t-(t0+beta[1]*e)

def forecast_timing_uncertainty(cycle,covariance):
    """Propagate fitted T0,P covariance to future transit epochs."""
    e=np.asarray(cycle,dtype=float)
    cov=np.asarray(covariance,dtype=float)
    if cov.shape!=(2,2) or not np.all(np.isfinite(cov)) or not np.allclose(cov,cov.T) or np.linalg.eigvalsh(cov).min() < -1e-12 or not np.all(np.isfinite(e)):
        raise ValueError("invalid covariance or cycle")
    var=cov[0,0]+2*e*cov[0,1]+e**2*cov[1,1]
    return np.sqrt(np.maximum(var,0))
