import numpy as np
from scipy.optimize import curve_fit


def fit_biGaussian(x,y, bounds=None, p0=None):
    """
    Derive the best fit curve for the flux-time distribution
    """

    if bounds is None:
        bounds = ([-np.inf,-1e-5,-1e-5,0,-1e5],[np.inf,np.inf,np.inf,np.inf,np.inf])

    if p0 is None:
        p0 = (x[np.argmax(y)],np.std(x)/3,np.std(x),1, 0)
    popt, pcov = curve_fit(biGaussian,x,y,p0=p0,bounds=bounds)
    return popt, pcov


def biGaussian(x,x0,sig1,sig2,A,const):
    # combine 2 gaussian:
    return A*np.exp(-0.5*((x-x0)/
        (sig1*(x<x0)+sig2*(x>=x0)))**2)+const