#!/usr/bin/env python
# coding: utf-8

# In[1]:

import os
import matplotlib
matplotlib.use('PDF')
import dadi
import pylab
import matplotlib.pyplot as plt
import numpy as np
from numpy import array
from dadi import Misc,Spectrum,Numerics,PhiManip,Integration


# In[2]:



os.chdir('/storage/work/ceb6313/LDpruned')
data= dadi.Spectrum.from_file("Poplar_3pops_diff_topo_20chromosomes_per_pop_LDpruned_unfolded.sfs")




ns =data.sample_sizes



pts_1 = [40,50,60]




def ADMIX(params, ns, pts):
    """
    Model with split between pop 1 and 2, gene flow does occur followed by a period of isolation
    Population 3 is admixed from pop 1 and 2 at T3 and asymGF follows to present.
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after origin.
    mA: ancestral migration (symmetrical)
    m12: Migrations from pop 2 to pop 1 (m21 is vice versa)
    m32: Migration rate from pop 2 to pop3 (2*Na*m)
    m31: Migration rate from pop 1 to pop3
    T1: The scaled time between the split of pops 1 vs 2 with gene flow and origin of 3
    (in units of 2*Na generations).
    T2: The scaled time between the split of pops 1 vs 2 withOUT gene flow and origin of 3
    T3: The scaled time between the origin of pop 3 and the present with asym gene flow
    (in units of 2*Na generations).
    f: Fraction of pop 3 derived from pop 1
    (with fraction 1-f derived from pop 2).
    """
    # 8 parameters
    nu1, nu2, nu3, T1, T2, T3, mA, mCB, mBC, mBI, mIB, mCI, mIC, f = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, T1, nu1=nu1, nu2=nu2, m12=mA, m21=mA)
    phi = Integration.two_pops(phi, xx, T2, nu1=nu1, nu2=nu2, m12=0, m21=0)
    phi = PhiManip.phi_2D_to_3D_admix(phi, f, xx, xx, xx)
    phi = Integration.three_pops(phi, xx, T3, nu1=nu1, nu2=nu2, nu3=nu3, m12=mCB, m21=mBC, m23=mBI, m32=mIB, m13=mCI, m31=mIC)
    fs = Spectrum.from_phi(phi, ns, (xx, xx, xx))
    return fs
    
func = ADMIX





from dadi import Inference
from dadi.Spectrum_mod import Spectrum

import numpy

from dadi import Godambe

cache = {}



two_pt_deriv_test = False

def hessian_elem(func, f0, p0, ii, jj, eps, args=(), one_sided=None):
    """
    Calculate element [ii][jj] of the Hessian matrix, a matrix
    of partial second derivatives w.r.t. to parameters ii and jj
        
    func: Model function
    f0: Evaluation of func at p0
    p0: Parameters for func
    eps: List of absolute step sizes to use for each parameter when taking
         finite differences.
    args: Additional arguments to func
    one_sided: Optionally, pass in a sequence of length p0 that determines
               whether a one-sided derivative will be used for each parameter.
    """
    # Note that we need to specify dtype=float, to avoid this being an integer
    # array which will silently fail when adding fractional eps.
    if one_sided is None:
        one_sided = [False]*len(p0)

    pwork = numpy.array(p0, copy=True, dtype=float)
    if ii == jj:
        if pwork[ii] != 0 and not one_sided[ii]:
            pwork[ii] = p0[ii] + eps[ii]
            fp = func(pwork, *args)
            
            pwork[ii] = p0[ii] - eps[ii]
            fm = func(pwork, *args)
            
            element = (fp - 2*f0 + fm)/eps[ii]**2
        else:
            pwork[ii] = p0[ii] + 2*eps[ii]
            fpp = func(pwork, *args)
            
            pwork[ii] = p0[ii] + eps[ii]
            fp = func(pwork, *args)

            element = (fpp - 2*fp + f0)/eps[ii]**2
    else:
        if pwork[ii] != 0 and pwork[jj] != 0 and not one_sided[ii] and not one_sided[jj]:
            # f(xi + hi, xj + h)
            pwork[ii] = p0[ii] + eps[ii]
            pwork[jj] = p0[jj] + eps[jj]
            fpp = func(pwork, *args)
            
            # f(xi + hi, xj - hj)
            pwork[ii] = p0[ii] + eps[ii]
            pwork[jj] = p0[jj] - eps[jj]
            fpm = func(pwork, *args)
            
            # f(xi - hi, xj + hj)
            pwork[ii] = p0[ii] - eps[ii]
            pwork[jj] = p0[jj] + eps[jj]
            fmp = func(pwork, *args)
            
            # f(xi - hi, xj - hj)
            pwork[ii] = p0[ii] - eps[ii]
            pwork[jj] = p0[jj] - eps[jj]
            fmm = func(pwork, *args)

            element = (fpp - fpm - fmp + fmm)/(4 * eps[ii]*eps[jj])
        else:
            # f(xi + hi, xj + h)
            pwork[ii] = p0[ii] + eps[ii]
            pwork[jj] = p0[jj] + eps[jj]
            fpp = func(pwork, *args)
            
            # f(xi + hi, xj)
            pwork[ii] = p0[ii] + eps[ii]
            pwork[jj] = p0[jj]
            fpm = func(pwork, *args)
            
            # f(xi, xj + hj)
            pwork[ii] = p0[ii]
            pwork[jj] = p0[jj] + eps[jj]
            fmp = func(pwork, *args)
            
            element = (fpp - fpm - fmp + f0)/(eps[ii]*eps[jj])
    return element





def get_godambe(func_ex, grid_pts, all_boot, p0, data, eps, log=False,
                just_hess=False, boot_theta_adjusts=[]):
    """
    Godambe information and Hessian matrices

    func_ex: Model function
    grid_pts: Number of grid points to evaluate the model function
    all_boot: List of bootstrap frequency spectra
    p0: Best-fit parameters for func_ex.
    data: Original data frequency spectrum
    eps: Fractional stepsize to use when taking finite-difference derivatives
         Note that if eps*param is < 1e-6, then the step size for that parameter
         will simply be eps, to avoid numerical issues with small parameter
         perturbations.
    log: If True, calculate derivatives in terms of log-parameters
    just_hess: If True, only evaluate and return the Hessian matrix
    boot_theta_adjusts: Factors by which to adjust theta for each bootstrap
                        sample, relative to full data theta.
    """
    ns = data.sample_sizes
    if not boot_theta_adjusts:
        boot_theta_adjusts = numpy.ones(len(all_boot))

    # Cache evaluations of the frequency spectrum inside our hessian/J 
    # evaluation function
    def func(params, data, theta_adjust=1):
        key = (func_ex.__hash__(), tuple(params), tuple(ns), tuple(grid_pts))
        if key not in cache:
            cache[key] = func_ex(params, ns, grid_pts)
        # theta_adjust deals with bootstraps that need  different thetas
        fs = theta_adjust*cache[key]
        return Inference.ll(fs, data)
    def log_func(logparams, data, theta_adjust=1):
        return func(numpy.exp(logparams), data, theta_adjust)

    # First calculate the observed hessian.
    # theta_adjust defaults to 1.
    if not log:
        hess = -get_hess(func, p0, eps, args=[data])
    else:
        hess = -get_hess(log_func, numpy.log(p0), eps, args=[data])

    if just_hess:
        return hess





def FIM_uncert(func_ex, grid_pts, p0, data, log=False, multinom=True, eps=0.01):
    """
    Parameter uncertainties from Fisher Information Matrix

    Returns standard deviations of parameter values.

    func_ex: Model function
    all_boot: List of bootstrap frequency spectra
    p0: Best-fit parameters for func_ex
    data: Original data frequency spectrum
    eps: Fractional stepsize to use when taking finite-difference derivatives.
         Note that if eps*param is < 1e-6, then the step size for that parameter
         will simply be eps, to avoid numerical issues with small parameter
         perturbations.
    log: If True, assume log-normal distribution of parameters. Returned values
         are then the standard deviations of the *logs* of the parameter values,
         which can be interpreted as relative parameter uncertainties.
    multinom: If True, assume model is defined without an explicit parameter for
              theta. Because uncertainty in theta must be accounted for to get
              correct uncertainties for other parameters, this function will
              automatically consider theta if multinom=True. In that case, the
              final entry of the returned uncertainties will correspond to
              theta.
    """
    if multinom:
        func_multi = func_ex
        model = func_multi(p0, data.sample_sizes, grid_pts)
        theta_opt = Inference.optimal_sfs_scaling(model, data)
        p0 = list(p0) + [theta_opt]
        func_ex = lambda p, ns, pts: p[-1]*func_multi(p[:-1], ns, pts)
    H = get_godambe(func_ex, grid_pts, [], p0, data, eps, log, just_hess=True)
    return numpy.sqrt(numpy.diag(numpy.linalg.inv(H)))





def get_hess(func, p0, eps, args=()):
    """
    Calculate Hessian matrix of partial second derivatives. 
    Hij = dfunc/(dp_i dp_j)
    
    func: Model function
    p0: Parameter values to take derivative around
    eps: Fractional stepsize to use when taking finite-difference derivatives
         Note that if eps*param is < 1e-6, then the step size for that parameter
         will simply be eps, to avoid numerical issues with small parameter
         perturbations.
    args: Additional arguments to func
    """
    # Calculate step sizes for finite-differences.
    eps_in = eps
    eps = numpy.empty([len(p0)])
    one_sided = [False]*len(p0)
    for i, pval in enumerate(p0):
        if pval != 0:
            # Account for floating point arithmetic issues
            print('type pval=', type(pval))
            print('type eps_in =', type(eps_in))
            if pval*eps_in < 1e-6:
                eps[i] = eps_in
                one_sided[i] = True
            else:
                eps[i] = eps_in*pval
        else:
            # Account for parameters equal to zero
            eps[i] = eps_in

    f0 = func(p0, *args)
    hess = numpy.empty((len(p0), len(p0)))
    for ii in range(len(p0)):
        for jj in range(ii, len(p0)):
            hess[ii][jj] = hessian_elem(func, f0, p0, ii, jj, eps, args=args, one_sided=one_sided)
            hess[jj][ii] = hess[ii][jj]
    return hess






func_ex = dadi.Numerics.make_extrap_func(func)



#### parameter estimates from the best supported (AIC) model of divergence history
popt_fold= array([2.60129866, 7.12370083, 3.83476941, 2.9930843,  0.22958452, 1.81570141,  5.3387714,  1.02041246, 0.36376898, 0.25480264, 0.87861879, 1.42509928,  1.28423296, 0.34808687])



Uncert_FIM = FIM_uncert(func_ex, pts_1, popt_fold, data,log=False, multinom=True, eps=0.001)



print('Estimated parameter standard deviations from FIM_eps0.001: {0}'.format(Uncert_FIM))



### changing eps= to get a confidence interval for some params
Uncert_FIM2 = FIM_uncert(func_ex, pts_1, popt_fold, data,log=False, multinom=True, eps=0.00001)



print('Estimated parameter standard deviations from FIM_eps0.00001: {0}'.format(Uncert_FIM2))



