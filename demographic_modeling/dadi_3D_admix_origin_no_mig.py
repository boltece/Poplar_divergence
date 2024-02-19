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

os.chdir('/storage/work/ceb6313/LDpruned/')
data = dadi.Spectrum.from_file("Poplar_3pops_diff_topo_20chromosomes_per_pop_LDpruned_unfolded.sfs")


# In[3]:


ns=data.sample_sizes



# In[4]:


pts_1 = [40,50,60]




def ADMIX(params, ns, pts):
    """
    Model with split between pop 1 and 2, gene flow (m) does not occur.
    Population 3 is admixed from pop 1 and 2.
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after origin.
    T1: The scaled time between the split of pops 1 vs 2
    T2: The scaled time between the origin of pop 3 and the present
    (in units of 2*Na generations).
    f: Fraction of pop 3 derived from pop 1
    (with fraction 1-f derived from pop 2).
    """
    # 8 parameters
    nu1, nu2, nu3, T1, T2, f = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, T1, nu1=nu1, nu2=nu2, m12=0, m21=0)
    phi = PhiManip.phi_2D_to_3D_admix(phi, f, xx, xx, xx)
    phi = Integration.three_pops(
        phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=0, m21=0, m23=0, m32=0, m13=0, m31=0)

    fs = Spectrum.from_phi(phi, ns, (xx, xx, xx))
    return fs
    
func = ADMIX




#### changed upper and lower bounds to reflect first run being closer to lower bound 
upper_bound = [10, 10, 10, 10, 10, 0.90]
lower_bound = [0.001, 0.001, 0.001, 0.001, 0.001, 0.10]
p0 = array([1, 1, 1, 1, 1, 0.50])




p0 = dadi.Misc.perturb_params(p0, upper_bound=upper_bound,lower_bound=lower_bound)
print('starting set: {0}'.format(p0))
func_ex = dadi.Numerics.make_extrap_log_func(func)
popt = dadi.Inference.optimize_log(p0, data, func_ex, pts_1,lower_bound=lower_bound,upper_bound=upper_bound,verbose=len(p0), maxiter=20)



print('Best-fit parameters: {0}'.format(popt))



model = func_ex(popt, ns, pts_1)
ll_model = dadi.Inference.ll_multinom(model, data)
print('Maximum log composite likelihood: {0}'.format(ll_model))



theta = dadi.Inference.optimal_sfs_scaling(model, data)
print('Optimal value of theta: {0}'.format(theta))



plt=dadi.Plotting.plot_3d_comp_multinom(model, data, vmin=1, resid_range=3, pop_ids =("coastal_tricho","dbalsam", "interior_tricho"))

matplotlib.pyplot.savefig("SI_admix_origin_unfolded_run1.pdf")


