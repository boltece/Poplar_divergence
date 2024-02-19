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
    Model with split between pop 1 and 2, ancestral gene flow (mA) occurs.
    Population 3 is admixed from pop 1 and 2.
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after origin.
    m12: Migrations from pop 2 to pop 1
    m32: Migration rate from pop 2 to pop3 (2*Na*m)
    m31: Migration rate from pop 1 to pop3
    T1: The scaled time between the split of pops 1 vs 2 and origin of 3
    (in units of 2*Na generations).
    T2: The scaled time interval associated with isolation (no gene flow)
    T3: The scaled time between the origin of pop 3  and the present
    (in units of 2*Na generations) with asymmetrical gene flow
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




#### run 1
upper_bound = [10, 10, 10, 3, 3, 3, 10, 10, 10, 10, 10, 10, 10, 0.90]
lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.10]
p0 = array([1, 1, 1, 0.1, 0.1, 0.1, 1, 1, 1, 1, 1, 1, 1, 0.50])




p0 = dadi.Misc.perturb_params(p0, upper_bound=upper_bound,lower_bound=lower_bound)
print('starting set: {0}'.format(p0))
func_ex = dadi.Numerics.make_extrap_func(func)
popt = dadi.Inference.optimize_log(p0, data, func_ex, pts_1,lower_bound=lower_bound,upper_bound=upper_bound,verbose=len(p0), maxiter=20)



print('Best-fit parameters: {0}'.format(popt))



model = func_ex(popt, ns, pts_1)
ll_model = dadi.Inference.ll_multinom(model, data)
print('Maximum log composite likelihood: {0}'.format(ll_model))



theta = dadi.Inference.optimal_sfs_scaling(model, data)
print('Optimal value of theta: {0}'.format(theta))



plt=dadi.Plotting.plot_3d_comp_multinom(model, data, vmin=1, resid_range=3, pop_ids =("coastal_tricho","dbalsam", "interior_tricho"))

matplotlib.pyplot.savefig("Admix_origin_T1_GF_T2_noGF_T3_asymGF_unfolded_run1.pdf")


