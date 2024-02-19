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
data = dadi.Spectrum.from_file("Poplar_3pops_20chromosomes_per_pop_LDpruned_unfolded.sfs")


# In[3]:


ns=data.sample_sizes


# # recurring gene flow model (GF) with two categories of loci experiencing heterogeneous migration rates throughout the genome

# In[4]:


pts_1 = [40,50,60]

def SI(params, ns, pts):
	nuB,nuA,nuC,nuI,T1,T2 = params
	"""
	3-populations, strict isolation.

	nuB:  Current size of balsam population, after split.
	nuC:  Current size of coastal tricho population.
	nuI: Current size of interior tricho population
	nuA: Ancestral Ne of Coastal and Interior
	T1:   Time for divergence between balsam and tricho species with no gene flow
	T2: Time interval with gene flow and divergence of coastal and interior
	ns:   Size of fs to generate.
	pts:  Number of points to use in grid for evaluation.
	"""
	xx = dadi.Numerics.default_grid(pts)
	### Calculate the neutral spectrum
	phi = dadi.PhiManip.phi_1D(xx)
	phi = dadi.PhiManip.phi_1D_to_2D(xx,phi)
	phi = dadi.Integration.two_pops(phi,xx,T1,nu1=nuB, nu2=nuA, m12=0, m21=0)
	phi = dadi.PhiManip.phi_2D_to_3D_split_2(xx,phi)
	phi = dadi.Integration.three_pops(phi,xx,T2,nu1=nuB, nu2=nuC, nu3=nuI,m12=0, m21=0,m23=0, m32=0, m13=0, m31=0)
	fs = dadi.Spectrum.from_phi(phi,ns,(xx,xx,xx))
	return fs

func = SI





upper_bound = [10, 10, 10, 10, 10, 10]
lower_bound = [1e-2, 1e-2, 1e-2, 1e-2, 1e-2, 1e-2]
p0 = array([1, 1, 1, 1, 1, 1])




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



plt=dadi.Plotting.plot_3d_comp_multinom(model, data, vmin=1, resid_range=3, pop_ids =("balsam","coastal_tricho", "interior_tricho"))

matplotlib.pyplot.savefig("LDpruned_SI_3pop_unfolded_maxiter20_run1.pdf")

