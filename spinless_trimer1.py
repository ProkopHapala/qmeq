## Simulation of a molecular trimer using master equation approach
## as implemented in the QmeQ package (github.com/gedaskir/qmeq)
## Vladislav Pokorny; 2024; pokornyv@fzu.cz

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline as IUS

from time import time,ctime
from sys import path,argv,exit
#path.insert(0, '/home/pokorny/bin/usr/local/lib64/python3.12/site-packages/')
path.insert(0, '/home/prokop/bin/home/prokop/venvs/ML/lib/python3.12/site-packages/qmeq/')
#path.insert(0, '/home/pokorny/bin/qmeq-1.1/')
import qmeq

print('# '+argv[0]+' start, '+ctime())

NSingle = 3 ## number of impurity states
NStates = 2**NSingle
NLeads  = 2 ## number of leads

## all parmeters are in meV
## quantum impurities
eps1 = eps2 = eps3 = -10.0
t    = 0.0   ## direct hopping
U    = 220.0 ## on-site Coulomb interaction, useless for spinless case
W    = 20.0  ## inter-site coupling

## leads
#VBias = float(argv[1])
muS   = 0.0    ## substrate chemical potential
muT   = 0.0    ## scanning tip chemical potential
Temp  = 0.224  ## (2.6K) temperature in meV, 1meV = 11.6K, 1K = 0.0862meV
DBand = 1000.0 ## lead bandwidth

GammaS = 0.20 ## coupling to substrate
GammaT = 0.05 ## coupling to scanning tip
## tunneling amplitudes weighted by lead DOS
VS = np.sqrt(GammaS/np.pi)
VT = np.sqrt(GammaT/np.pi)



VBiasMin = 0.0
VBiasMax = 60.0
dVBias   = 1.0





VBias_A  = np.arange(VBiasMin,VBiasMax,dVBias)
NPoints  = len(VBias_A)
I_A      = np.zeros(NPoints)
Eigen_A  = np.zeros([NPoints,NStates])

print(f'# eps1: {eps1} eps2: {eps2} eps3: {eps3} t: {t} U: {U} W: {W}')
print(f'# GammaS: {GammaS} GammaT: {GammaT} VS: {VS:.4f} VT: {VT:.4f}')
print(f'# muS: {muS:.4f} muT: {muT:.4f} Temp: {Temp:.4f} DBand: {DBand:.1f}')
print(f'# VBiasMin: {VBiasMin:.1f} VBiasMax: {VBiasMax:.1f} dVBias: {dVBias:.3f} NPoints: {NPoints}')

## coefficients that simulate indirectly the position of the tip
coeffE = 0.4
coeffT = 0.3

for k in range(NPoints):
	VBias = VBias_A[k]
	## one-particle Hamiltonian
	H1p = {(0,0): eps1-coeffE*VBias, (0,1): t, (0,2): t,
     	  (1,1): eps2, (1,2): t,
     	  (2,2): eps3}

	## two-particle Hamiltonian: inter-site coupling
	H2p = {(0,1,1,0): W,
	       (1,2,2,1): W,
	       (0,2,2,0): W}

	## leads: substrate (S) and scanning tip (T)
	mu_L   = {0: muS,  1: muT + VBias}
	Temp_L = {0: Temp, 1: Temp}

	## coupling between leads (1st number) and impurities (2nd number)
	TLeads = {(0,0): VS, # S <-- 1
	          (0,1): VS, # S <-- 2
	          (0,2): VS, # S <-- 3
	          (1,0): VT, # T <-- 1
	          (1,1): coeffT*VT, # T <-- 2
	          (1,2): coeffT*VT} # T <-- 3

	#Builder(nsingle=0, hsingle={}, coulomb={}, nleads=0, tleads={}, mulst={}, tlst={}, dband={}, 
	#indexing=None, kpnt=None, kerntype='Pauli', symq=True, norm_row=0, solmethod=None, itype=0, 
	#dqawc_limit=10000, mfreeq=False, phi0_init=None, mtype_qd=<class 'complex'>, mtype_leads=<class 'complex'>, 
	#symmetry=None, herm_hs=True, herm_c=False, m_less_n=True)

	kerntype = 'Pauli'
	#kerntype = '1vN'
	#system = qmeq.Builder(NSingle, H1p, H2p, NLeads, TLeads, mu_L, Temp_L, DBand, kerntype=kerntype)
	system = qmeq.Builder(NSingle, H1p, H2p, NLeads, TLeads, mu_L, Temp_L, DBand, kerntype=kerntype, indexing='Lin', itype=0, symq=True, solmethod='lsqr', mfreeq=0)

	system.solve()
	#print(system.phi0)
	#system = qmeq.Builder(NSingle, H1p, H2p, NLeads, TLeads, mu_L, Temp_L, DBand, kerntype='2vN', kpnt=2**8)
	#system.solve(niter=7)
	I_A[k] = system.current[1]
	Eigen_A[k,:] = system.Ea

## calculate dI/dV from current
I      = IUS(VBias_A,I_A)
dIdV   = I.derivative()
dIdV_A = dIdV(VBias_A)

plt.plot(VBias_A,I_A)
plt.plot(VBias_A,dIdV_A)
plt.show()





# for k in range(NPoints):
# 	print('{0: .5f} {1: 3e} {2: 3e}'.format(VBias_A[k],I_A[k],dIdV_A[k]))
# print('\n\n')

# for k in range(NPoints):
# 	print('{0: .5f}'.format(VBias_A[k]),*np.around(Eigen_A[k],6),sep='\t')

# if 0:
# 	print('Reduced density matrix:')
# 	## for Pauli ME, DM is diagonal, we get 2**NSingle values
# 	## for 1vN the first 2**NSingle values are the diagonal elements, rest are real and imag part of
# 	## the off-diagonal elements within given charge blocks, e.g., for NSingle=3 we get 8+6+6=20 elements
# 	print(system.phi0)
# 	print('Many−body tunneling amplitudes:')
# 	print(system.Tba)

## spinless_trimer1.py END ##

