import numpy as np
import qmeq

import qmeq
# from qmeq.approach.base.pauli import PauliSolver as QmeqPauliSolver
# from qmeq.builder_base import StateIndexingDM
# from qmeq.specfunc import Func1

np.set_printoptions(linewidth=256, suppress=True)

from pauli_solver_lib import PauliSolverLib

def run_solvers(bRunQmeq=True, bRunCpp=True):
    eps1 = -5.0
    eps2 = 0.0
    eps3 = 5.0
    t = 0.1
    W = 1.0
    VBias = 0.0
    GammaS = 0.01
    GammaT = 0.01
    coeffE = 1.0
    coeffT = 1.0
    Temp = 0.5
    muS = 0.0
    muT = 0.0
    
    NSingle = 3
    NStates = 2**NSingle
    NLeads = 2
    
    qmeq_res = None
    if bRunQmeq:
        print("\nRunning QMeq solver...")
        
        indexing = StateIndexingDM(NSingle)
        pauli = QmeqPauliSolver(indexing)
        
        pauli.add_parameters({
            'e1': eps1, 'e2': eps2, 'e3': eps3,
            't12': t, 't23': t*coeffT, 't13': 0.0,
            'gammaL1': GammaS, 'gammaR1': 0.0,
            'gammaL2': 0.0, 'gammaR2': 0.0,
            'gammaL3': 0.0, 'gammaR3': GammaT,
            'temp': Temp,
            'muL': muS, 'muR': muT + VBias,
            'dband': W
        })
        
        pauli.solve()
        
        qmeq_res = {
            'current': pauli.current[1],
            'energies': pauli.energies,
            'probabilities': pauli.si.get_vec_dm(pauli.si.dvec),
            'kernel': pauli.kern
        }
        print("\nQMeq probabilities:", qmeq_res['probabilities'])
        print("QMeq current:", qmeq_res['current'])
    
    cpp_res = None
    if bRunCpp:
        print("\nRunning C++ solver...")
        
        print( "\n\n######################################################################" )
        print( "######################################################################" )
        print( "######################################################################" )
        print( "######################################################################" )
        print( "######################################################################" )
        print( "### Running C++ solver /home/prokop/git_SW/qmeq/cpp/pauli_solver.hpp " )
        print("\nDEBUG: System parameters:")
        print(f"NSingle={NSingle}, NLeads={NLeads}")
        print(f"eps1={eps1}, eps2={eps2}, eps3={eps3}")
        print(f"t={t}, W={W}, VBias={VBias}")
        print(f"GammaS={GammaS}, GammaT={GammaT}")
        print(f"coeffE={coeffE}, coeffT={coeffT}")
        
        solver = PauliSolverLib(verbosity=1)
        cpp_probs, cpp_kernel, cpp_rhs, cpp_factors = solver.solve(
            nsingle=NSingle,
            eps1=eps1, eps2=eps2, eps3=eps3,
            t=t, W=W, VBias=VBias,
            GammaS=GammaS, GammaT=GammaT,
            coeffE=coeffE, coeffT=coeffT,
            muS=muS, muT=muT,
            Temp=Temp
        )
        
        cpp_res = {
            'current': cpp_rhs[1],
            'energies': cpp_kernel[0],  # First row contains energies
            'probabilities': cpp_probs,
            'kernel': cpp_kernel,
        }
        print("\nC++ states:", [bin(i)[2:].zfill(NSingle) for i in range(NStates)])
        print("C++ energies:", cpp_res['energies'])
    
    return qmeq_res, cpp_res

def main():
    qmeq_res, cpp_res = run_solvers()
    
    if qmeq_res is not None and cpp_res is not None:
        print("\nComparing results:")
        print("Energy difference:", np.max(np.abs(qmeq_res['energies'] - cpp_res['energies'])))
        print("Probability difference:", np.max(np.abs(qmeq_res['probabilities'] - cpp_res['probabilities'])))
        print("Kernel difference:", np.max(np.abs(qmeq_res['kernel'] - cpp_res['kernel'])))
        print("Current difference:", np.abs(qmeq_res['current'] - cpp_res['current']))
#!/usr/bin/env python3

import numpy as np
from sys import path
path.insert(0, '/home/prokop/bin/home/prokop/venvs/ML/lib/python3.12/site-packages/qmeq/')
import qmeq

# setup numpy print options to infinite line length
np.set_printoptions(linewidth=256, suppress=True)
#np.set_printoptions(linewidth=256, suppress=False)

from pauli_solver_lib import PauliSolverLib

def run_solvers( bRunQmeq=True, bRunCpp=True ):

    # System parameters
    NSingle = 3  # number of impurity states
    NLeads  = 2   # number of leads

    # Parameters (in meV)
    eps1 = eps2 = eps3 = -10.0
    t     = 0.0      # direct hopping
    W     = 20.0     # inter-site coupling
    VBias = 0.0  # bias voltage

    # Lead parameters
    muS    = 0.0    # substrate chemical potential
    muT    = 0.0    # tip chemical potential
    Temp   = 0.224 # temperature in meV
    DBand  = 1000.0 # lead bandwidth
    GammaS = 0.20 # coupling to substrate
    GammaT = 0.05 # coupling to tip

    # Tunneling amplitudes
    VS = np.sqrt(GammaS/np.pi)  # substrate
    VT = np.sqrt(GammaT/np.pi)  # tip

    # Position-dependent coefficients
    coeffE = 0.4
    coeffT = 0.3
    
    # One-particle Hamiltonian
    hsingle = {(0,0): eps1-coeffE*VBias, (0,1): t, (0,2): t,
               (1,1): eps2, (1,2): t,
               (2,2): eps3}
    
    # Two-particle Hamiltonian: inter-site coupling
    coulomb = {(0,1,1,0): W,
               (1,2,2,1): W,
               (0,2,2,0): W}
    
    # Leads: substrate (S) and scanning tip (T)
    mu_L   = {0: muS, 1: muT + VBias}
    Temp_L = {0: Temp, 1: Temp}
    
    # Coupling between leads (1st number) and impurities (2nd number)
    TLeads = {(0,0): VS,         # S <-- 1
              (0,1): VS,         # S <-- 2
              (0,2): VS,         # S <-- 3
              (1,0): VT,         # T <-- 1
              (1,1): coeffT*VT,  # T <-- 2
              (1,2): coeffT*VT}  # T <-- 3


    verbosity = 1
    
    # Run QmeQ solver
    if bRunQmeq:
        print( "\n\n######################################################################" )
        print( "######################################################################" )
        print( "######################################################################" )
        print( "######################################################################" )
        print( "######################################################################" )
        print( "\n### Running QmeQ Pauli solver /home/prokop/git_SW/qmeq/qmeq/approach/base/pauli.py " )
        system = qmeq.Builder(NSingle, hsingle, coulomb, NLeads, TLeads, mu_L, Temp_L, DBand,   kerntype='Pauli', indexing='Lin', itype=0, symq=True,   solmethod='lsqr', mfreeq=0)
        system.appr.verbosity = verbosity  # Set verbosity after instance creation
        system.verbosity = verbosity
        system.solve()
        qmeq_res = {
            'current': system.current[1],
            'energies': system.Ea,
            'probabilities': system.phi0,
            'kernel': system.kern,
        }
        print("\nQmeQ hsingle:", hsingle)
        print("QmeQ coulomb:", coulomb)
        print("QmeQ states:", [bin(i)[2:].zfill(NSingle) for i in range(2**NSingle)])
        print("QmeQ energies:", system.Ea)

    # Run C++ solver
    if bRunCpp:
        print( "\n\n######################################################################" )
        print( "######################################################################" )
        print( "######################################################################" )
        print( "######################################################################" )
        print( "######################################################################" )
        print( "### Running C++ solver /home/prokop/git_SW/qmeq/cpp/pauli_solver.hpp " )
        print("\nDEBUG: System parameters:")
        print(f"NSingle={NSingle}, NLeads={NLeads}")
        print(f"eps1={eps1}, eps2={eps2}, eps3={eps3}")
        print(f"t={t}, W={W}, VBias={VBias}")
        print(f"GammaS={GammaS}, GammaT={GammaT}")
        print(f"VS={VS}, VT={VT}")
        print(f"coeffE={coeffE}, coeffT={coeffT}")
        
        solver = PauliSolverLib(verbosity=1)
        cpp_probs, cpp_kernel, cpp_rhs, cpp_factors = solver.solve(
            nsingle=NSingle,
            eps1=eps1, eps2=eps2, eps3=eps3,
            t=t, W=W, VBias=VBias,
            GammaS=GammaS, GammaT=GammaT,
            coeffE=coeffE, coeffT=coeffT,
            muS=muS, muT=muT,
            Temp=Temp
        )
        
        cpp_res = {
            'current': cpp_rhs[1],
            'energies': cpp_kernel[0],  # First row contains energies
            'probabilities': cpp_probs,
            'kernel': cpp_kernel,
        }
        print("\nC++ states:", [bin(i)[2:].zfill(NSingle) for i in range(2**NSingle)])
        print("C++ energies:", cpp_res['energies'])
    
    return qmeq_res, cpp_res

def compare_results(qmeq_res, cpp_res, tol=1e-8):
    """Compare results from both solvers"""
    print("\n\n#### Comparing QmeQ vs C++ results:")

    diff = np.max(np.abs(qmeq_res['energies'] - cpp_res['energies']))
    if diff > tol:
        print("Energies:   diff:", diff)
        print("QmeQ:", qmeq_res['energies'])
        print("C++: ", cpp_res['energies'])
    else:
        print(f"Energies:   OK (diff({diff}) < tol({tol}))")
    
    diff = np.max(np.abs(qmeq_res['probabilities'] - cpp_res['probabilities']))
    if diff > tol:
        print("Probabilities:   diff:", diff)
        print("QmeQ:", qmeq_res['probabilities'])
        print("C++: ", cpp_res['probabilities'])
    else:
        print(f"Probabilities:   OK (diff({diff}) < tol({tol}))")
    
    diff = np.max(np.abs(qmeq_res['kernel'] - cpp_res['kernel']))
    if diff > tol:
        print("Kernel:   diff:", diff)
        print("QmeQ:", qmeq_res['kernel'])
        print("C++: ", cpp_res['kernel'])
    else:
        print(f"Kernel:   OK (diff({diff}) < tol({tol}))")
    
    diff = np.max(np.abs(qmeq_res['current'] - cpp_res['current']))
    if diff > tol:
        print("Current:   diff:", diff)
        print("QmeQ:", qmeq_res['current'])
        print("C++: ", cpp_res['current'])
        print("Relative diff:", abs(qmeq_res['current'] - cpp_res['current'])/abs(qmeq_res['current']))
    else:
        print(f"Current:   OK (diff({diff}) < tol({tol}))")

def main():
    qmeq_res, cpp_res = run_solvers()
    compare_results(qmeq_res, cpp_res)

if __name__ == "__main__":
    main()
