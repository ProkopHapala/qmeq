"""Module containing python functions, which generate first order Pauli kernel."""

import numpy as np
import itertools

from ...wrappers.mytypes import doublenp

from ...specfunc.specfunc import func_pauli
from ..aprclass import Approach

from ...config import debug_print, verb_print, verb_print_

# ---------------------------------------------------------------------------------------------------
# Pauli master equation
# ---------------------------------------------------------------------------------------------------
class ApproachPauli(Approach):

    kerntype = 'pyPauli'

    def __init__(self, qd):
        """
        Initialization of the Pauli approach.

        Parameters
        ----------
        qd : QuantumDot
            Quantum dot system.
        """
        Approach.__init__(self, qd)
        self.paulifct = np.zeros((self.si.nleads, self.si.ndm1, 2), dtype=float)
        self.verbosity = 4  # Add verbosity control

    def get_kern_size(self):
        return self.si.npauli

    def prepare_arrays(self):
        Approach.prepare_arrays(self)
        nleads, ndm1 = self.si.nleads, self.si.ndm1
        self.paulifct = np.zeros((nleads, ndm1, 2), dtype=doublenp)

    def clean_arrays(self):
        Approach.clean_arrays(self)
        self.paulifct.fill(0.0)

    def generate_fct(self):
        """
        Make factors used for generating Pauli master equation kernel.
        """
        print("ApproachPauli.generate_fct(self):")
        #exit()

        if self.verbosity > 3:
            print("\nQmeQ ApproachPauli::generate_fct() in ", __file__, " verbosity = ", self.verbosity, " inputs: \n")
            print(f"State energies (E):", self.qd.Ea)
            #print(f"Tunneling amplitudes (Tba):", self.leads.Tba)
            print(f"Tunneling amplitudes (Tba) in pauli.py of QmeQ:\nLead 0:\n", (self.leads.Tba[0,:,:]).real, "\nLead 1:\n", (self.leads.Tba[1,:,:]).real) 
            #print(f"Tunneling amplitudes (Tba) in pauli.py of QmeQ:\n", (self.leads.Tba).real) 
            print(f"Chemical potentials (mulst):", self.leads.mulst)
            print(f"Temperatures (tlst):", self.leads.tlst)
            print(f"Band parameters (dlst):", self.leads.dlst)
            print(f"Number of charge states:", self.si.ncharge)
            print(f"States by charge (statesdm):", self.si.statesdm)
        
        #raise NotImplementedError("DEBUG: we exit here to make the debugging easier")
            
        E, Tba, si = self.qd.Ea, self.leads.Tba, self.si
        mulst, tlst, dlst = self.leads.mulst, self.leads.tlst, self.leads.dlst
        ncharge, nleads, statesdm = si.ncharge, si.nleads, si.statesdm

        itype = self.funcp.itype
        paulifct = self.paulifct
        #verb_print_(3,"ApproachPauli.generate_fct(): paulifct.shape = ", paulifct.shape)
        for charge in range(ncharge-1):
            ccharge = charge+1
            bcharge = charge
            for c, b in itertools.product(statesdm[ccharge], statesdm[bcharge]):
                cb = si.get_ind_dm1(c, b, bcharge)
                Ecb = E[c]-E[b]
                for l in range(nleads):
                    xcb = (Tba[l, b, c]*Tba[l, c, b]).real
                    rez = func_pauli(Ecb, mulst[l], tlst[l], dlst[l, 0], dlst[l, 1], itype)
                    paulifct[l, cb, 0] = xcb*rez[0]  # Forward
                    paulifct[l, cb, 1] = xcb*rez[1]  # Backward
                    #verb_print_(3,f"ApproachPauli.generate_fct() l: {l} i: {c} j: {b} cb: {cb} E_diff: {Ecb:.6f} coupling: {xcb:.6f} fermi: {rez[0]/(2*np.pi):.6f} factors:[{paulifct[l,cb,0]:.6f}, {paulifct[l,cb,1]:.6f}]")

        #raise NotImplementedError("DEBUG: we exit here to make the debugging easier")

    def generate_kern(self):
        """
        Generate Pauli master equation kernel.

        Parameters
        ----------
        kern : array
            (Modifies) Kernel matrix for Pauli master equation.
        """
        print("ApproachPauli.generate_kern()")
        if self.verbosity > 3:
            print("ApproachPauli.generate_kern() Building kernel matrix...")
            
        verb_print_(0,"ApproachPauli.generate_kern() ncharge: {}  statesdm: {}".format(self.si.ncharge, self.si.statesdm))
        si, kh = self.si, self.kernel_handler
        ncharge, statesdm = si.ncharge, si.statesdm

        verb_print_(0,"ApproachPauli.generate_kern() ncharge: {}  statesdm: {}".format(ncharge, statesdm))

        self.generate_fct()
        #if self.verbosity > 3:print("DEBUG: ApproachPauli.generate_kern() after generate_fct() kh.kern:\n", kh.kern)
        #exit(0)
        if self.verbosity > 3:  
            #print("ApproachPauli.generate_kern().1 after generate_fct() kh.kern:\n", kh.kern)
            print("ApproachPauli.generate_kern()paulifct  : \n", self.paulifct )
            print("lenlst    : ", self.si.lenlst )
            print("dictdm    : ", self.si.dictdm )
            print("shiftlst0 : ", self.si.shiftlst0)
            print("shiftlst1 : ", self.si.shiftlst1)
            print("mapdm0    : ", self.si.mapdm0 )
            raise ValueError("DEBUG TERMINATION")

        for bcharge in range(ncharge):
            for b in statesdm[bcharge]:
                if not kh.is_unique(b, b, bcharge):
                    continue
                self.generate_coupling_terms(b, b, bcharge)
        #debug_print("DEBUG: ApproachPauli.generate_kern() kh.kern:\n", kh.kern)

        if self.verbosity > 0:  
            print("ApproachPauli.generate_kern().2 after generate_coupling_terms() kh.kern:\n", kh.kern)
            np.set_printoptions(precision=15)
            print(kh.kern)
            np.set_printoptions(precision=5)  # Reset to default
        print("===== ApproachPauli.generate_kern() DONE ====="  )
        raise ValueError("DEBUG TERMINATION")

    def generate_coupling_terms(self, b, bp, bcharge):
        print(f"ApproachPauli.generate_coupling_terms() b: {b} bp: {bp}  bcharge: {bcharge} statesdm: {self.si.statesdm}", self.verbosity)
        """Generate coupling terms for the Pauli master equation."""
        verb = 3
        #if self.verbosity > 3:  print(f"\nDEBUG: generate_coupling_terms() state b: {b} bp: {bp}  bcharge: {bcharge} statesdm: {self.si.statesdm}")
        if self.verbosity > verb:  print(f"\nQmeQ pauli.py ======== generate_coupling_terms() state b: {b} Q: {bcharge}")
            
        #debug_print(f"DEBUG: ApproachPauli.generate_coupling_terms() b: {b}  bp: {bp}  bcharge: {bcharge} statesdm: {self.si.statesdm}")
        Approach.generate_coupling_terms(self, b, bp, bcharge)
        paulifct = self.paulifct

        #if self.verbosity > 3:  print(f"DEBUG generate_coupling_terms(): paulifct:{ paulifct }")

        si, kh = self.si, self.kernel_handler
        nleads, statesdm = si.nleads, si.statesdm

        acharge = bcharge-1
        ccharge = bcharge+1

        bb = si.get_ind_dm0(b, b, bcharge)
        
        # Handle transitions from lower charge states
        if self.verbosity > verb:  print(f"ApproachPauli.generate_coupling_terms() Q-1 states: ", statesdm[acharge] )
        for a in statesdm[acharge]:   # Loop over states with charge acharge = bcharge-1
            aa = si.get_ind_dm0(a, a, acharge)
            ba = si.get_ind_dm1(b, a, acharge)
            if self.verbosity > verb:  print(f"LOWER: state:{b} other:{a} aa:{aa} ba:{ba}")
            fctm, fctp = 0, 0
            for l in range(nleads):
                #if self.verbosity > 3:  print(f"DEBUG:   lead:{l} idx:{l},{ba}")
                fctm -= paulifct[l, ba, 1]  # Electron leaving
                fctp += paulifct[l, ba, 0]  # Electron entering
            #verb_print("LOWER  [%i,%i] fctm: %.6f fctp: %.6f     bb: %i aa: %i" %(b, a, fctm, fctp, bb, aa ) ); 
            kh.set_matrix_element_pauli(fctm, fctp, bb, aa)
            if self.verbosity > verb:   print(f"ApproachPauli.generate_coupling_terms() state:{b} other:{a} rate:{fctp:.6f}")
        
        # Handle transitions to higher charge states
        if self.verbosity > verb:  print(f"ApproachPauli.generate_coupling_terms() Q+1 states", statesdm[ccharge] )
        for c in statesdm[ccharge]: # Loop over states with charge ccharge = bcharge+1
            cc = si.get_ind_dm0(c, c, ccharge)
            cb = si.get_ind_dm1(c, b, bcharge)
            if self.verbosity > verb: print(f"HIGHER: state:{b} other:{c} cc:{cc} cb:{cb}")
            fctm, fctp = 0, 0
            for l in range(nleads):
                #if self.verbosity > 3: print(f"DEBUG:   lead:{l} idx:{l},{cb}")
                fctm -= paulifct[l, cb, 0]  # Electron entering
                fctp += paulifct[l, cb, 1]  # Electron leaving
            #verb_print("HIGHER [%i,%i] fctm: %.6f fctp: %.6f      bb: %i cc: %i" %(b, c, fctm, fctp,    bb, cc ) ); 
            kh.set_matrix_element_pauli(fctm, fctp, bb, cc)
            if self.verbosity > verb: print(f"ApproachPauli.generate_coupling_terms() state:{b} other:{c} rate:{fctp:.6f}")

        #if self.verbosity > 3: print(f"DEBUG: generate_coupling_terms() in {__file__}, kh.kern:\n", kh.kern)

    def generate_current(self):
        """
        Calculates currents using Pauli master equation approach.

        Parameters
        ----------
        current : array
            (Modifies) Values of the current having nleads entries.
        energy_current : array
            (Modifies) Values of the energy current having nleads entries.
        heat_current : array
            (Modifies) Values of the heat current having nleads entries.
        """
        if self.verbosity > 3:
            print("\nDEBUG: generate_current() Calculating currents...")
            
        debug_print("DEBUG: ApproachPauli.generate_current()")
        phi0, E, paulifct, si = self.phi0, self.qd.Ea, self.paulifct, self.si
        ncharge, nleads, statesdm = si.ncharge, si.nleads, si.statesdm

        current = self.current
        energy_current = self.energy_current

        for charge in range(ncharge-1):
            ccharge = charge+1
            bcharge = charge
            for c in statesdm[ccharge]:
                cc = si.get_ind_dm0(c, c, ccharge)
                for b in statesdm[bcharge]:
                    bb = si.get_ind_dm0(b, b, bcharge)
                    cb = si.get_ind_dm1(c, b, bcharge)
                    for l in range(nleads):
                        fct1 = +phi0[bb]*paulifct[l, cb, 0]
                        fct2 = -phi0[cc]*paulifct[l, cb, 1]
                        current[l] += fct1 + fct2
                        energy_current[l] += -(E[b]-E[c])*(fct1 + fct2)

        self.heat_current[:] = energy_current - current*self.leads.mulst
# ---------------------------------------------------------------------------------------------------
