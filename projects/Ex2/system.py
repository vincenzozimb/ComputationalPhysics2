import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import scipy as sp

from basis import *
        
class SYSTEM(STO4G): 
    def __init__(self, name, N, Z, orbitals):
        self.name=name
        self.N=N
        self.Z=Z
        self.orbitals = orbitals
        self.basis = STO4G(self)
        self.h= self.kinetic() + self.ext_pot()
        self.DE = self.de()
    
    def kinetic(self):
        """
        Compute kinetic energy matrix.
        INPUT:
            SYS: system
        OUTPUT:
            T: Kinetic energy matrix
        """
        T = np.zeros((len(self.basis.STO4G_alpha),len(self.basis.STO4G_alpha)))
        for p in range(len(self.basis.STO4G_alpha)):
            for q in range(len(self.basis.STO4G_alpha)):
                T[p,q] = 3*(self.basis.STO4G_alpha[p]*self.basis.STO4G_alpha[q]*(np.pi**1.5))/(self.basis.STO4G_alpha[p]+self.basis.STO4G_alpha[q])**(2.5)
        return T

    def ext_pot(self):
        """
        Compute external potential matrix.
        INPUT:
            SYS: system
        OUTPUT:
            V: External potential matrix
        """
        V = np.zeros((len(self.basis.STO4G_alpha),len(self.basis.STO4G_alpha)))
        for p in range(len(self.basis.STO4G_alpha)):
            for q in range(len(self.basis.STO4G_alpha)):
                V[p,q] = (2*np.pi/(self.basis.STO4G_alpha[p]+self.basis.STO4G_alpha[q]))
        return V
    def de(self):
        """
        Compute the direct term in Hartree-Fock.
        INPUT:
            SYS: system
        OUTPUT:
            H: Direct term in Hartree-Fock V(p,r,q,s)=<pr|V|qs> + the exchange term  V(p,r,q,s)=<pr|V|sq>
        """
        alpha = self.basis.STO4G_alpha
        D = np.zeros((len(alpha),len(alpha),len(alpha),len(alpha)))
        E = D 
        ## initialize the diagonal terms
        for p in range(0,len(alpha)):
            for q in range(0,len(alpha)):
                for r in range(0,len(alpha)):
                    for s in range(0,len(alpha)):
                        D[p,r,q,s] = (2*(np.pi)**(2.5))/((alpha[p]+alpha[q])*(alpha[r]+alpha[s])*(np.sqrt(alpha[p]+alpha[q]+alpha[r]+alpha[s])))
                        E[p,r,s,q] = (2*(np.pi)**(2.5))/((alpha[p]+alpha[q])*(alpha[r]+alpha[s])*(np.sqrt(alpha[p]+alpha[q]+alpha[r]+alpha[s])))
        return D + E 

    def fill_density_matrix(self, C): ####### CHECK THIS
        """
        Compute the density matrix.
        INPUT:
            SYS: system
            C: AO coefficients
        OUTPUT:
            D: Density matrix
        """
        density = np.outer(np.conjugte(C), C)
        return density

    def fill_fock_matrix(self, C): ####### CHECK THIS
        """
        Compute the Fock matrix.
        INPUT:
            SYS: system
            C: AO coefficients
        OUTPUT:
            F: Fock matrix
        """
        F = self.h
        F += 0.5*np.einsum('pqrs,rs->pq', self.DE, C) 
        F = np.einsum('pqrs,rs->pq', self.h, C) + np.einsum('pqrs,pr->qs', self.h, C)
        return F


    def solve_hf(self): ##### CHECK THIS
        for i in range (0,round(self.N/2)):
            epsilon, C = sp.linalg.eigh(np.identity(len(self.basis.STO4G_alpha)))




sys = SYSTEM("H",1,1,["1s"])