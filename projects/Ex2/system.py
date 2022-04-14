import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import scipy.linalg
import argparse
import sys
from basis import *
        
#instruction how to run the program 
##################################
#SELECTOR FROM COMMAND LINE
##################################
ap = argparse.ArgumentParser(description="The program performs the calculation of the ground state energy and the ground state wave function for a given atom.")
ap.add_argument("-a", "--Atom choice", choices=["H", "He", "Be"], type=str, required=True, help=" H--> Hartree fock for close shell Hydrogen atom; \n \
#     He--> Hartree fock for Helium atom; \
#     Be--> Hartree fock for Berillium atom.")
args = vars(ap.parse_args())
#print(sys.argv)
selector = args["Atom choice"]
#selector =str(sys.argv[1])
print(selector)

#########################################
#definition of functions
#########################################
class SYSTEM(STO4G): 
    def __init__(self, name):
        self.name=name
        self.imax = 10000
        self.basis = STO4G(name)
        self.orbitals = self.basis.orbitals
        self.N = self.basis.N
        self.Z = self.basis.Z
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
        E = np.zeros((len(alpha),len(alpha),len(alpha),len(alpha)))
        ## initialize the diagonal terms
        for p in range(0,len(alpha)):
            for q in range(0,len(alpha)):
                for r in range(0,len(alpha)):
                    for s in range(0,len(alpha)):
                        D[p,r,q,s] = (2*(np.pi)**(2.5))/((alpha[p]+alpha[q])*(alpha[r]+alpha[s])*(np.sqrt(alpha[p]+alpha[q]+alpha[r]+alpha[s])))
                        E[p,r,s,q] = (2*(np.pi)**(2.5))/((alpha[p]+alpha[s])*(alpha[r]+alpha[q])*(np.sqrt(alpha[p]+alpha[q]+alpha[r]+alpha[s])))
        return D - E 

    def fill_density_matrix(self, C): 
        """
        Compute the density matrix.

        Returns the density matrix evaluated using the coefficient matrix C, according to
        
        :math:`M_{\alpha\beta} = \sum_i^{occ} C_{\alpha,i}^* C_{\beta,i}`
        i runs over the occupied orbitals 
        alpha and beta over the STO4G coefficients

        INPUT:
            SYS: system
            C: AO coefficients
        OUTPUT:
            D: Density matrix
        """
        density = np.zeros( (4,4) , dtype=np.complex128)
        for i in range(round(self.N/2+0.1)):
            density += np.outer( np.conjugate(C[:,i]), C[:,i])
        return density

    def fill_fock_matrix(self, C): 
        """
        Compute the Fock matrix.
        INPUT:
            SYS: system
            C: AO coefficients
        OUTPUT:
            F: Fock matrix
        """
        f = np.zeros(self.h.shape, dtype=np.complex128)    
        density = self.fill_density_matrix(C)
        f += np.einsum('ij,aibj->ab', density, self.DE, dtype=np.complex128)
        f += self.h
        return f


    def solve_hf(self): ##### CHECK THIS
        r"""
        Solve the Hartree-Fock equation, with initial guess on the coefficients.

        INPUT:
            SYS: system
        OUTPUT:
            C: AO coefficients
            Energy: Hartree-Fock energy
        """
        # Initial guess
        r = round(self.N/2+0.1)  
        C =np.zeros((4, r))
        epsilon = np.zeros((r))
        for i in range(0,r):
                C[0,i] = 1
                epsilon[i] = 1
        # density = self.fill_density_matrix(C)
        # F = self.fill_fock_matrix(C)
        # print(F)
        #epsilon, C = scipy.linalg.eigh(np.identity(r))
        epsilon_old = np.zeros(epsilon.shape)

        #density = self.fill_density_matrix(C)
        l=10e-6
        i=0
        while np.linalg.norm(epsilon-epsilon_old) > l and i<self.imax:
            epsilon_old = epsilon
            F = self.fill_fock_matrix(C)
            epsilon, C = scipy.linalg.eig(F, self.basis.overlap(),check_finite=True)
            i +=1
        ## Print coefficients and energies
        print("i=",i)
        print("Overlap:")
        print(self.basis.overlap())
        print("Coefficients:")
        print(C,"\n")
        print("Energies:")
        print(epsilon)



#################################################################
#MAIN
#################################################################
#print(selector)
sys = SYSTEM(selector)
sys.solve_hf()