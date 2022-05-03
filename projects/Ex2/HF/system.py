import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import scipy.linalg
import argparse
from basis import *

#### Define parameters 
beta=10e-2
"""
Hartree-Fock for closed shell atoms 
----------------------------------------------------
"""
#########################################
#definition of functions and classes
#########################################
class SYSTEM(STO4G): 
    """_summary_

    Args:
        STO4G (basis set): gaussian basis set
    """
    def __init__(self, name):
        """Constructor for the system class.

        Args:
            name (str): name of the atom
        Instance variables:
            :h (matrix): one-electron integrals
            :DE (matrix): direct and exchange term
            _basis (basis set): basis set
            _Z (float): atomic number
            -name (str): name of the atom
        """
        self.name=name
        self.imax = 50000
        self.basis = STO4G(name)
        self.orbitals = self.basis.orbitals
        self.N = self.basis.N
        self.Z = self.basis.Z
        self.h= self.kinetic() + self.ext_pot()
        self.DE = self.de()
        self.threshold=10e-5
    
    def kinetic(self):
        """
        Compute kinetic energy matrix.
        INPUT:
            SYS: system
        OUTPUT:
            T: Kinetic energy matrix
        """
        print("Basis set:",self.basis.alpha)
        T = np.zeros((len(self.basis.alpha),len(self.basis.alpha)))
        for p in range(len(self.basis.alpha)):
            for q in range(len(self.basis.alpha)):
                T[p,q] = self.basis.alpha[p]*self.basis.alpha[q]/((self.basis.alpha[p]+self.basis.alpha[q])**(2.5))
        T = (3*(np.pi**1.5))*T
        return T

    def ext_pot(self):
        """This funxctiommhlnl

        Returns:
            _type_: _description_
        """
        V = np.zeros((len(self.basis.alpha),len(self.basis.alpha)))
        for p in range(len(self.basis.alpha)):
            for q in range(len(self.basis.alpha)):
                V[p,q] = 1/(self.basis.alpha[p]+self.basis.alpha[q])
        V = -(self.Z*(2*np.pi))*V
        return V

    def de(self): 
        """
        Compute the direct term in Hartree-Fock.
        INPUT:
            SYS: system
        OUTPUT:
            D-E: Direct term in Hartree-Fock V(p,r,q,s)=<pr|V|qs> + the exchange term  V(p,r,q,s)=<pr|V|sq> 
        """
        alpha = self.basis.alpha
        D = np.zeros((len(alpha),len(alpha),len(alpha),len(alpha)))
        E = np.zeros((len(alpha),len(alpha),len(alpha),len(alpha)))
        ## initialize the direct and exchange terms for only occupied orbitals and elements different from the hydrogen atom
        if self.name!="H":
            for p in range(0,len(alpha)):
                for q in range(0,len(alpha)):
                    for r in range(0,len(alpha)):
                        for s in range(0,len(alpha)):
                            D[p,q,r,s] = 1/((alpha[p]+alpha[r])*(alpha[q]+alpha[s])*(np.sqrt(alpha[p]+alpha[q]+alpha[r]+alpha[s])))
                            E[p,q,r,s] = 1/((alpha[p]+alpha[s])*(alpha[r]+alpha[q])*(np.sqrt(alpha[p]+alpha[q]+alpha[r]+alpha[s])))
            D = (2*(np.pi)**(2.5))*D
            E = (2*(np.pi)**(2.5))*E
        return (2*D-E)

    def fill_density_matrix(self,C): 
        r"""
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
        density = np.zeros( (len(self.basis.alpha),len(self.basis.alpha)) , dtype=np.complex128)
        for i in range(len(self.orbitals)):
            density += np.outer( np.conjugate(C[:,i]), C[:,i])
        return density

    def fill_fock_matrix(self,C): 
        """Compute the Fock matrix.

        Args:
            C (array): AO coefficients

        Returns:
            matrix: Fock matrix
        """
        f = np.zeros(self.h.shape, dtype=np.complex128)   
        density = self.fill_density_matrix(C)
        f += self.h
        f += np.einsum('qs,pqrs->pr', density, self.DE)
        return f

    def evaluate_ground_state_energy(self):
        """
        Evaluate the ground state energy.

        INPUT(self):
            SYS: system
            The eigenvalue problem of the Fock matrix should be already solved.
            The system must have a value epsilon (in principle the minimum eigenvalue) and the coefficients are already implemented in the system.
            C: AO coefficients

        OUTPUT:
            Energy: Hartree-Fock energy
        """
        if(self.name=="H"):
            energy_groundstate=self.epsilon_min
            #energy_groundstate = np.einsum('pr,pr', self.fill_density_matrix(self.C), self.h)
            print('Energy of the ground state (matrix element):',energy_groundstate, "compared with energy eignvalue:" , self.epsilon_min)
        else:
            energy_groundstate = 2*np.einsum('pr,pr', self.fill_density_matrix(self.C), self.h) ## Lowest eigenvalue of the Fock matrix
            f = np.einsum('bd,abcd->ac', self.fill_density_matrix(self.C), self.DE)
            energy_groundstate += np.einsum('pr,pr', self.fill_density_matrix(self.C), f)## Direct and exchange term contribution
            f1 = np.einsum('qs,pqrs->pr', self.fill_density_matrix(self.C), self.DE)
            energy_groundstate2 = 2*self.epsilon_min.sum() - np.einsum('pr,pr', self.fill_density_matrix(self.C), f1)
            print("energy_groundstate=",energy_groundstate, "energy_groundstate2=",energy_groundstate2)


    def solve_hf(self):
        r"""
        Solve the Hartree-Fock equation, with initial guess on the coefficients.
        
        The coefficient matrix is initialized as follow:
            column index = run over the orbitals
            row index = run over the coefficients on each gaussian basis function (STO4G) 
        As initial guess we state that the orbital wave function is equal to the first gaussian in the basis expansion,
         i.e. :math:'\phi_i =C_{1i}exp(alpha_1 r^2)' and C_{1i}=1 for the first gaussian.

        INPUT:
            SYS: system
        OUTPUT:
            self.C: AO coefficients
            self.energy_min: minimum energy eigenvalue after the Hartree-Fock procedure
        """
        #First guess on the coefficients
        r = len(self.orbitals)  
        C =np.zeros((len(self.basis.alpha), r))
        epsilon = np.ones((len(self.basis.alpha)))
        for i in range(0,r):
                C[i*4,i] = 1
                C[:,i] = C[:,i]/np.sqrt(np.einsum('i,i',np.transpose(C[:,i]), np.einsum('ij,j->i',self.basis.overlap1(),C[:,i])))
        epsilon_old = np.zeros(epsilon.shape)
        i=0
        self.C = C
        converg = 1
        F = self.fill_fock_matrix(self.C)
        #print("C=",self.C,"size=",self.C.shape)
        while  converg > self.threshold and i < self.imax:
            converg = 0
            #Store the old eigenvalues and coefficients
            epsilon_old = epsilon
            self.C_old = self.C
            epsilon , C = self.Rothan_eq(F,self.basis.overlap1())
            self.C=C[:,0:r]
            print("self.C=",self.C,"size=",self.C.shape)
            #print("epsilon=",epsilon,"shape=",epsilon.shape)
            epsilon = epsilon[0:r]
            #Vary slowly the coefficients (! then they are not normalized anymore)
            self.C = self.C*beta + (1-beta)*self.C_old
            for k in range (0,r):
                #normalize the coefficients
                converg += (abs(epsilon[k]-epsilon_old[k]))
                self.C[:,k] = self.C[:,k]/np.sqrt(np.einsum('i,i',np.transpose(self.C[:,k]), np.einsum('ij,j->i',self.basis.overlap1(),self.C[:,k])))
            F = self.fill_fock_matrix(self.C)
            i +=1
        self.epsilon_min=epsilon
        #Compute the ground state energy
        print("number of cycles:",i)
        self.evaluate_ground_state_energy()


    def Rothan_eq(self,F, S):
        ei, C = scipy.linalg.eig(F, S)
    # sort eigvalue and eigvector from lower to higher
        idx = ei.argsort()[::1]   
        ei = ei[idx]
        C = C[:,idx]
    # Normalize the coefficients
        #for i in range(len(self.orbitals)):
         #   C[:,i] = C[:,i]/np.sqrt(np.einsum('k,k',C[:,i], np.einsum('ij,j->i',S,C[:,i])))
        return ei, C
    def plot_wf(self):
        r"""
        Plot the wavefunction.
            INPUT:
                SYS: system
        """
        r = np.linspace(0.001,6,1000)
        s = self.basis.gaussian(r, self.basis.alpha)
        wf = -np.einsum('ij,ik->jk', self.C, s)
        plt.plot(r,wf[0,:])
        plt.legend(["Radial wavefunction of "+self.name])
        plt.xlabel(r'$r/a_0$')
        plt.ylabel(r'$R(r)/a_0^{3}$')
        plt.show()



#################################################################
#MAIN
#################################################################
####################################
#instruction how to run the program 
##################################
#SELECTOR FROM COMMAND LINE
##################################
if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="The program performs the calculation of the ground state energy and the ground state wave function for a given atom.")
    ap.add_argument("-a", "--Atom choice", choices=["H", "He", "Be"], type=str, required=True, help=" H--> Hartree fock for close shell Hydrogen atom; \n \
        He--> Hartree fock for Helium atom; \
        Be--> Hartree fock for Berillium atom.")
    args = vars(ap.parse_args())
    selector = args["Atom choice"]
    syst = SYSTEM(selector)
    syst.solve_hf()
    #syst.plot_wf()
    #selector =str(sys.argv[1])
####################################