import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import scipy.linalg
import argparse
from basis import *

#### Define parameters 
beta=10e-3


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
        self.imax = 10000
        self.basis = STO4G(name)
        self.orbitals = self.basis.orbitals
        self.N = self.basis.N
        self.Z = self.basis.Z
        self.h= self.kinetic() + self.ext_pot()
        self.DE = self.de()
        self.threshold=10e-8
    
    def kinetic(self):
        """
        Compute kinetic energy matrix.
        INPUT:
            SYS: system
        OUTPUT:
            T: Kinetic energy matrix
        """
        print("Basis set:",self.basis.STO4G_alpha)
        T = np.zeros((len(self.basis.STO4G_alpha),len(self.basis.STO4G_alpha)))
        for p in range(len(self.basis.STO4G_alpha)):
            for q in range(len(self.basis.STO4G_alpha)):
                T[p,q] = self.basis.STO4G_alpha[p]*self.basis.STO4G_alpha[q]/((self.basis.STO4G_alpha[p]+self.basis.STO4G_alpha[q])**(2.5))
        T = (3*(np.pi**1.5))*T
        return T

    def ext_pot(self):
        """This funxctiommhlnl

        Returns:
            _type_: _description_
        """
        V = np.zeros((len(self.basis.STO4G_alpha),len(self.basis.STO4G_alpha)))
        for p in range(len(self.basis.STO4G_alpha)):
            for q in range(len(self.basis.STO4G_alpha)):
                V[p,q] = 1/(self.basis.STO4G_alpha[p]+self.basis.STO4G_alpha[q])
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
        alpha = self.basis.STO4G_alpha
        D = np.zeros((len(alpha),len(alpha),len(alpha),len(alpha)))
        E = np.zeros((len(alpha),len(alpha),len(alpha),len(alpha)))
        ## initialize the direct and exchange terms for only occupied orbitals and elements different from the hydrogen atom
        if self.name!="H":
            for p in range(0,len(alpha)):
                for q in range(0,len(alpha)):
                    for r in range(0,len(alpha)):
                        for s in range(0,len(alpha)):
                            if self.name == "He": #### ??????? Also exchange element in helium?
                                D[p,q,r,s] = 1/((alpha[p]+alpha[r])*(alpha[q]+alpha[s])*(np.sqrt(alpha[p]+alpha[q]+alpha[r]+alpha[s])))
                                E[p,q,r,s] = 1/((alpha[p]+alpha[s])*(alpha[r]+alpha[q])*(np.sqrt(alpha[p]+alpha[q]+alpha[r]+alpha[s])))

                            else: 
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
        density = np.zeros( (len(self.basis.STO4G_alpha),len(self.basis.STO4G_alpha)) , dtype=np.complex128)
        for i in range(round(self.N/2+0.1)):
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
        f += np.einsum('qs,pqrs->pr', density, self.DE, dtype=np.complex128)
        return f


    def norm(self):
        """
        Compute the norm of the density matrix.
        INPUT:
            SYS: system
            wf: wave function
        OUTPUT:
            norm: Norm of the wavefunction
        """

        return 0
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
            f = np.einsum('bd,abcd->ac', self.fill_density_matrix(self.C), self.DE, dtype=np.complex128)
            energy_groundstate += np.einsum('pr,pr', self.fill_density_matrix(self.C), f)## Direct and exchange term contribution
        print("energy_groundstate=",energy_groundstate)


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
        r = round(self.N/2+0.1)  
        C =np.zeros((len(self.basis.STO4G_alpha), r))
        epsilon = np.ones((len(self.basis.STO4G_alpha)))
        for i in range(0,r):
                C[0,i] = 1
        epsilon_old = np.zeros(epsilon.shape)
        i=0
        self.C = C
        F = self.fill_fock_matrix(C)
        #print("C=",self.C,"size=",self.C.shape)
        while np.linalg.norm(min(epsilon)-min(epsilon_old)) > self.threshold and i < self.imax:
            #Store the old eigenvalues and coefficients
            epsilon_old = epsilon
            C_old=self.C
            #Compute the new eigenvalues and coefficients
            epsilon, C = scipy.linalg.eig(F, self.basis.overlap(),check_finite=True)
            #Select only the lowest eigenvalue and get the index of the lowest eigenvector
            self.epsilon_min=min(epsilon)
            self.epsilon_index=np.where(epsilon==self.epsilon_min)[0]
            #Store the new coefficients in the coefficient instance of the class
            self.C=C[:,self.epsilon_index]
            #normalize the wf coefficients
            self.C=self.C/np.sqrt(np.einsum('ij,jk->ik',np.transpose(self.C), np.einsum('ij,jk->ik',self.basis.overlap(),self.C)))
            #Vary slowly the coefficients (! then they are not normalized anymore)
            self.C = self.C*beta + (1-beta)*C_old
            #self.C=self.C/np.sqrt(np.einsum('ij,jk->ik',np.transpose(self.C), np.einsum('ij,jk->ik',self.basis.overlap(),self.C)))
            #Compute the new Fock matrix
            F = self.fill_fock_matrix(self.C)
            i +=1
        #get the minimun eigenvalue in order to calculate the ground state energy: 
        self.epsilon_min=min(epsilon)
        self.epsilon_index=np.where(epsilon==self.epsilon_min)[0]
        self.C=C[:,self.epsilon_index]
        #normalize the wf
        self.C=self.C/np.sqrt(np.einsum('ij,jk->ik',np.transpose(self.C), np.einsum('ij,jk->ik',self.basis.overlap(),self.C)))
        #Compute the ground state energy
        self.evaluate_ground_state_energy()


    def plot_wf(self):
        r"""
        Plot the wavefunction.
            INPUT:
                SYS: system
        """
        r = np.linspace(0.001,6,1000)
        s = self.basis.gaussian(r, self.basis.STO4G_alpha)
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
    syst.plot_wf()
    #selector =str(sys.argv[1])
####################################