import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import scipy as sp


class STO4G:
    def __init__(self,name):
        self.name = name
        self.K = 4
        if self.name == "H":
            self.alpha = np.array([13.00773,1.962079,0.444529,0.1219492])
            self.orbitals = ["1s"]
            self.N = 1
            self.Z = 1
        if self.name == "He":
            self.alpha = np.array([14.899983,2.726485,0.757447,0.251390])
            self.orbitals = ["1s"]
            self.N=2
            self.Z=2
        if self.name == "Be":
            self.alpha =np.array([70.64859542,12.92782254,3.591490662,1.191983464,3.072833610,0.6652025433,0.2162825386,0.08306680972])
            self.orbitals = ["1s","2s"]
            self.N = 4
            self.Z = 4
        self.STO4G_alpha= []
        for o in self.orbitals:
            if o == "1s":
                self.STO4G_alpha = np.concatenate((self.STO4G_alpha, self.alpha[0:self.K]))
            if o == "2s":
                self.STO4G_alpha = np.concatenate((self.STO4G_alpha, self.alpha[self.K:2*self.K]))
    
    def overlap(self):
        """
        Compute overlap matrix S.
        INPUT:
            BASIS: basis set
        OUTPUT:
            S: Overlap matrix
        """
        S = np.zeros((len(self.STO4G_alpha),len(self.STO4G_alpha)))
        for i in range(len(self.STO4G_alpha)):
            for j in range(len(self.STO4G_alpha)):
                S[i,j] = (np.pi/(self.STO4G_alpha[i]+self.STO4G_alpha[j]))**1.5
        return S
    
    def gaussian(self, x, alpha):
        """
        Compute vector of gaussian function"""
        self.gg = []
        for a in alpha:
            self.gg.append(np.exp(-a*x**2))
        return np.array(self.gg)

    def info(self):
        """
        Print informations about the bais set.
        """
        print("########################")
        print("STO-4G MINIMAL BASIS SET")
        print("########################\n")
###############################################################################
#           MAIN                                                              #
###############################################################################