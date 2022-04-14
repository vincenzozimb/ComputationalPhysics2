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
            self.alpha = np.array([0.5,0.5,0.5,0.5])
            self.orbitals = ["1s"]
            self.N=2
            self.Z=2
        if self.name == "Be":
            self.alpha =np.array([0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5])
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