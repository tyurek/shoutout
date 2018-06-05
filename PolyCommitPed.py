from charm.toolbox.pairinggroup import PairingGroup,ZR,G1,G2,GT,pair
from base64 import encodestring, decodestring
import random
from helperfunctions import *
from math import floor

#group = PairingGroup('SS512')
#group = PairingGroup('MNT159')
#group = PairingGroup('MNT224')


class PolyCommitPed:
    def __init__ (self, t, pk, group, symflag, minimizepairings=False, seed=None):
        self.group = group
        self.pk = pk
        self.g = pk[0]
        #if str(self.group.groupType()) == 'SS512' or str(self.group.groupType()) == 'SS1024' or str(self.group.groupType()) == 'SS1536':
        if symflag:
            self.h = pk[t+1]
            self.symmetric = True
        #elif str(self.group.groupType()) == 'MNT159' or str(self.group.groupType()) == 'MNT224' or str(self.group.groupType()) == 'BN256':
        elif not symflag:
            self.h = pk[t+3]
            self.symmetric = False
        else:
            print "Error: Invalid Curve Specified"
            return
        self.t = t
        self.seed = seed
        self.ONE = group.random(ZR)*0+1
        self.c = self.ONE
        self.minimizepairings = minimizepairings
        if minimizepairings:
            self.lhss = {}
        if symflag:
            self.gg = self.group.pair_prod(self.g, self.g)
            self.gh = self.group.pair_prod(self.g, self.h)
        else:
            self.gg = self.group.pair_prod(self.g, pk[t+1])
            self.gh = self.group.pair_prod(self.h, pk[t+1])

    def commit (self, poly, secretpoly):
        #secretpoly is the polynomial phi hat which is used to make the polynomial commitment
        if self.symmetric:
            fudge = 0
        else:
            fudge = 2
        c = self.ONE
        i = 0
        for item in self.pk[:self.t+1]:
            c *= item ** poly[i]
            i += 1
        i = 0
        for item in self.pk[self.t+1+fudge:]:
            c *= item ** secretpoly[i]
            i += 1
        #c should be equivalent to (self.g **(f(poly, self.alpha))) * (self.h **(f(self.secretpoly, self.alpha)))
        return c

    #def open (self):
    #    return {'c': self.c, 'poly': self.poly, 'secretpoly': self.secretpoly} 

    def verify_poly (self, c, poly, secretpoly):
        if self.symmetric:
            fudge = 0
        else:
            fudge = 2
        tempc = self.ONE
        i = 0
        for item in self.pk[:self.t+1]:
            tempc *= item ** poly[i]
            i += 1
        i = 0
        for item in self.pk[self.t+1+fudge:]:
            tempc *= item ** secretpoly[i]
            i += 1
        return c == tempc

    def create_witness(self, poly, secretpoly, i):
        if self.symmetric:
            psi = polynomial_divide([poly[0] - f(poly,i)] + poly[1:], [self.ONE*i*-1,self.ONE])
            psihat = polynomial_divide([secretpoly[0] - f(secretpoly,i)] + secretpoly[1:], [self.ONE*i*-1,self.ONE])
            witness = self.ONE
            j = 0
            for item in self.pk[:self.t]:
                witness *= item ** psi[j]
                j += 1
            j = 0
            #funky indexing is due to the structure of pk and that g^alpha^t and h^alpha^t aren't needed
            for item in self.pk[self.t+1:self.t+1 + self.t]:
                witness *= item ** psihat[j]
                j += 1
            #witness should be equivalent to (self.g **(f(psi, self.alpha))) * (self.h **(f(psihat, self.alpha)))
            return witness
        else:
            psi = polynomial_divide([poly[0] - f(poly,i)] + poly[1:], [self.ONE*i*-1,self.ONE])
            psihat = polynomial_divide([secretpoly[0] - f(secretpoly,i)] + secretpoly[1:], [self.ONE*i*-1,self.ONE])
            witness = self.ONE
            j = 0
            for item in self.pk[:self.t]:
                witness *= item ** psi[j]
                j += 1
            j = 0
            for item in self.pk[self.t+3:self.t+3 + self.t]:
                witness *= item ** psihat[j]
                j += 1
            return witness
    
    #If reusing the same commitment, the lhs of the comparison will be the same. Take advantage of this to save pairings
    def verify_eval(self, c, i, polyeval, secretpolyeval, witness):
        if self.symmetric:
            if self.minimizepairings:
                if str(c) in self.lhss:
                    lhs = self.lhss[str(c)]
                else:
                    lhs =  self.group.pair_prod(c, self.g)
                    self.lhss[str(c)] = lhs
            else:
                lhs =  self.group.pair_prod(c, self.g)
            #rhs = self.group.pair_prod(witness, self.pk[1] / (self.g ** i)) * self.group.pair_prod(self.g**polyeval * self.h**secretpolyeval, self.g)
            #super awesome pairing optimization
            rhs = self.group.pair_prod(witness, self.pk[1] / (self.g ** i)) * self.gg**polyeval * self.gh**secretpolyeval
            return lhs == rhs
        else:
            #self.pk[self.t + 1] is ghat in G2
            #lhs =  self.group.pair_prod(c, self.pk[self.t + 1])
            if self.minimizepairings:
                if str(c) in self.lhss:
                    lhs = self.lhss[str(c)]
                else:
                    lhs =  self.group.pair_prod(c, self.pk[self.t + 1])
                    self.lhss[str(c)] = lhs
            else:
                lhs =  self.group.pair_prod(c, self.pk[self.t + 1])
            #rhs = self.group.pair_prod(witness, self.pk[self.t + 2] / (self.pk[self.t + 1] ** i)) * self.group.pair_prod(self.g**polyeval * self.h**secretpolyeval, self.pk[self.t + 1])
            rhs = self.group.pair_prod(witness, self.pk[self.t + 2] / (self.pk[self.t + 1] ** i)) * self.gg**polyeval * self.gh**secretpolyeval
            return  lhs == rhs

    #This could be easily rewritten to account for multiple i's, but was not needed in this usecase
    def batch_verify_eval(self, cs, i, polyevals, secretpolyevals, witnesses):
        cprod = cs[0] / cs[0]
        for c in cs:
            cprod = cprod * c
        witnessprod = witnesses[0] / witnesses[0]
        for witness in witnesses:
            witnessprod = witnessprod * witness
        #Oh Charm...
        witnesspownegiprod = (witnessprod ** (i)) ** (-1)
        polyevalsum = polyevals[0] * 0
        for polyeval in polyevals:
            polyevalsum = polyevalsum + polyeval
        secretpolyevalsum = secretpolyevals[0] * 0
        for secretpolyeval in secretpolyevals:
            secretpolyevalsum = secretpolyevalsum + secretpolyeval
        if self.symmetric:
            lhs =  self.group.pair_prod(cprod, self.g)
            #rhs = self.group.pair_prod(witness, self.pk[1] / (self.g ** i)) * self.group.pair_prod(self.g**polyeval * self.h**secretpolyeval, self.g)
            #super awesome pairing optimization
            rhs = self.group.pair_prod(witnessprod, self.pk[1]) * self.group.pair_prod(witnesspownegiprod, self.g) * self.gg**polyevalsum * self.gh**secretpolyevalsum
            return lhs == rhs
        else:
            lhs =  self.group.pair_prod(cprod, self.pk[self.t + 1])
            #rhs = self.group.pair_prod(witness, self.pk[self.t + 2] / (self.pk[self.t + 1] ** i)) * self.group.pair_prod(self.g**polyeval * self.h**secretpolyeval, self.pk[self.t + 1])
            rhs = self.group.pair_prod(witnessprod, self.pk[self.t + 2]) * self.group.pair_prod(witnesspownegiprod, self.pk[self.t + 1]) * self.gg**polyevalsum * self.gh**secretpolyevalsum
            return  lhs == rhs

    #Returns an array of Trues and Falses corresponding to if each input set was valid or not
    def find_valid_evals(self, cs, i, polyevals, secretpolyevals, witnesses):
        assert len(cs) == len(polyevals) and len(cs) == len(secretpolyevals) and len(cs) == len(witnesses)
        validsets = [True]*len(cs)
        queue = []
        indexes = []
        for j in range(len(cs)):
            indexes.append(j)
        queue.append([cs, polyevals, secretpolyevals, witnesses, indexes])
        #Would have been easy to do this recursively, but this is *hopefully* more performant
        while len(queue) > 0:
            if self.batch_verify_eval (queue[0][0], i, queue[0][1], queue[0][2], queue[0][3]):
                queue.pop(0)
            else:
                if len(queue[0][0]) == 1:
                    validsets[queue[0][4][0]] = False
                    queue.pop(0)
                else:
                    #index at which to divide the lists in two
                    sind = floor(len(queue[0][0]))
                    queue.append([queue[0][0][0:sind], queue[0][1][0:sind], queue[0][2][0:sind], queue[0][3][0:sind], queue[0][4][0:sind]])
                    queue.append([queue[0][0][sind:], queue[0][1][sind:], queue[0][2][sind:], queue[0][3][sind:], queue[0][4][sind:]])
                    queue.pop(0)
        return validsets
                

    #If you're going to use the same PolyCommit object for multiple AVSS runs, use this after each run so you aren't carying around extra pairings
    def flush(self):
        self.lhss = {}
