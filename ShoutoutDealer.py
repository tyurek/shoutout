from PolyCommitPed import *
import math
from helperfunctions import *
from charm.toolbox.pairinggroup import PairingGroup,ZR,G1,G2,GT,pair

class ShoutoutDealer:
    #def __init__ (self, pk, secret, position, serverids, group, matrixdim, recv_function, send_function):
    def __init__ (self, pk, secret, position, serverids, group, matrixdim, send_function, recv_function):
        #secret is a list of field elements that can be mapped to a message
        #The database is assumed to be a matrixdim X matrixdim*r matrix
        n = len(serverids)
        t = int(math.floor((n - 1)/3))
        if type(secret) is list or type(secret) is tuple:
            r = len(secret)
        else:
            r = 1
            secret = [secret]
        g = pk[0]
        ghat = pk[t+1]
        h = pk[t+3]
        m = secret
        vpolys, wpolys, spolys, vcommits, wcommits, scommits, vhidingpolys, whidingpolys, shidingpolys = ([] for i in range(9))
        pc = PolyCommitPed(t=t, pk=pk, group=group, symflag=False)
        #ZERO and ONE represent 0 and 1 in ZR respectively. This is necessary as field elements have their own datatype.
        ZERO = group.random(ZR) * 0
        ONE = ZERO + 1
        for i in range(matrixdim):
            #Polynomials are represented as a list of coefficients
            #The following statement generates a list of random coefficients
            W_i = list(group.random(ZR, count=t+1))
            #The first entry in the list is the polynomial's constant, followed by the coefficient for x, and so on
            W_i[0] = ONE if i == position else ZERO
            wpolys.append(W_i)
            W_i_hat = list(group.random(ZR, count=t+1))
            wcommits.append(pc.commit(W_i, W_i_hat))
            whidingpolys.append(W_i_hat)
            for j in range(r):
                V_i = list(group.random(ZR, count=t+1))
                V_i[0] = secret[j] if i == position else ZERO
                vpolys.append(V_i)
                V_i_hat = list(group.random(ZR, count=t+1))
                vcommits.append(pc.commit(V_i, V_i_hat))
                vhidingpolys.append(V_i_hat)

                S_i = polynomial_subtract(V_i, polynomial_multiply_constant(W_i, secret[j]))
                spolys.append(S_i)
                S_i_hat = list(group.random(ZR, count=t+1))
                scommits.append(pc.commit(S_i, S_i_hat))
                shidingpolys.append(S_i_hat)
        #Prove that wcommits represents commitments to a standard basis vector
        prove_pc_vector_std_basis(wcommits, wpolys, g, h, group.random(GT), position, group)
        msg = ["payload", vcommits, wcommits, scommits]#, mcommits]
        for pid in serverids:
            send_function(pid, msg)







