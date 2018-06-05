from charm.toolbox.pairinggroup import PairingGroup,ZR,G1,G2,GT,pair
from base64 import encodestring, decodestring
import random
import hashlib

#unknown results will occur if denominator does not perfectly divide numerator
#an input of [a,b,c] corresponds to cx^2 + bx + a
def polynomial_divide(numerator, denominator):
    temp = numerator
    factors = []
    while len(temp) >= len(denominator):
        diff = len(temp) - len(denominator)
        factor = temp[len(temp) - 1] / denominator[len(denominator) - 1]
        factors.insert(0, factor)
        for i in range(len(denominator)):
            temp[i+diff] = temp[i+diff] - (factor * denominator[i])
        temp = temp[:len(temp)-1]
    return factors

def polynomial_multiply_constant(poly1, c):
    #myzero will be appropriate whether we are in ZR or G
    #myzero = poly1[0] - poly1[0]
    product = [None] * len(poly1)
    for i in range(len(product)):
        product[i] = poly1[i] * c
    return product

def polynomial_multiply(poly1, poly2):
    #if group == None:
    myzero = poly1[0] - poly1[0]
    #myzero = group.random(ZR)*0
    product = [myzero] * (len(poly1) + len(poly2) -1)
    for i in range(len(poly1)):
        temp = polynomial_multiply_constant(poly2, poly1[i])
        while i > 0:
            temp.insert(0,myzero)
            i -= 1
        product = polynomial_add(product, temp)
    return product

def polynomial_add(poly1, poly2):
    #if group == None:
    #myzero = poly2[0] - poly2[0]
    if len(poly1) >= len(poly2):
        bigger = poly1
        smaller = poly2
    else:
        bigger = poly2
        smaller = poly1
    polysum = [None] * len(bigger)
    for i in range(len(bigger)):
        polysum[i] = bigger[i]
        if i < len(smaller):
            polysum[i] = polysum[i] + smaller[i]
    return polysum

def polynomial_subtract(poly1, poly2):
    negpoly2 = polynomial_multiply_constant(poly2, -1)
    return polynomial_add(poly1, negpoly2)

# Polynomial projection (evaluate the bivariate polynomial at a given y to get a univariate polynomial)
def projf(poly, y):
#Note that ZERO ** 0 evaluates to 0 rather than 1, so this function requires some light tweaking to work.
    ZERO = poly[0][0] * 0
    ONE = poly[0][0] * 0 + 1
    y = ONE * y
    t = len(poly)
    out = [ZERO] * t
    for i in range(t):
        for j in range(t):
            if j != 0:
                out[i] += (poly[i][j]) * (y ** (j))
            else:
                out[i] += (poly[i][j])
    return out

# Polynomial evaluation
def f(poly, x, group=None):
    if type(poly) is not list:
        return "UNDEFINED"
    if group == None:
        ZERO = poly[0] - poly[0]
        ONE = poly[0]/poly[0]
    else:
        ONE = group.random(G1)
        ONE = ONE/ONE
        ZERO = ONE - ONE
    y = ZERO
    xx = ONE
    for coeff in poly:
        y += coeff * xx
        xx *= x
    return y

def f_old(poly, x):
    if type(poly) is not list:
        return "UNDEFINED"
    ZERO = poly[0] - poly[0]
    ONE = poly[0]/poly[0]
    y = ZERO
    xx = ONE
    for coeff in poly:
        y += coeff * xx
        xx *= x
    return y

#interpolates a list of cordinates of the form [x,y] and evaulates at given x
def interpolate_at_x(coords, x, group, order=-1):
    ONE = group.random(ZR)*0 + 1
    if order == -1:
        order = len(coords)
    xs = []
    sortedcoords = sorted(coords, key=lambda x: x[0])
    for coord in sortedcoords:
        xs.append(coord[0])
    S = set(xs[0:order])
    #The following line makes it so this code works for both members of G and ZR
    #out = ONE * (coords[0][1] - coords[0][1])
    out = coords[0][1] - coords[0][1]
    for i in range(order):
        out = out + (lagrange_at_x(S,xs[i],x,group) * sortedcoords[i][1])
    return out

#Turns out a separate interpolation function for commitments is unnecessary
#but I'll leave it here in case it's useful later on
def interpolate_commitment_at_x(coords, x, group, order = -1):
    ONE = coords[0][1]/coords[0][1]
    if order == -1:
        order = len(coords)
    xs = []
    sortedcoords = sorted(coords, key=lambda x: x[0])
    for coord in sortedcoords:
        xs.append(coord[0])
    S = set(xs[0:order])
    out = ONE
    for i in range(order):
        out = out * (sortedcoords[i][1] ** (lagrange_at_x(S,xs[i],x, group)))
    return out

def lagrange_at_x(S,j,x,group):
    ONE = group.random(ZR)*0 + 1
    S = sorted(S)
    assert j in S
    mul = lambda a,b: a*b
    num = reduce(mul, [x - jj  for jj in S if jj != j], ONE)
    den = reduce(mul, [j - jj  for jj in S if jj != j], ONE)
    return num / den

def interpolate_poly(coords, group=None):
    myone = coords[0][1] / coords[0][1]
    myzero = coords[0][1] - coords[0][1]
    if group is not None:
        myone = group.random(ZR)*0 + 1
    #print "IT'SA ME " + str(myzero) + ", THE IDENTITY ELEMENT!"
    #print "Before: " + str(coords[0][1]) + " After: " + str(myzero + coords[0][1])
    poly = [myzero] * len(coords)
    for i in range(len(coords)):
        temp = [myone]
        for j in range(len(coords)):
            if i == j:
                continue
            temp = polynomial_multiply(temp, [ -1 * (coords[j][0] * myone), myone])
            temp = polynomial_divide(temp, [myone * coords[i][0] - myone * coords[j][0]])
        poly = polynomial_add(poly, polynomial_multiply_constant(temp,coords[i][1]))
    return poly

def interpolate_poly_old(coords):
    myone = coords[0][1] / coords[0][1]
    myzero = coords[0][1] - coords[0][1]
    #print "IT'SA ME " + str(myzero) + ", THE IDENTITY ELEMENT!"
    #print "Before: " + str(coords[0][1]) + " After: " + str(myzero + coords[0][1])
    poly = [myzero] * len(coords)
    for i in range(len(coords)):
        temp = [myone]
        for j in range(len(coords)):
            if i == j:
                temp  = polynomial_multiply(temp, [coords[j][1]])
                continue
            temp = polynomial_multiply(temp, [myzero -(coords[j][0] * myone), myone])
            temp = polynomial_divide(temp, [coords[i][0] *myone -coords[j][0] *myone])
        poly = polynomial_add(poly, temp)
    return poly

#this is necessary because ints have a limited size
def hexstring_to_ZR(string, group):
    ZERO = group.random(ZR)*0
    ONE = group.random(ZR)*0 + 1
    i = len(string) - 1
    out = ZERO
    while i >= 0:
        if i > 0:
            temp = ONE * 2
            temp = temp ** (i*4)
        else:
            temp = ONE
        out = out + temp*int(string[i],16)
        i = i - 1
    return out

def intstring_to_ZR(string, group):
    i = len(string) - 1
    ZERO = group.random(ZR)*0
    ONE = group.random(ZR)*0 + 1
    out = ZERO
    while i >= 0:
        if i > 0:
            temp = ONE * 10
            temp = temp ** (i)
        else:
            temp = ONE
        out = out + temp*int(string[i])
        i = i - 1
    return out

#Check that a subset of t+1 points will correctly interpolate to the polynomial which contains all points
def check_commitment_integrity(commitments, t, group):
    points = []
    for i in commitments:
        points.append([i, commitments[i]])
    out = True
    for i in range(t+1,len(commitments)):
        out = out and (interpolate_at_x(points[:t+1], points[i][0], group) == (points[i][1]))
    return out

##### Hereafter be ZK Proofs ############

### NIZK that two group elements are raised to the same exponent ###
#Let g1=e1^exp and g2=e2^exp. Generate a ZK proof that 
#g1 and g2 are the result of raising e1 and e2 to the same exponent
def prove_same_exponent(e1,e2,exp,group):
    blind = group.random(ZR)
    k1 = e1 ** blind
    k2 = e2 ** blind
    challenge = hexstring_to_ZR(hashlib.sha256(str(k1) + str(k2)).hexdigest(), group)
    s = blind - challenge*exp
    return [k1, k2, s]

#"proof" is a list of the form [k1, k2, challenge, s]
def check_same_exponent_proof(proof, e1, e2, g1, g2, group):
    challenge = hexstring_to_ZR(hashlib.sha256(str(proof[0]) + str(proof[1])).hexdigest(), group)
    eq1 = str(proof[0]) == str(e1 ** proof[2] * g1 ** challenge)
    eq2 = str(proof[1]) == str(e2 ** proof[2] * g2 ** challenge)
    return eq1 and eq2

### NIZK that a vector of polynomial commitments correspond to a standard basis vector ###
def prove_pc_vector_std_basis(pcvector, pcpolys, g, h, gt, j, group):
    #This proof is formalized in Appendix B of Henry's CCS11 paper
    assert len(pcvector) == len(pcpolys)
    r = len(pcvector)
    #Start by making the first part noninteractive by using Fiat-Shamir for the challenge vector (a)
    a = []
    f_a = [0]
    for i in range(r):
        a_i = hexstring_to_ZR(hashlib.sha256(str(pcvector[i])).hexdigest(), group)
        f_a = polynomial_add(f_a, polynomial_multiply_constant(pcpolys[i], a_i))
        a.append(a_i)
    #evaluate f_a at zero
    y = f(f_a, 0)
    blind = group.random(ZR)
    Y = gt ** (y * blind)
    #TODO: Insert ZKPOK that f_a(0) != 0
    v = h ** blind
    hlist = []
    for i in range(r):
        hlist.append(v**a[i])
    ## Next line should be gt, not g
    prove_one_of_r_same_exponent(g, h, [g ** (y * blind)] * r, hlist, j, y * blind, group)

def check_pc_vector_std_basis_proof(pcvector, proof):
    return True

### NIZK that 1 of many pairs of group elements share the same exponent ###
### let g_1 = g^x_1, g_2 = g^x_2, h_1 = h^y_1, etc. Prove that x_j = y_j for at least one j
def prove_one_of_r_same_exponent(g, h, glist, hlist, j, exp, group, k=40):
    assert len(glist) == len(hlist)
    r = len(glist)
    assert j < r and j > 0
    modulus = 2 ** 40
    blinds = []
    cs = []
    g_icommits = []
    h_icommits = []
    #challenge = int(hexstring_to_ZR(hashlib.sha256(str(glist) + str(hlist)).hexdigest(), group)) %modulus
    #This line pains me
    challenge = int(str(hexstring_to_ZR(hashlib.sha256(str(glist) + str(hlist)).hexdigest(), group))[:10]) %modulus 
    c_j = challenge
    for i in range(r):
        blind_i = group.random(ZR)
        #c_i = group.random(ZR)%modulus
        c_i = random.randint(0,modulus-1)
        g_icommit = g**blind_i if i == j else g**blind_i * glist[i] ** c_i
        h_icommit = h**blind_i if i == j else h**blind_i * hlist[i] ** c_i
        c_j = c_j if i == j else (c_j - c_i)%modulus
        blinds.append(blind_i)
        cs.append(c_i)
        g_icommits.append(g_icommit)
        h_icommits.append(h_icommit)
    cs[j] = c_j
    
    #Verification part. This should be elsewhere
    bs = []
    v = 0
    for i in range(r):
        b_i = random.randint(0,modulus-1)
        v = v + b_i*blinds[i]
        bs.append(b_i)
    p1lhs = 1
    p1rhs = g ** v
    p2lhs = 1
    p2rhs = h ** v
    p3rhs = 0
    for i in range(r):
        p1lhs = p1lhs * g_icommits[i] ** bs[i]
        p1rhs = p1rhs * glist[i] ** ((bs[i] * cs[i])%modulus)
        #p1rhs = p1rhs * glist[i] ** ((bs[i] * cs[i]))
        p2lhs = p2lhs * h_icommits[i] ** bs[i]
        p2rhs = p2rhs * hlist[i] ** ((bs[i] * cs[i])%modulus)
        #p2rhs = p2rhs * hlist[i] ** ((bs[i] * cs[i]))
        p3rhs = (p3rhs + cs[i])%modulus
    print p1lhs == p1rhs
    print p2lhs == p2rhs 
    print challenge == p3rhs
    return [g_icommits, h_icommits]





