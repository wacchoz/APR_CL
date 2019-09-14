# -*- coding: utf-8 -*-

import copy
import time
from math import gcd  # version >= 3.5



# primality test by trial division
def isprime_slow(n):
    if n<2:
        return False
    elif n==2 or n==3:
        return True
    elif n%2==0:
        return False
    else:
        i = 3
        while i*i <= n:
            if n%i == 0:
                return False
            i+=2
    return True        


# v_q(t): how many time is t divided by q
def v(q, t):
    ans = 0
    while(t % q == 0):
        ans +=1
        t//= q
    return ans


def prime_factorize(n):
    ret = []
    p=2
    while p*p <= n:
        if n%p==0:
            num = 0
            while n%p==0:
                num+=1
                n//= p
            ret.append((p,num))
        p+= 1

    if n!=1:
        ret.append((n,1))

    return ret


# calculate e(t)
def e(t):
    s = 1
    q_list = []
    for q in range(2, t+2):
        if t % (q-1) == 0 and isprime_slow(q):
            s *= q ** (1+v(q,t))
            q_list.append(q)
    return 2*s, q_list

# Jacobi sum
class JacobiSum(object):
    def __init__(self, p, k, q):
        self.p = p
        self.k = k
        self.q = q
        self.m = (p-1)*p**(k-1)
        self.pk = p**k
        self.coef = [0]*self.m

    # 1
    def one(self):
        self.coef[0] = 1
        for i in range(1,self.m):
            self.coef[i] = 0
        return self


    # product of JacobiSum
    # jac : JacobiSum
    def mul(self, jac):
        m = self.m
        pk = self.pk
        j_ret=JacobiSum(self.p, self.k, self.q)
        for i in range(m):
            for j in range(m):
                if (i+j)% pk < m:
                    j_ret.coef[(i+j)% pk] += self.coef[i] * jac.coef[j]
                else:
                    r = (i+j) % pk - self.p ** (self.k-1)                    
                    while r>=0:
                        j_ret.coef[r] -= self.coef[i] * jac.coef[j]
                        r-= self.p ** (self.k-1)

        return j_ret


    def __mul__(self, right):
        if type(right) is int:
            # product with integer
            j_ret=JacobiSum(self.p, self.k, self.q)
            for i in range(self.m):
                j_ret.coef[i] = self.coef[i] * right
            return j_ret
        else:
            # product with JacobiSum
            return self.mul(right)
    

    # power of JacobiSum（x-th power mod n）
    def modpow(self, x, n):
        j_ret=JacobiSum(self.p, self.k, self.q)
        j_ret.coef[0]=1
        j_a = copy.deepcopy(self)
        while x>0:
            if x%2==1:
                j_ret = (j_ret * j_a).mod(n)
            j_a = j_a*j_a
            j_a.mod(n)
            x //= 2
        return j_ret
    

    # applying "mod n" to coefficient of self
    def mod(self, n):
        for i in range(self.m):
            self.coef[i] %= n
        return self
    

    # operate sigma_x
    # verification for sigma_inv
    def sigma(self, x):
        m = self.m
        pk = self.pk
        j_ret=JacobiSum(self.p, self.k, self.q)
        for i in range(m):
            if (i*x) % pk < m:
                j_ret.coef[(i*x) % pk] += self.coef[i]
            else:
                r = (i*x) % pk - self.p ** (self.k-1)                    
                while r>=0:
                    j_ret.coef[r] -= self.coef[i]
                    r-= self.p ** (self.k-1)
        return j_ret
    
                
    # operate sigma_x^(-1)
    def sigma_inv(self, x):
        m = self.m
        pk = self.pk
        j_ret=JacobiSum(self.p, self.k, self.q)
        for i in range(pk):
            if i<m:
                if (i*x)%pk < m:
                    j_ret.coef[i] += self.coef[(i*x)%pk]
            else:
                r = i - self.p ** (self.k-1)
                while r>=0:
                    if (i*x)%pk < m:
                        j_ret.coef[r] -= self.coef[(i*x)%pk]
                    r-= self.p ** (self.k-1)

        return j_ret
    

    # Is self p^k-th root of unity (mod N)
    # if so, return h where self is zeta^h
    def is_root_of_unity(self, N):
        m = self.m
        p = self.p
        k = self.k

        # case of zeta^h (h<m)
        one = 0
        for i in range(m):
            if self.coef[i]==1:
                one += 1
                h = i
            elif self.coef[i] == 0:
                continue
            elif (self.coef[i] - (-1)) %N != 0:
                return False, None
        if one == 1:
            return True, h

        # case of zeta^h (h>=m)
        for i in range(m):
            if self.coef[i]!=0:
                break
        r = i % (p**(k-1))
        for i in range(m):
            if i % (p**(k-1)) == r:
                if (self.coef[i] - (-1))%N != 0:
                    return False, None
            else:
                if self.coef[i] !=0:
                    return False, None

        return True, (p-1)*p**(k-1)+ r
            

# find primitive root
def smallest_primitive_root(q):
    for r in range(2, q):
        s = set({})
        m = 1
        for i in range(1, q):
            m = (m*r) % q
            s.add(m)
        if len(s) == q-1:
            return r
    return None   # error


# calculate f_q(x)
def calc_f(q):
    g = smallest_primitive_root(q)
    m = {}
    for x in range(1,q-1):
        m[pow(g,x,q)] = x
    f = {}
    for x in range(1,q-1):
        f[x] = m[ (1-pow(g,x,q))%q ]

    return f


# sum zeta^(a*x+b*f(x))
def calc_J_ab(p, k, q, a, b):
    j_ret = JacobiSum(p,k,q)
    f = calc_f(q)
    for x in range(1,q-1):
        pk = p**k
        if (a*x+b*f[x]) % pk < j_ret.m:
            j_ret.coef[(a*x+b*f[x]) % pk] += 1
        else:
            r = (a*x+b*f[x]) % pk - p**(k-1)
            while r>=0:
                j_ret.coef[r] -= 1
                r-= p**(k-1)
    return j_ret


# calculate J(p,q)（p>=3 or p,q=2,2）
def calc_J(p, k, q):
    return calc_J_ab(p, k, q, 1, 1)

            
# calculate J_3(q)（p=2 and k>=3）
def calc_J3(p, k, q):
    j2q = calc_J(p, k, q)
    j21 = calc_J_ab(p, k, q, 2, 1)
    j_ret = j2q * j21
    return j_ret


# calculate J_2(q)（p=2 and k>=3）
def calc_J2(p, k, q):
    j31 = calc_J_ab(2, 3, q, 3, 1)
    j_conv = JacobiSum(p, k, q)
    for i in range(j31.m):
        j_conv.coef[i*(p**k)//8] = j31.coef[i]
    j_ret = j_conv * j_conv
    return j_ret


# in case of p>=3
def APRtest_step4a(p, k, q, N):

    print("Step 4a. (p^k, q = {0}^{1}, {2})".format(p,k,q))
    
    J = calc_J(p, k, q)
    # initialize s1=1
    s1 = JacobiSum(p,k,q).one()
    # J^Theta
    for x in range(p**k):
        if x % p == 0:
            continue
        t = J.sigma_inv(x)
        t = t.modpow(x, N)
        s1 = s1 * t
        s1.mod(N)

    # r = N mod p^k
    r = N % (p**k)

    # s2 = s1 ^ (N/p^k)
    s2 = s1.modpow(N//(p**k), N)

    # J^alpha
    J_alpha = JacobiSum(p,k,q).one()
    for x in range(p**k):
        if x % p == 0:
            continue
        t = J.sigma_inv(x)
        t = t.modpow((r*x)//(p**k), N)
        J_alpha = J_alpha * t
        J_alpha.mod(N)

    # S = s2 * J_alpha
    S = (s2 * J_alpha).mod(N)

    # Is S root of unity
    exist, h = S.is_root_of_unity(N)

    if not exist:
        # composite!
        return False, None
    else:
        # possible prime
        if h%p!=0:
            l_p = 1
        else:
            l_p = 0
        return True, l_p


# in case of p=2 and k>=3
def APRtest_step4b(p, k, q, N):

    print("Step 4b. (p^k, q = {0}^{1}, {2})".format(p,k,q))

    J = calc_J3(p, k, q)
    # initialize s1=1
    s1 = JacobiSum(p,k,q).one()
    # J3^Theta
    for x in range(p**k):
        if x % 8 not in [1,3]:
            continue
        t = J.sigma_inv(x)
        t = t.modpow(x, N)
        s1 = s1 * t
        s1.mod(N)

    # r = N mod p^k
    r = N % (p**k)

    # s2 = s1 ^ (N/p^k)
    s2 = s1.modpow(N//(p**k), N)

    # J3^alpha
    J_alpha = JacobiSum(p,k,q).one()
    for x in range(p**k):
        if x % 8 not in [1,3]:
            continue
        t = J.sigma_inv(x)
        t = t.modpow((r*x)//(p**k), N)
        J_alpha = J_alpha * t
        J_alpha.mod(N)

    # S = s2 * J_alpha * J2^delta
    if N%8 in [1,3]:
        S = (s2 * J_alpha ).mod(N)
    else:
        J2_delta = calc_J2(p,k,q)
        S = (s2 * J_alpha * J2_delta).mod(N)

    # Is S root of unity
    exist, h = S.is_root_of_unity(N)

    if not exist:
        # composite 
        return False, None
    else:
        # possible prime
        if h%p!=0 and (pow(q,(N-1)//2,N) + 1)%N==0:
            l_p = 1
        else:
            l_p = 0
        return True, l_p


# in case of p=2 and k=2
def APRtest_step4c(p, k, q, N):

    print("Step 4c. (p^k, q = {0}^{1}, {2})".format(p,k,q))

    J2q = calc_J(p, k, q)

    # s1 = J(2,q)^2 * q (mod N)
    s1 = (J2q * J2q * q).mod(N)

    # s2 = s1 ^ (N/4)
    s2 = s1.modpow(N//4, N)

    if N%4 == 1:
        S = s2
    elif N%4 == 3:
        S = (s2 * J2q * J2q).mod(N)
    else:
        print("Error")

    # Is S root of unity
    exist, h = S.is_root_of_unity(N)

    if not exist:
        # composite
        return False, None
    else:
        # possible prime
        if h%p!=0 and (pow(q,(N-1)//2,N) + 1)%N==0:
            l_p = 1
        else:
            l_p = 0
        return True, l_p


# in case of p=2 and k=1
def APRtest_step4d(p, k, q, N):

    print("Step 4d. (p^k, q = {0}^{1}, {2})".format(p,k,q))

    S2q = pow(-q, (N-1)//2, N)
    if (S2q-1)%N != 0 and (S2q+1)%N != 0:
        # composite
        return False, None
    else:
        # possible prime
        if (S2q + 1)%N == 0 and (N-1)%4==0:
            l_p=1
        else:
            l_p=0
        return True, l_p


# Step 4
def APRtest_step4(p, k, q, N):

    if p>=3:
        result, l_p = APRtest_step4a(p, k, q, N)
    elif p==2 and k>=3:
        result, l_p = APRtest_step4b(p, k, q, N)
    elif p==2 and k==2:
        result, l_p = APRtest_step4c(p, k, q, N)
    elif p==2 and k==1:
        result, l_p = APRtest_step4d(p, k, q, N)
    else:
        print("error")

    if not result:
        print("Composite")

    return result, l_p


def APRtest(N):
    t_list = [
        2,
        12,
        60,
        180,
        840,
        1260,
        1680,
        2520,
        5040,
        15120,
        55440,
        110880,
        720720,
        1441440,
        4324320,
        24504480,
        73513440
        ]

    print("N=", N)

    if N<=3:
        print("input should be greater than 3")
        return False
 
    # Select t
    for t in t_list:
        et, q_list = e(t)
        if N < et*et:
            break
    else:
        print("t not found")
        return False

    print("t=", t)
    print("e(t)=", et, q_list)

    # Step 1
    print("=== Step 1 ===")
    g = gcd(t*et, N)
    if g > 1:
        print("Composite")
        return False

    # Step 2
    print("=== Step 2 ===")
    l = {}
    fac_t = prime_factorize(t)
    for p, k in fac_t:
        if p>=3 and pow(N, p-1, p*p)!=1:
            l[p] = 1
        else:
            l[p] = 0
    print("l_p=", l)

    # Step 3 & Step 4
    print("=== Step 3&4 ===")
    for q in q_list:
        if q == 2:
            continue
        fac = prime_factorize(q-1)
        for p,k in fac:

            # Step 4
            result, l_p = APRtest_step4(p, k, q, N)

            if not result:
                # composite
                print("Composite")
                return False
            elif l_p==1:
                l[p] = 1

    # Step 5          
    print("=== Step 5 ===")
    print("l_p=", l)
    for p, value in l.items():
        if value==0:
            # try other pair of (p,q)
            print("Try other (p,q). p={}".format(p))
            count = 0
            i = 1
            found = False
            # try maximum 30 times
            while count < 30:
                q = p*i+1
                if N%q != 0 and isprime_slow(q) and (q not in q_list):
                    count += 1

                    k = v(p, q-1)
                    # Step 4
                    result, l_p = APRtest_step4(p, k, q, N)

                    if not result:
                        # composite
                        print("Composite")
                        return False
                    elif l_p == 1:
                        found = True
                        break
                i += 1
            if not found:
                print("error in Step 5")
                return False

    # Step 6
    print("=== Step 6 ===")
    r = 1
    for t in range(t-1):
        r = (r*N) % et
        if r!=1 and r!= N and N % r == 0:
            print("Composite", r)
            return False
    # prime!!
    print("Prime!")
    return True


if __name__ == '__main__':


    start_time = time.time()

    APRtest(2**521-1)   # 157 digit, 18 sec
#    APRtest(2**1279-1)  # 386 digit, 355 sec
#    APRtest(2074722246773485207821695222107608587480996474721117292752992589912196684750549658310084416732550077)

    end_time = time.time()
    print(end_time - start_time, "sec")
