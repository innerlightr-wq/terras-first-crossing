"""Exact follow-up research checks for Terras's coefficient stopping-time problem.

Run with Python 3.8+: python3 terras_shortest_failure_probe.py

Input for the seed floor: Rozier and Terracol, Corollary 5.4,
https://arxiv.org/html/2502.00948v5 . All words in the finite example have
length 485 and are already excluded by the published length theorem.

For a shortest-length failure at n, with x_j=T^j(n), the endpoint is strictly
below every interior x_j. Indeed, start at the least interior x_j. Relative
to that start the coefficient is below one by the endpoint since rho_k/rho_j
<1. This shorter first crossing must descend. No subsequent interior state
is below the least interior x_j, so its descent must occur at x_k. All
interior x_j are eligible since x_j>=n>=2 before a first crossing.

Put p=3^d_j, q=2^j, P=3^d, Q=2^k, A=Q-P and C_j/C the respective
prefix/whole-word affine carries. Eliminate n=(C-Q*h)/A from x_j>n+h:
  h < ((p-q)*C + A*C_j)/(p*Q-q*P).
This compares prefixes with a common endpoint. It applies to a shortest
failure only. The modular eligible window taken alone instead bounds
h <= floor((C-2*A)/Q). For n>B replace 2 by B+1.

Small-memory range: write alpha=log_2(3), Lambda=k*ln(2)-d*ln(3)>0.
Failure forces M>n*Lambda>B*Lambda. For adjacent fractions p/q<alpha<P/Q
with P*q-p*Q=1, write (k,d)=u(P,Q)+v(p,q), u=q*k-p*d>=1. If d<Q+q
then v<=0 and Lambda >= P ln(2)-Q ln(3). This code checks both signs
and the decimal bounds with an explicit positive rational log series.
"""
from fractions import Fraction as F


def require(test, message):
    if not test:
        raise ArithmeticError(message)


def log_interval(x, terms=120):
    t = F(x-1, x+1)
    low = 2 * sum(t**(2*i+1)/(2*i+1) for i in range(terms))
    remainder = 2*t**(2*terms+1)/((2*terms+1)*(1-t*t))
    return low, low+remainder


def check_small_memory_range():
    P, Q = 630_118_245_525_664_765, 397_560_349_370_386_783
    p, q = 6_724_555_128_221_608_268, 4_242_721_909_926_539_673
    B = 28*10**18
    ln2, ln3 = log_interval(2), log_interval(3)
    require(P*q-p*Q == 1, 'adjacent fractions')
    lower = P*ln2[0]-Q*ln3[1]
    require(B*lower > F(424722566, 10**8), 'upper approximant gap')
    require(p*ln2[1]-q*ln3[0] < 0, 'lower approximant sign')
    low_next = (P+p)*ln2[0]-(Q+q)*ln3[1]
    high_next = (P+p)*ln2[1]-(Q+q)*ln3[0]
    require(F(75511, 10**5) < B*low_next <= B*high_next < F(75512, 10**5),
            'endpoint of uniform range')
    print(f'If d < {Q+q:,}, a failure requires M > 4.24722566.')


def affine(bits):
    C = d = 0
    for j, bit in enumerate(bits):
        if bit == '1':
            C = 3*C + (1 << j)
            d += 1
    return C, d


def values(bits):
    prefix_odd = 0
    for j, bit in enumerate(bits, start=1):
        prefix_odd += int(bit)
        if j < len(bits):
            require(3**prefix_odd >= 2**j, 'earlier coefficient crossing')
    C, d = affine(bits)
    P, Q = 3**d, 2**len(bits)
    require(P<Q, 'endpoint is a crossing')
    return C, d, P, Q


def perturb(word, sites):
    result = list(word)
    for i in sites:
        require(''.join(result[i:i+4]) == '1011', 'swap site')
        result[i:i+4] = '1101'
    return ''.join(result)


def check_example():
    k, d = 485, 306
    pos = {(3**i).bit_length()-1 for i in range(d)}
    mechanical = ''.join('1' if j in pos else '0' for j in range(k))
    runs = k-d
    low_memory = '1'*(d-runs+1)+'01'*(runs-1)+'0'
    probe = perturb(mechanical, [101, 410])
    for word in (mechanical, low_memory, probe):
        C, odd_count, P, Q = values(word)
        require(len(word)==k and odd_count==d and 1+word.count('01')==runs,
                'same shell and run count')
        print(f'd={d}, k={k}, runs={runs}, M={float(F(C,P)):.12f}')
    A = Q-P
    require(A % 929 == 0, 'factor')
    C, _ = affine(probe)
    j = 84
    Cj, dj = affine(probe[:j])
    p, q = 3**dj, 2**j
    U = F((p-q)*C + A*Cj, p*Q-q*P)
    h0 = C*pow(Q, -1, 929) % 929
    cap = (C-2*A)//Q
    require(cap==73 and h0==66 and 53 < U < F(53453,1000), 'test result')
    print(f'Swapped word: h <= {cap}, h = {h0} mod 929. Shortest-failure prefix gives h < {float(U):.12f}.')


if __name__ == '__main__':
    check_small_memory_range()
    check_example()
