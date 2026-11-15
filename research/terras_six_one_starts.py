"""Candidate extension of the eight-occurrence family theorem for independent audit.
Run python3 terras_six_one_starts.py (Python 3 + NumPy).
Uses integer residue arithmetic and Fraction bounds, without expanding the word.
"""
from fractions import Fraction as F
from itertools import product
import numpy as np

K = 1_412_737_287_532_845_639_425
D = 891_337_988_684_307_508_298
N = (K-17)//8
HEAD, TAIL = "11111110", "111100000"
L = 112


def log_bounds(base):
    x = F(base-1, base+1)
    lo = sum(2*x**(2*j+1)/(2*j+1) for j in range(120))
    return lo, lo+2*x**241/(241*(1-x*x))


def alphabet(c, size):
    return ["".join(w)+"0" for w in product("01",repeat=7) if w.count("1")==c
            and (size==21 or w[:2]!=("0","0"))]


def carry(w):
    return sum(2**j*3**w[j+1:].count("1") for j in range(len(w)) if w[j]=="1")


def span_certificate(letters, depth, forced):
    # Last depth bytes all have exactly five ones. In 2^(8*depth)*x
    # the incoming ancestor is multiplied by 3^(5*depth) and vanishes
    # modulo mod. A common forward window then subtracts a common dyadic
    # residue before dividing by 2^112. This only translates the set.
    mod=3**(5*depth)
    assert 2*mod < 2**63
    residues=np.array([0], dtype=np.int64)
    carries=[carry(w) for w in letters]
    for j in range(depth):
        scale = pow(2,8*j-8*depth-L,mod)*pow(3,5*(depth-1-j),mod)%mod
        offsets=np.array([scale*c%mod for c in carries], dtype=np.int64)
        residues=np.unique((residues[:,None]+offsets[None,:])%mod)
    cyclic=np.concatenate((residues,residues[:forced-1]+mod))
    return len(residues), int(np.min(cyclic[forced-1:]-residues))


def check(size, past, marker, expected_count, expected_span):
    lo2,hi2=log_bounds(2)
    lo3,hi3=log_bounds(3)
    lam_lo,lam_hi=K*lo2-D*hi3,K*hi2-D*lo3
    assert 0<lam_lo<lam_hi<F(1,1000)
    beta_lo, beta_hi = lo2/hi3, hi2/lo3
    eta=F(8*D,K)-5
    assert 0<eta<F(1,20) and 8+8*N==K-9
    s=lambda t:t*D//K+2
    six=s(K-9)-s(8)-5*N
    assert s(8)==7 and s(K-9)==D-4 and six==8_377_183_976_278_983_657
    A5,A6=alphabet(5,size),alphabet(6,size)
    assert len(A5)==size and len(A6)==7
    words=A5+A6
    # Any h+14 adjacent middle bytes have <=1 six since (h+14)*eta<1.
    # A window starts with a six, the following thirteen and the prior
    # h bytes are fives. Delete up to h+13 sixes for edge starts.
    assert (past+14)*eta<1
    available=six-past-13
    labels=7*size**13
    forced=(available+labels-1)//labels
    assert forced==({20:15,21:8}[size])
    # At a boundary t, rho=3^(2-{tD/K})*exp(-Lambda*t/K).
    # A six next implies {tD/K}>=1-eta.
    assert 3*(1-lam_hi)>F("2.99") and 3**21<F("3.17")**20
    # All middle boundary coefficients exceed 2.99, including those
    # before a five. Each middle byte adds carry(w)/(3^c*rho_start)
    # to memory, <=c/3 even if it crosses inside the byte.
    assert all(F(carry(w),3**w.count("1"))/F("2.99") < F(w.count("1"),3)
               for w in words)
    # Proper head/tail inputs have rho>=1. The tail lower bound uses
    # S(K)-beta*K=-Lambda/log(3), then subtracts the remaining ones.
    for j in range(1,len(HEAD)+1):
        assert HEAD[:j].count("1")-j*beta_hi > 0
    for j in range(len(TAIL)):
        assert (len(TAIL)-j)*beta_lo-TAIL[j:].count("1")-lam_hi/lo3>0
    if size==20:
        # At a boundary S(t)-beta*t>1-Lambda/log(3)>.999.
        assert lam_hi/lo3<F(".001")
        assert all(F(".999")+w[:j].count("1")-j*beta_hi > 0
                   for w in words for j in range(9))
    # A marker of length <=10 in the interior fits in three consecutive
    # bytes. Head and tail junctions must be checked as well.
    assert all(marker not in a+b+c for a,b,c in product(words,repeat=3))
    assert all(marker not in HEAD+a+b for a,b in product(words,repeat=2))
    for a,b in product(words,repeat=2):
        ending=a+b+TAIL
        assert ending.find(marker)==len(ending)-len(marker)
    assert K-9-L <= K-len(marker)
    # Determinism: equality of two selected states repeats the unique
    # terminal marker earlier. Thus integer states are distinct.
    M=F(D,3)
    Q=M/lam_lo
    # If x_K>=n, n<=Q, and all chosen states share the same n and
    # lie in (2.99 n, 3.17(n+M)). Divide their width by 2^112.
    width=(F(".18")*Q+F("3.17")*M)/2**L
    count,span=span_certificate(A5,past,forced)
    assert count==expected_count and span==expected_span
    assert width<101468<span<3**(5*past)
    print(f"{size}/7: >= {available} starts, {labels} labels, >= {forced} repeats")
    print(f"  {count} residues, span {span}, available width < 101468")


if __name__=="__main__":
    check(20,4,"00000",128_000,107_996)
    check(21,5,"0111100000",3_111_696,152_838)
