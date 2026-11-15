"""Verify the 21/7 prefix comparison, using the span checked by span.cpp."""
from fractions import Fraction as F
from itertools import product
from math import gcd

K = 1_412_737_287_532_845_639_425
D = 891_337_988_684_307_508_298
HEAD, TAIL = "11111110", "111100000"
MARKER = "0111100000"
SPAN = 3_611_678_508

def log_interval(base):
    x = F(base - 1, base + 1)
    lower = sum(2 * x**i / i for i in range(1, 240, 2))
    return lower, lower + 2 * x**241 / (241 * (1 - x*x))

def count(t):
    return t*D//K + 2

def multiplier(bits):
    return F(3**bits.count("1"), 2**len(bits))

def memory_cost(bits):
    return sum(F(2**i, 3**(bits[:i].count("1") + 1))
               for i, bit in enumerate(bits) if bit == "1")

lo2, hi2 = log_interval(2)
lo3, hi3 = log_interval(3)
lam_lo, lam_hi = K*lo2-D*hi3, K*hi2-D*lo3
eta = F(8*D, K)-5
assert 0 < lam_lo < lam_hi < F(1, 10**18)
assert 0 < 20*eta < 1
nbytes = (K-17)//8
end = K-len(TAIL)
assert 8+8*nbytes == end and count(8)==7 and count(end)==D-4
sixes = D-11-5*nbytes
labels = 7*21**12
assert 162*labels < sixes-19 <= 163*labels

alphabet = ["".join(w)+"0" for w in product("01", repeat=7) if w.count("1") in (5,6)]
regular = [w for w in alphabet if w.count("1")==5 and w!="00111110"]
six = [w for w in alphabet if w.count("1")==6]
bad = "00111110"
assert len(regular) == 20 and len(six) == 7
assert max(map(memory_cost,regular)) == F(682,243)
assert max(map(memory_cost,six)) == F(1330,729)
assert memory_cost(bad) == F(844,243)
assert all(multiplier(w[:j])*F("2.99")>1 for w in regular+six for j in range(1,9))
assert all(multiplier(bad[:j])*F("2.99")>1 for j in (1,3,4,5,6,7,8))
assert multiplier(bad[:2]) == F(1,4)
assert all(multiplier(HEAD[:j])>1 for j in range(1,9))
for j in range(len(TAIL)):
    remaining=TAIL[j:]
    assert len(remaining)*lo2 - remaining.count("1")*hi3 > lam_hi
for a,b,c in product(alphabet,repeat=3):
    assert MARKER not in a+b+c
for a,b in product(alphabet,repeat=2):
    assert MARKER not in HEAD+a+b
    ending = a+b+TAIL
    assert ending.find(MARKER)==len(ending)-len(MARKER)
assert memory_cost(HEAD) < F(7,3)
assert 3*(1-lam_hi) > F("2.999")
assert F(256,81)-F("2.999") < F(".162")
for w in alphabet:
    assert memory_cost(w)/F("2.999") < F(w.count("1"),3)

# On an earlier first-crossing prefix any completed bad byte starts at rho>=4.
assert max(map(memory_cost,regular))/F("2.999") < F(".94")
assert max(map(memory_cost,six))/F("2.999") < F(".94")
assert memory_cost(bad)/4 < F(".94")

floor = (2**104 * SPAN - F(256,81)*F(D,3)) / F(".162")
assert floor > 452 * 10**39
assert floor*lam_lo > F(D,3)  # Also settles the full endpoint.
# If a middle crossing fails, 0<E=j ln2-a ln3<M_j/n<=D/(3*floor).
# h=K*a-D*j therefore lies strictly between -K*D/(3*floor*ln3)
# and K*Lambda/ln3. Use exact bounds to cover all such integers.
assert -K*D/(3*floor*lo3) > -1
assert K*lam_hi/lo3 < 131
assert gcd(K,D)==1
inv=pow(D,-1,K)
possible=[]
for h in range(131):
    j=(-h*inv)%K
    t=j-2
    if not (8<=t<end and t%8==0):
        continue
    a=count(t)
    if K*a-D*j!=h or count(t+8)!=a+5:
        continue
    elo, ehi = j*lo2-a*hi3, j*hi2-a*lo3
    if ehi < 0:
        continue
    assert elo > 0
    if floor*elo > F(a,3):
        continue
    possible.append((h,j))
    cap = F(7,3) + F(".94")*((t-8)//8)
    assert cap/(elo*floor) < F(".952")
    print("Conservative upper bound on failure seed:", cap//elo + 1)
assert possible == [(60, 651_838_711_708_843_936_730)]
print("Conservative lower bound on every full-word seed:", int(floor))
print("Unique site after coarse memory test:", possible[0])
print("All possible first crossings covered on this shell.")
