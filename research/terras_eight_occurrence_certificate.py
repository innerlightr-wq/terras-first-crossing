"""An eight-occurrence obstruction on the specified eight-step Terras family.

Run with Python 3.9+ and assertions enabled. No input files or nonstandard
packages are needed. This proves pointwise descent at the first coefficient crossing for every
integer seed >=2 realizing a word in this specified family. Its unique
terminal 00000 gives the distinctness needed; global shortestness is unused.

Conventions
-----------
For a shortcut parity prefix with s odd bits in t positions, rho=3^s/2^t.
Then x_t=rho*(n+M_t), where at an odd position M increases by 1/(3*rho).
The first crossing is at K. Failure implies n<=M_K/(exp(Lambda)-1),
Lambda=K log(2)-D log(3). Each proper prefix has rho>=1, so M_t<=D/3.
For any integer realizer, equality x_i=x_j with i<j<=K-5 would force the
five parity bits at K-5-(j-i) to match those at K-5 by determinism. This
would be an earlier occurrence of 00000. Each middle byte has at most three
zero bits and at most one initial zero. Thus a zero run can have length at
most four across a middle boundary. The explicit head and tail give a single
00000 in the complete word, at K-5. All starts used below are <=K-5.

The words below start with 11111110, end with 111100000, and have S(t) at
every interior eight-step boundary given by boundary_ones(). This certificate
only needs the first-crossing property as a hypothesis; count, positive
height-range, incoming residue classes, and final exclusion are certified here.

If W is a fixed outgoing 120-bit window, equality of parity histories forces
x = r_W mod 2^120. Writing x=r_W+2^120*z makes z integral. Four preceding
bytes have at least 20 odd bits, so 2^32*x = 3^s*y+C forces x = 2^-32*C
mod 3^20. Hence z modulo 3^20 lies in the same translate of our enumerated
set for every occurrence of W. An interval shorter than the modulus cannot
contain two distinct integers with the same residue. The span of the actual
z's is too short to hold eight distinct residues from any translate.

Two independent carry constructions and a cyclic scan check the finite
certificate. The depth-20 proof only needs the same seed at each boundary.
An optional depth-18 check also uses the phase of fifteen five-one bytes.
"""
from fractions import Fraction as F
from itertools import product

K = 1_412_737_287_532_845_639_425
D = 891_337_988_684_307_508_298
END = K-9
N = (END-8)//8
L = 120


def log_bounds(base):
    """Positive atanh series with an upper bound for its tail."""
    u = F(base-1, base+1)
    partial = sum(2*u**(2*j+1)/(2*j+1) for j in range(120))
    return partial, partial+2*u**241/(241*(1-u*u))


def boundary_ones(t):
    return (t*D)//K+2


def alphabet(c):
    out = []
    for letters in product("01", repeat=8):
        w = "".join(letters)
        if w[-1]=="0" and w.count("1")==c and all(
            1000*w[:r].count("1")-631*r >= -800 for r in range(9)
        ):
            out.append(w)
    return out


def byte_carry(w):
    # For a word of length 8 and weight c, 256*x_next=3^c*x+C.
    return sum(2**i * 3**w[i+1:].count("1") for i,b in enumerate(w) if b=="1")


def history_residues():
    mod = 3**20
    choices = {c:alphabet(c) for c in (5,6)}
    assert [len(choices[c]) for c in (5,6)] == [18,7]
    carries = {c:[byte_carry(w) for w in choices[c]] for c in (5,6)}
    total = set()
    for six in (None,0,1,2,3):
        counts = tuple(6 if i==six else 5 for i in range(4))
        # Closed form for the carry of the concatenation.
        factors = [2**(8*i)*3**sum(counts[i+1:]) for i in range(4)]
        scale = pow(2, -(32+L), mod)
        direct = {
            (scale*sum(a*b for a,b in zip(blocks,factors))) % mod
            for blocks in product(*(carries[c] for c in counts))
        }
        # Independent check: propagate an unknown ancestor starting at 0.
        # Its coefficient vanishes modulo 3^20 at the end.
        mixed = {0}
        inv2 = pow(2,-1,mod)
        for c in counts:
            next_set = set()
            # Affine maps are reconstructed from the individual shortcut steps.
            byte_maps = []
            for word in choices[c]:
                p,b=1,0
                for bit in word:
                    if bit=="1":
                        p,b = 3*p,3*b+1
                    p,b = p*inv2%mod,b*inv2%mod
                byte_maps.append((p,b))
            for p,b in byte_maps:
                next_set.update((p*x+b)%mod for x in mixed)
            mixed=next_set
        assert direct == {x*pow(2,-L,mod)%mod for x in mixed}
        total.update(direct)
    assert len(total)==166_104
    return total


def min_span(residues, mod, number):
    ordered=sorted(residues)
    cyclic=ordered + [v+mod for v in ordered[:number-1]]
    return min(cyclic[i+number-1]-v for i,v in enumerate(ordered))


def main():
    assert __debug__, "Run without -O."
    assert END==8+8*N
    phase_step=F(8*D,K)-5
    assert 0<phase_step< F(1,15)
    # Any r consecutive bytes contain 5r plus floor/ceil(r*phase_step) ones;
    # in four consecutive bytes at most one byte has weight six.
    six_count=boundary_ones(END)-boundary_ones(8)-5*N
    assert six_count==8_377_183_976_278_983_657
    # N-14 forward-window starts, remove four for the incoming history.
    # A six can spoil at most 15 potential starts, including near the ends.
    starts_lower=N-18-15*six_count
    assert starts_lower==50_934_401_297_420_950_053
    assert 7*18**15 < starts_lower
    head, tail = "11111110", "111100000"
    words = alphabet(5)+alphabet(6)
    assert all(w.count("0")<=3 and not w.startswith("00") for w in words)
    # Exact junction checks supplement the zero-run argument in the docstring.
    assert all("00000" not in head+w for w in words)
    assert all("00000" not in w+v for w in words for v in words)
    assert all((w+tail).count("00000")==1 and
               (w+tail).find("00000")==len(w+tail)-5 for w in words)
    assert END-120 < K-5

    lo2,hi2=log_bounds(2)
    lo3,hi3=log_bounds(3)
    lam_lo,lam_hi=K*lo2-D*hi3,K*hi2-D*lo3
    assert 0<lam_lo<lam_hi<F(1,300)
    # E=exp(Lambda)-1 > Lambda. In a failure n<=Q_upper.
    memory_upper=F(D,3)
    Q_upper=memory_upper/lam_lo
    # At a boundary t, rho_t=3^(2-{tD/K})*exp(-Lambda*t/K),
    # so 3*(1-Lambda)<rho_t<=9. This uses no assumption on the bytes.
    assert 3*(1-lam_hi)>F("2.99")
    scaled_width=(F("6.01")*Q_upper+9*memory_upper)/2**L
    assert scaled_width<13_234
    residues=history_residues()
    span20=min_span(residues,3**20,8)
    assert span20==14_966
    assert scaled_width < span20 < 3**20

    # For the 15 all-five bytes {tD/K}+15*phase_step < 1.
    # Thus rho_t>3^(1+15*phase_step)*exp(-Lambda)>6.5.
    assert 1+15*phase_step>F(171,100)
    assert lam_hi<F(1,1000)
    assert 3**171*(1-F(1,1000))**100>F(13,2)**100
    scaled_phase_width=(F(5,2)*Q_upper+9*memory_upper)/2**L
    assert scaled_phase_width < 5_506
    residues18={x%3**18 for x in residues}
    span18=min_span(residues18,3**18,8)
    assert len(residues18)==62_964 and span18==5_663
    assert scaled_phase_width<span18<3**18

    print("Eligible starts, conservative lower bound:",starts_lower)
    print("Seven times number of forward words:",7*18**15)
    print("At least eight distinct input states must share a 120-bit window.")
    print("Depth 20: residues",len(residues),"minimum span of eight",span20,
          "available span < 13234")
    print("Depth 18 with phase refinement: residues",len(residues18),
          "minimum span of eight",span18,"available span < 5506")
    print("Unique terminal 00000: every eligible integer realizer descends at K.")


if __name__=="__main__":
    main()
