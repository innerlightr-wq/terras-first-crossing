"""Quantify the obstruction, and test the only refinement with geometric leverage."""
import sys; sys.path.insert(0,'scripts')
from fractions import Fraction as F
AL=F(317,200)                     # alpha = log2 3 < 317/200 (certified in the release suite)
def G(s,r=AL):
    t,p=F(0),F(1)
    for _ in range(s): t+=p; p*=r
    return t
print("=== A. why the argument stops at 75 ===")
print("   G_s*L_max must stay below d_min.  L_max is pinned at 70 and d_min is")
print("   essentially constant in s, while G_s is multiplied by alpha each run.")
print("   %4s %13s %13s %9s"%("s","G_s*70","d_min","ratio"))
dmin={70:1.2960e17,75:1.2939e17,76:1.2934e17,80:1.2917e17,85:1.2896e17,91:1.2870e17}
for s in (70,75,76,80,85,91):
    lhs=float(G(s))*70
    print("   %4d %13.4e %13.4e %9.3f %s"%(s,lhs,dmin[s],lhs/dmin[s],"OK" if lhs<dmin[s] else "FAILS"))
m=dmin[75]/(float(G(75))*70)
print("   margin at s=75: %.4f  (only %.1f%% of slack)"%(m,100*(m-1)))
print("   each extra run costs a factor alpha = %.4f, so reaching"%float(AL))
for tgt in (76,80,83,92):
    print("     s=%-3d needs the bound improved by %.4g x"%(tgt,float(AL)**(tgt-75)/m))
print()
print("=== B. the only lever with geometric leverage is the run recurrence ===")
print("   d_min cannot move: it is ~q_next of a convergent of log2 3 with q <~ B ln2/s,")
print("   and the search over ALL convergents q<1e19 and ALL splits found nothing better.")
print("   L_max enters LINEARLY: even L_max=1 would buy only 70x, against the 1.1e3x")
print("   needed for s=92.  So only a smaller growth rate in G_s can work.")
print()
print("=== C. refinement: maximal runs are separated by >=1 even step ===")
print("   Current: X_{i+1} <= alpha X_i      (the even steps between runs are discarded)")
print("   Since runs are MAXIMAL there is at least one even step between consecutive runs,")
print("   and an even step sends y=x+1 to (y+1)/2 <= 0.75y for y>=2, i.e. X -> X-c with")
print("   c = log2(4/3) = 0.415.  So  X_{i+1} <= alpha X_i - c.")
c=F(415,1000)
def Gref(s,L):
    """sum of X_i under X_1=L, X_{i+1} <= alpha X_i - c"""
    tot,X=F(0),F(L)
    for _ in range(s):
        tot+=X; X=AL*X-c
    return tot
for s in (75,76):
    plain=float(G(s))*70; ref=float(Gref(s,70))
    print("   s=%d:  plain sum = %.4e   refined sum = %.4e   gain = %.4f%%"
          %(s,plain,ref,100*(1-ref/plain)))
    print("        refined < d_min = %.4e ?  %s"%(dmin[s], ref<dmin[s]))
print()
print("   The refinement is PROVED but subtractive: it shifts L by c/(alpha-1) = %.3f,"%float(c/(AL-1)))
print("   i.e. ~%.1f%% off a 70-bit budget, against the %.0f%% needed for one more run."
      %(100*float(c/(AL-1))/70, 100*(float(AL)-1)))
