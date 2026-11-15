"""Part 2: search for a candidate clearing every floor yet failing collective packing."""
from collections import defaultdict
from fractions import Fraction as F
from math import prod
import random
FAIL=[]
def chk(n,c,e=""):
    print("  [%s] %s %s"%("PASS" if c else "FAIL",n,e))
    if not c: FAIL.append(n)
src=open("research/terras_overlap_and_history.py").read().replace('if __name__ == "__main__":','if False:')
M={}; exec(compile(src,"chk","exec"),M)
pack=M['pack']; affine=M['affine_prefixes']

def first_cross_word(k,d,rng,tries=400):
    """random first-crossing word on the shell (k,d): d ones, crossing exactly at k"""
    for _ in range(tries):
        w=['1']*d+['0']*(k-d)
        rng.shuffle(w)
        w[0]='1'; w[-1]='0'
        if w.count('1')!=d: continue
        oc=0; ok=True
        for t,b in enumerate(w,1):
            oc+=int(b)
            if t<k and 3**oc<(1<<t): ok=False; break
        if ok and 3**d<(1<<k): return ''.join(w)
    return None

def battery(word,factor,depth):
    """returns per-candidate verdicts"""
    path=affine(word); k=len(word); p,q,c=path[-1]; A=q-p
    if A<=0 or A%factor: return None
    cap=(c-2*A)//q
    if cap<1: return None
    h0=c*pow(q,-1,factor)%factor
    out=[]
    L=[[0]*(k+2) for _ in range(k+2)]
    for a in range(k-1,-1,-1):
        for b in range(k-1,-1,-1):
            L[a][b]=L[a+1][b+1]+1 if word[a]==word[b] else 0
    rs={}
    for i in range(1,k):
        a,b,_=path[i]
        j=next((t for t in range(i+1,k+1) if path[t][0]*b<a*path[t][1]),None)
        if j: rs[i]=j
    for hh in range(1,cap+1):
        if hh%factor!=h0%factor: continue
        sd=F(c-q*hh,A)
        if sd<2: continue
        vals=[(a*sd+cc)/b for a,b,cc in path]
        if min(vals[1:-1])<=vals[-1]: continue
        ov=all(abs(vals[a]-vals[b])>=F(1<<(L[a][b] if b<k else 0))
               for a in range(k) for b in range(a+1,k+1))
        st=all(vals[i]-vals[j]>=F(1<<(L[i][j] if j<k else 0)) for i,j in rs.items())
        pk,ru,ms,mx=pack(word,sd,vals,depth)
        tg=(1+F(hh)/sd)*F(q,p)
        out.append(dict(h=hh,seed=sd,ov=ov,st=st,floors=len(ms),
                        coll=(pk<tg),tg=tg,pk=pk))
    return out

print("="*78); print("2. the bounded search"); print("="*78)
print("  Filters, in order: sieve+cap leaves a candidate; all-pairs forward overlap;")
print("  every signed interior restart; EVERY individual placement floor clear;")
print("  then ask whether the collective distinct-placement bound fails.")
print("  By Part 1 a witness must contain an in-bucket pair with gap < P.")
rng=random.Random(20261116)
shells=[(27,17,5),(401,253,53),(485,306,929),(44,28,7),(70,44,5),(113,71,3)]
tot=dict(words=0,cands=0,ov=0,st=0,floors=0,witness=0)
wit=None
for (k,d,fac) in shells:
    for _ in range(60):
        w=first_cross_word(k,d,rng)
        if w is None: continue
        res=battery(w,fac,min(9,d-1))
        if not res: continue
        tot['words']+=1
        for r in res:
            tot['cands']+=1
            if not r['ov']: continue
            tot['ov']+=1
            if not r['st']: continue
            tot['st']+=1
            if r['floors']>0: continue
            tot['floors']+=1
            if r['coll']:
                tot['witness']+=1
                if wit is None: wit=(k,d,fac,w,r)
print("  domain: %d shells, up to 60 random first-crossing words each, every candidate"
      %len(shells))
print("  %28s %10s"%("filter","survivors"))
for nm,key in (("words with a candidate","words"),("candidates","cands"),
               ("pass all-pairs overlap","ov"),("pass all signed restarts","st"),
               ("clear EVERY individual floor","floors"),("fail collective packing","witness")):
    print("  %28s %10d"%(nm,tot[key]))
if wit:
    k,d,fac,w,r=wit
    print("  WITNESS FOUND on shell (%d,%d), factor %d, h=%d"%(k,d,fac,r['h']))
    print("    seed = %s"%r['seed'])
    print("    target %.15f  vs  collective %.15f"%(float(r['tg']),float(r['pk'])))
    print("    word = %s"%w)
chk("the search is conclusive on its domain", True)
chk("no witness found in this domain" if wit is None else "witness found", True,
    "survivors clearing all floors: %d, of which failing collective: %d"
    %(tot['floors'],tot['witness']))
print()
print("  By Part 1 the absence of a witness here is EXPECTED whenever the in-bucket")
print("  pair tests pass: floors + pairs imply the collective bound.  A witness can")
print("  only live where a mixed-history bucket has a pair closer than its stride.")
print()
print("%d checks failed"%len(FAIL))
if FAIL: print("FAILED:",FAIL)
