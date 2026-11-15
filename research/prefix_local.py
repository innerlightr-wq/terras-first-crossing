"""Phase 3: can the seed floor be localized before an earlier crossing?"""
from fractions import Fraction as F
from decimal import Decimal, getcontext
import math
getcontext().prec=260
AL=Decimal(3).ln()/Decimal(2).ln()
LN2=F(int(Decimal(2).ln()*Decimal(10)**220),10**220)
LN3=F(int(Decimal(3).ln()*Decimal(10)**220),10**220)
K=1412737287532845639425; D=891337988684307508298
SIGMA163=3611678508; L=104; PAST=7; FW=12; labels=7*21**12
j=651838711708843936730; t0=j-2
def a_of(t): return t*D//K+2
a=a_of(t0)
lam=F(K)*LN2-F(D)*LN3
print("="*78); print("PHASE 3.  prefix-localization at the surviving site"); print("="*78)
print("  The 21/7 theorem quantifies over seeds realizing a COMPLETE word: the 163")
print("  anchors are spread across the whole middle region.  A prefix-local theorem")
print("  would need anchors lying entirely before the earlier crossing j.")
print("  surviving site: j = %d  (j/K = %.4f)"%(j,j/K))
nb_pref=(j-17)//8
sixes_pref=a-11-5*nb_pref
print("  six-one bytes before j: %.4e  (of %.4e in the whole word)"
      %(float(sixes_pref), float(D-11-5*((K-17)//8))))
m_loc=(sixes_pref-PAST-FW-1+labels-1)//labels
print("  forced multiplicity before j:  m_local = %d   (vs 163 for the full word)"%m_loc)
elo=F(j)*F(int(Decimal(2).ln()*Decimal(10)**220),10**220)-F(a)*F(int(Decimal(3).ln()*Decimal(10)**220),10**220)
cap=F(7,3)+F(".94")*((t0-8)//8)
Qloc=cap/elo
print("  failure ceiling at j:  M_j/E <= %.6e"%float(Qloc))
Mloc=F(a,3)
Vloc=F("0.162")*Qloc+F(256,81)*Mloc
need=Vloc/2**L
print("  state width before j:  V_local = %.6e, so V_local/2^104 = %.6e"
      %(float(Vloc),float(need)))
print("  required: sigma_%d > %.6e"%(m_loc,float(need)))
print("  known:    sigma_163 = %d, and sigma_m is NONDECREASING in m, so"%SIGMA163)
print("            sigma_%d <= sigma_163 = %.6e"%(m_loc,float(SIGMA163)))
print("  => prefix-localization needs sigma_%d to reach %.2f%% of sigma_163."
      %(m_loc,100.0*float(need)/SIGMA163))
print("     Under equidistribution sigma_m ~ (m-1) 3^t/|R|, so the expected ratio is")
print("     (%d-1)/(163-1) = %.3f -- a factor %.2f short."
      %(m_loc,(m_loc-1)/162.0,float(need)/SIGMA163/((m_loc-1)/162.0)))
print("  The exact sigma_%d is computed separately (span75)."%m_loc)
print()
print("  WHAT THE SUFFIX SUPPLIES.  The 163-anchor floor uses six-one starts spread over")
print("  the ENTIRE middle region; before j only %.1f%% of them are available"
      %(100.0*float(sixes_pref)/float(D-11-5*((K-17)//8))))
print("  (j/K = %.3f), and the multiplicity drops from 163 to %d.  The suffix is"%(j/K,m_loc))
print("  therefore essential in a quantified way: it supplies the extra 163 - %d = %d"%(m_loc,163-m_loc))
print("  forced occurrences without which the span inequality has no margin.")
print()
print("  MINIMUM CONTINUATION NEEDED.  Solve for the prefix length f*K at which the")
print("  localized multiplicity alone would suffice, under the equidistribution model:")
for fr in (0.5,0.6,0.7,0.8,0.9,1.0):
    jj=int(fr*K); nbp=(jj-17)//8; sx=a_of(jj-2)-11-5*nbp
    mm=(sx-PAST-FW-1+labels-1)//labels
    est=SIGMA163*(mm-1)/162.0
    print("    prefix = %.0f%% of K : m_local = %4d, sigma estimate %.4e vs need %.4e -> %s"
          %(100*fr,mm,est,float(need),"OK" if est>float(need) else "short"))
print("  (the 'need' column is fixed by the surviving site; a longer prefix raises the")
print("   multiplicity but the ceiling M_j/E is what it is.)")
