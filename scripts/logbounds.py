"""Certified rational enclosures of ln 2 and ln 3, by the note's own recipe.

For x in {2,3} put z = (x-1)/(x+1).  Then ln x = 2 artanh(z) and
    L_x = 2 * sum_{h=0}^{95} z^{2h+1}/(2h+1),      U_x = L_x + 2 z^193 / (193 (1 - z^2)),
with L_x < ln x < U_x.  Everything is an exact Fraction; no float decides anything.
"""
from fractions import Fraction as F

def log_bounds(x, terms=96):
    z = F(x - 1, x + 1)
    L = F(0)
    for h in range(terms):
        L += z ** (2 * h + 1) / (2 * h + 1)
    L *= 2
    tail = 2 * z ** (2 * terms + 1) / ((2 * terms + 1) * (1 - z * z))
    return L, L + tail

L2, U2 = log_bounds(2)
L3, U3 = log_bounds(3)
