"""Chain K=O, b1=0 -> chi(O) -> chi_top (Noether) -> tau (Hirzebruch) -> (b2+,b2-) (Hodge index). All exact Fractions.
Every step uses a declared input from decl.py; nothing typed as a target."""
from fractions import Fraction
def chain():
    c1sq = 0                       # K = O => c1 = 0 (K_trivial)
    b1 = 0                         # K_trivial
    h0 = 1                         # hodge_h0_O
    h1 = Fraction(b1, 2)           # hodge_h1_O_from_b1
    h2 = h0                        # serre_duality with K = O
    chi_O = h0 - h1 + h2           # computed: 2
    NOETHER = 12                   # noether_formula (declared)
    c2 = NOETHER * chi_O - c1sq    # chi_top
    tau = Fraction(c1sq - 2 * c2, 3)   # hirzebruch_signature
    b3 = b1; b0 = b4 = 1
    b2 = c2 - b0 - b4 + b1 + b3    # chi_top = b0-b1+b2-b3+b4
    bp = (b2 + tau) / 2; bm = (b2 - tau) / 2   # hodge_index
    assert bp.denominator == 1 and bm.denominator == 1
    # second Hodge-diamond route to tau: p_g = h^2(O)=1, h11 = b2 - 2 p_g
    pg = h2; h11 = b2 - 2 * pg
    diamond = {(0,0):1,(0,1):0,(0,2):pg,(1,0):0,(1,1):h11,(1,2):0,(2,0):pg,(2,1):0,(2,2):1}
    tau2 = sum((-1) ** p * v for (p, q), v in diamond.items())
    return {"chi_O": chi_O, "chi_top": int(c2), "b2": int(b2), "tau": int(tau), "tau_diamond": int(tau2),
            "b2plus": int(bp), "b2minus": int(bm), "h11": int(h11), "p_g": int(pg)}
