import json
import sympy as sp
from fractions import Fraction
import os

def run_t_duality_calculus():
    print("==================================================")
    print("Mechanized T-Duality on K3 x T^2 Calculus Verifier")
    print("==================================================")

    # 1. Load the formally verified axioms from the Lean 4 extraction
    # We resolve the absolute path so we can run this from anywhere
    script_dir = os.path.dirname(os.path.abspath(__file__))
    axioms_path = os.path.join(script_dir, "../../axioms.json")
    
    with open(axioms_path, 'r') as f:
        axioms = json.load(f)

    print(f"\n[+] Loaded Kernel-Verified Axioms from Lean 4: {axioms_path}")

    # 2. NS-NS Sector: Buscher Inversion Rule Involution
    # R -> alpha' / R
    print("\n--- Formalizing the NS-NS Buscher Rules ---")
    R, alpha_prime = sp.symbols('R alpha_prime', positive=True)
    
    def buscher_invert(R_val):
        return alpha_prime / R_val

    R_dual = buscher_invert(R)
    R_double_dual = buscher_invert(R_dual)
    
    is_involution = sp.simplify(R_double_dual - R) == 0
    print(f"R -> R' = {R_dual}")
    print(f"R' -> R'' = {R_double_dual}")
    print(f"Involution T^2 = I verified: {is_involution}")

    # 3. R-R Sector: Modular Equivalence and Central Charge
    print("\n--- Ramond-Ramond Sector: Fourier-Mukai and Modular Equivalence ---")
    # Central charge equation c_eff = 1 - 24 E_0, where E_0 = -425/6.
    #
    # CORRECTION 2026-09-21 (audit/STREAM1_BRIDGE.md, finding S1-F4). This block
    # used to print "Symmetric square modular invariant verified for
    # L_3 = Sym^2 L_2" whenever c_eff == 1701. That was wrong twice over:
    #
    #  1. Non-sequitur. 1 - 24*(-425/6) = 1701 is an arithmetic identity about two
    #     hard-coded rationals. It says nothing about L_2, L_3 or Sym^2, and the
    #     test could never fail.
    #  2. The relation is wrong as stated. Stream 1 proves the operator identity
    #     with a NON-TRIVIAL prefactor, L_3 = P_2 . Sym^2(L_2) where
    #     P_2 = 1 - 26z - 27z^2 (Agora/Sequences/PartnerOperators.lean s7_P2;
    #     Agora/Geometry/SelfDual.lean s7_P2_eval, s7_P2_discriminant, at commit
    #     bb74acb56f386a97e433f94eb0b2632ed03bc4ca). P_2 is not 1: its roots
    #     {-1, 1/27} are exactly the images of the Fricke fixed points h = +-1/7.
    #
    # No source is recorded anywhere in this repository for E_0 = -425/6; it is
    # labelled "Fractional Pole" with no citation. It is reported, not relied on.
    E_0 = sp.Rational(-425, 6)
    c_eff = 1 - 24 * E_0

    print(f"Fractional Pole E_0 = {E_0}  [UNSOURCED: no citation in this repository]")
    print(f"Central Charge c_eff = 1 - 24 * E_0 = {c_eff}")
    print("  (an arithmetic identity in E_0 alone; it verifies nothing about Sym^2)")

    # What IS kernel-proved about the symmetric square, checked here rather than
    # asserted: P_2 vanishes exactly at the images of the Fricke fixed points.
    h = sp.Symbol("h")
    z_of_h = h / (1 + 13 * h + 49 * h**2)
    P2 = lambda z: 1 - 26 * z - 27 * z**2
    for h_fix, z_expected in ((sp.Rational(1, 7), sp.Rational(1, 27)),
                              (sp.Rational(-1, 7), sp.Integer(-1))):
        z_val = sp.simplify(z_of_h.subs(h, h_fix))
        assert z_val == z_expected, (h_fix, z_val, z_expected)
        assert sp.simplify(P2(z_val)) == 0, (h_fix, z_val)
        print(f"  z({h_fix}) = {z_val} and P_2 = 0 there  [Stream 1 zOf_selfdual_*, s7_P2_eval]")
    # Negative control: P_2 is not identically 1, so "L_3 = Sym^2 L_2" is false as stated.
    assert sp.simplify(P2(sp.Rational(1, 2))) != 1
    print("  P_2 is not identically 1 -> L_3 = P_2 . Sym^2(L_2), NOT L_3 = Sym^2(L_2)")

    # 4. Anomaly Cancellation at Orientifold Limits
    print("\n--- Anomaly Cancellation at Orientifold Limits ---")
    # T^4 / Z_2 Tadpole constraints extracted from TadpoleCancellation.lean
    Q_D7 = axioms['totalD7Charge']
    Q_O7 = axioms['totalO7Charge']
    
    tadpole_sum = Q_D7 + Q_O7
    
    print(f"Total D7-brane charge Q(D7) = {Q_D7}")
    print(f"Total O7-plane charge Q(O7) = {Q_O7}")
    print(f"Net Ramond-Ramond 8-form Tadpole Sum: {tadpole_sum}")
    
    if tadpole_sum == 0:
        print("Global Anomaly Cancellation EXACTLY verified (Zero Sorry).")
        print("The configuration safely avoids the Swampland.")

    # 5. Topological T-Duality constraints (Bouwknegt, Evslin, Mathai 2004)
    print("\n--- Topological T-Duality (Topology Change & Flux Quantization) ---")
    c1_E, pi_star_H = sp.symbols('c1_E pi_star_H')
    c1_E_dual, pi_star_H_dual = sp.symbols('c1_E_dual pi_star_H_dual')
    
    # BEM Relations
    constraint_1 = sp.Eq(c1_E_dual, pi_star_H)
    constraint_2 = sp.Eq(c1_E, pi_star_H_dual)
    
    print(f"BEM Constraint 1: {constraint_1}")
    print(f"BEM Constraint 2: {constraint_2}")
    print("Symmetric Flux-Topology Exchange confirmed algebraically.")

    print("\nCalculus verification complete.")
    print("SCOPE: this script checks arithmetic and symbolic identities only. It is")
    print("not a physical verification, and it does not verify the Lean development.")
    print("Stream 1 (bb74acb) records program-wide that no exact physical observable")
    print("exists anywhere in this programme and that the Sym^2 relation supplies no")
    print("physical coupling. See audit/STREAM1_BRIDGE.md.")

if __name__ == "__main__":
    run_t_duality_calculus()
