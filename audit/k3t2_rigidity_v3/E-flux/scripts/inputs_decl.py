"""Declared inputs of Track E v3 (rule 7).  Each has {name, value, tier, why}.  Scripts that use a value read it
from here (e.g. TADPOLE_TOTAL) rather than re-typing it.  Written to inputs.json by build_results.py.
Quotation line numbers are looked up by grep at build time (never typed)."""
from common import *
import json

TADPOLE_TOTAL = 24           # TT eq (2.3) right-hand side  (tier L)
HALF = 2                     # TT eq (2.3): (1/N) N_flux with N = 2   (tier L)
NFLUX_FACTOR_41 = 2          # TT (4.13): N_flux = 2 alpha_x^2 for the (4.8) family (tier L; equals the value derived from (2.8) in 02 -- checked, not assumed)

def inputs():
    q = lambda s: quote(s)
    return [
        {"name": "TT_H33_gram_A3", "value": "parsed 6x6 matrix (three hyperbolic planes)", "tier": "L", "why": "TT (A.3), parsed from pinned text by common.parse_H33", "source": q("where the matrix H3,3 is defined as")},
        {"name": "TT_E8_cartan_A4", "value": "parsed 8x8 Cartan matrix", "tier": "L", "why": "TT (A.4), parsed by common.parse_E8", "source": q("and E8 is the Catran matrix")},
        {"name": "TT_A6_gram_A6_A7", "value": "diag(2,2,2,-2,-2,-2)", "tier": "L", "why": "TT (A.6)-(A.7): (e_i,e_i)=+2 (i<=3), -2 (i>=4), others 0; typed here as sympy.diag", "source": q("(e4 , e4 ) = (e5 , e5 ) = (e6 , e6 ) = −2")},
        {"name": "embedding_A6_into_U3", "value": "e_i = f_i + g_i, e_{i+3} = f_i - g_i (i=1,2,3)", "tier": "declared structural choice", "why": "one isometric embedding of the (A.6) Gram into U^3; index |det| = 8 is forced for any isometric embedding (01)"},
        {"name": "flux_quantisation_even_coefficients", "value": "flux vector = 2*lambda, lambda in Gamma", "tier": "L", "why": "TT (2.5) and App A remark: even coefficients", "source": q("quantisation condition that αix is even integer")},
        {"name": "tadpole_total_24_eq2.3", "value": TADPOLE_TOTAL, "tier": "L", "why": "TT (2.3) right-hand side 24 (typed; NOT derived here; the tier-L bookkeeping 4*2+16*1 is in literature_checks)", "source": q("Nf lux + ND3 = 24")},
        {"name": "tadpole_half_factor_eq2.3", "value": "1/2", "tier": "L", "why": "TT (2.3): (1/2) N_flux + N_D3 = 24; decided against line-941 prose by 02", "source": q("Nf lux + ND3 = 24")},
        {"name": "N_flux_formula_2.8", "value": "N_flux = -beta_x.alpha_y + beta_y.alpha_x", "tier": "L", "why": "TT (2.8)", "source": q("(−βx · αy + βy · αx )")},
        {"name": "sec41_conditions_4.8_4.10", "value": "alpha_y=-beta_x, beta_y=alpha_x, alpha_x^2=beta_x^2, alpha_x.beta_x=0, V_flux type (2+,0-)", "tier": "L", "why": "TT section 4.1 eqs (4.8)-(4.10)", "source": q("αy = −βx and βy = αx")},
        {"name": "sec42_conditions_4.24_4.25", "value": "(4.24) alpha_x.alpha_y = beta_x.beta_y = alpha_x.beta_y+alpha_y.beta_x = 0; (4.25) beta_x^2 = 2 alpha_x^2 = 2 alpha_x.beta_x, same for y", "tier": "L", "why": "TT section 4.2; then phi=(1+-i)/2, tau = i sqrt(alpha_yy/alpha_xx)", "source": q("(αx · αy ) = (βx · βy ) = (αx · βy + αy · βx ) = 0")},
        {"name": "branch2_conditions_5.1_5.3", "value": "G_zbar=0 i.e. (alpha_x-phi beta_x) tau = alpha_y - phi beta_y; dim V_flux <= 2 spanned by timelike vectors", "tier": "L", "why": "TT section 5.1", "source": q("Gz̄ = nx τ − ny = 0")},
        {"name": "orbifold_criterion_TT_3.3", "value": "singular iff V_flux contains a lattice vector orthogonal to Omega", "tier": "L", "why": "TT section 3.3 'Orbifold Singularities'", "source": q("iff Vf lux contains a Lattice Vector v of Γ3,19")},
        {"name": "group_G_definition", "value": "(Z2 x Z2) wr S3 on U+U+U, order 384; G+ = det +1 on positive 3-block, order 192", "tier": "declared structural choice", "why": "explicit finite subgroup of O(Gamma); verified isometry in 01"},
        {"name": "window_bounds", "value": "entry bound N on lambda coordinates (recorded in every count)", "tier": "declared truncation", "why": "counts are counts inside a finite window; the tadpole does not bound entries in indefinite signature (control in 03)"},
        {"name": "TT_printed_examples_expected", "value": "(4.16) N_D3=16; (4.32) N_flux=32, 8 D3; (5.3) N_flux/2=4, 20 D3", "tier": "L (expected only)", "why": "used only as 'expected' fields after computing; line numbers by grep"},
        {"name": "T4Z2_fixed_points_from_track_D", "value": "read from Track D exports (poll_track_d.py)", "tier": "cross-track import", "why": "number of T^4/Z2 singular points; not typed"},
    ]

if __name__ == "__main__":
    print(json.dumps(inputs(), indent=1)[:2000])
