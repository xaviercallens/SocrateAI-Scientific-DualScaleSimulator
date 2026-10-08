"""Cross-check this repo's transcription of U + <2N> (leanflow.core.gamma0n_plus.gram_U_plus_2N) against
Stream 2's certificate K3xT2_READING_S.json (K3-DarkMatter, added at commit 5fe0f98, ruling D29', 2026-10-08).

What the certificate says (Stream 2, exact computation, its Tier B): for elliptic curves E, E' without CM joined by a cyclic
n-isogeny, the transcendental lattice of E x E' is isometric to each family's certified T, for n = 7 (cooper_s7) and
n = 10 (cooper_s10). This test checks only that our Gram matrix for U + <2n> is integrally isometric to the Gram matrix the
certificate stores, by an explicit unimodular change of basis found here, and that a wrong n is NOT isometric.

WHAT THIS DOES NOT DO (decision D2 of audit/DUALSCALE_TO_STREAMS_DECISIONS_2026_09_28.md stands): it does not choose a
level, it takes no parameter from E or E', and nothing here is a simulation input. Reading S is Stream 2's convention for
selection criteria only; the physical reading of K3 x T^2 stays Tier C.

Skips (does not fail) when the Stream 2 certificate is not present on this machine; set K3XT2_READING_S_JSON to its path.
"""
import itertools
import json
import os

import numpy as np
import pytest

from leanflow.core.gamma0n_plus import gram_G0N, gram_U_plus_2N

CERT = os.environ.get(
    "K3XT2_READING_S_JSON",
    os.path.expanduser("~/SocrateAI-Scientific-Agora-K3-DarkMatter/data/certificates/K3xT2_READING_S.json"),
)


def _cert():
    if not os.path.exists(CERT):
        pytest.skip(f"Stream 2 certificate not present at {CERT}")
    return json.load(open(CERT))


def _isometry(G_from, G_to, box=2):
    """unimodular C with C G_from C^T == G_to, by brute force over small entries; None if none found."""
    A, B = np.array(G_from, dtype=np.int64), np.array(G_to, dtype=np.int64)
    if round(np.linalg.det(A)) != round(np.linalg.det(B)):
        return None
    rows = [np.array(r) for r in itertools.product(range(-box, box + 1), repeat=3) if any(r)]
    by_norm = {}
    for r in rows:
        by_norm.setdefault(int(r @ A @ r), []).append(r)
    for r0 in by_norm.get(int(B[0, 0]), []):
        for r1 in by_norm.get(int(B[1, 1]), []):
            if int(r0 @ A @ r1) != B[0, 1]:
                continue
            for r2 in by_norm.get(int(B[2, 2]), []):
                C = np.array([r0, r1, r2])
                if abs(round(np.linalg.det(C))) == 1 and np.array_equal(C @ A @ C.T, B):
                    return C
    return None


def test_P1_certified_T_is_isometric_to_our_U_plus_2N():
    d = _cert()["result"]["generic"]
    for fam, v in d.items():
        n = v["n"]
        assert v["verdict"] == "ISOMETRIC_TO_CERTIFIED_T", (fam, v["verdict"])
        C = _isometry(v["T_gram"], gram_U_plus_2N(n))
        assert C is not None, f"{fam}: no isometry between the certificate's T and gram_U_plus_2N({n})"


def test_N1_wrong_level_is_not_isometric():
    d = _cert()["result"]["generic"]
    for fam, v in d.items():
        n = v["n"]
        for wrong in (n - 1, n + 1, 12):
            assert _isometry(v["T_gram"], gram_U_plus_2N(wrong)) is None, (fam, wrong)


def test_N2_the_other_rank3_lattice_is_not_isometric():
    """Stream 1's no_isometry_G0N_TN: <1> + U(2N) is a different lattice; the cross-check must not confuse them."""
    d = _cert()["result"]["generic"]
    for fam, v in d.items():
        assert _isometry(v["T_gram"], gram_G0N(v["n"])) is None, fam
