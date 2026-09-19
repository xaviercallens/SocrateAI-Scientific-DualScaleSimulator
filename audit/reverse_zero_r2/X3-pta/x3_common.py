"""Shared constants and helpers for X3 (NANOGrav 15-yr c4_pta_product bound).

Data-blind pulsar selection: uses only tim-file MJDs (no residuals, no cross-correlations).
"""
import pathlib
import re

HERE = pathlib.Path(__file__).resolve().parent
REPO = HERE.parents[2]  # .../audit/reverse_zero_r2/X3-pta -> repo root
CACHE = HERE / "cache"

# Data live outside the worktree (registration path data/real2/pulsar_timing/... is empty).
DATA_ROOT = pathlib.Path(
    "/mnt/disks/disk-socrateai-local-1/NANOGrav15yr/NANOGrav15yr_PulsarTiming_v2.1.0/narrowband"
)
EPHEM = "DE440"
NFREQ = 14                 # registration: 14 frequencies
GAMMA_CURN = 13.0 / 3.0    # registration: gamma = 13/3
# Fixed CURN amplitude. NOT ATTEMPTED: fetching/citing the exact NANOGrav 15-yr gamma=13/3 CURN
# amplitude from the paper this round; -14.62 is used as a representative fixed value only, not
# attributed to any specific NANOGrav table. A sensitivity run at a different log10_A
# (x3_persr.py --log10a) is reported alongside to show c4_pta_product is stable to this choice.
LOG10A_CURN = -14.62
MIN_SPAN_YR = 3.0
DUP_SUFFIXES = ("ao", "gbt")   # per-telescope splits of the same pulsar


def find_files(name):
    par = sorted(DATA_ROOT.glob(f"par/{name}_PINT_*.nb.par"))
    tim = sorted(DATA_ROOT.glob(f"tim/{name}_PINT_*.nb.tim"))
    assert len(par) == 1 and len(tim) == 1, (name, par, tim)
    return par[0], tim[0]


def all_names():
    names = sorted(p.name.split("_PINT_")[0] for p in DATA_ROOT.glob("par/*_PINT_*.nb.par"))
    return names


def is_split_duplicate(name, names):
    for s in DUP_SUFFIXES:
        if name.endswith(s) and name[: -len(s)] in names:
            return True
    return False


def tim_mjds(tim_path):
    out = []
    with open(tim_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("C ", "#", "FORMAT", "MODE", "TIME", "EFAC", "EQUAD")):
                continue
            t = line.split()
            if len(t) > 3:
                out.append(float(t[2]))
    return out


def chain_median_noise(name, burn=0.25):
    import numpy as np
    pars = (DATA_ROOT / "noise" / f"{name}.nb.pars.txt").read_text().split()
    ch = np.loadtxt(DATA_ROOT / "noise" / f"{name}.nb.chain_1.txt")
    ncol = len(pars)
    ch = ch[int(burn * len(ch)):, :ncol]
    return {p: float(np.median(ch[:, i])) for i, p in enumerate(pars)}
