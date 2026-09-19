"""Fetch the external linear-chromosome control for test 2 (run with the separate venv that has hic-straw):
  /mnt/disks/disk-socrateai-local-1/venv-tdaval/bin/python fetch_gm12878_chr1q.py
Rao et al. 2014 (GSE63525) GM12878 in-situ combined .hic, read remotely; chr1 q-arm 145,000,000-249,250,621 (hg19),
250 kb bins, KR-normalised observed contacts. Writes a dense .npy to the data directory (not git)."""
import json, os, hashlib
import numpy as np
import hicstraw
URL = "https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/combined.hic"
START, END, RES = 145_000_000, 249_250_621, 250_000
OUT = "/mnt/disks/disk-socrateai-local-1/dualscale-data-r3/tda_validation/genetics/gm12878"
os.makedirs(OUT, exist_ok=True)
region = f"1:{START}:{END}"
recs = hicstraw.straw("observed", "KR", URL, region, region, "BP", RES)
b0 = START // RES
nb = (END - 1) // RES - b0 + 1
M = np.full((nb, nb), np.nan)
for r in recs:
    i, j = r.binX // RES - b0, r.binY // RES - b0
    if 0 <= i < nb and 0 <= j < nb:
        M[i, j] = M[j, i] = r.counts
p = os.path.join(OUT, "gm12878_insitu_combined_chr1_145000000_249250621_250kb_KR_observed.npy")
np.save(p, M)
meta = {"url": URL, "hicstraw_call": f"hicstraw.straw('observed','KR',URL,'{region}','{region}','BP',{RES})",
        "n_records": len(recs), "shape": list(M.shape), "first_bin_start_bp": b0 * RES,
        "n_nan": int(np.isnan(M).sum()), "file": p,
        "sha256": hashlib.sha256(open(p, "rb").read()).hexdigest(),
        "genome": hicstraw.HiCFile(URL).getGenomeID()}
json.dump(meta, open(os.path.join(OUT, "fetch_meta.json"), "w"), indent=2)
print(meta)
