"""Build a random (unclustered) comparison catalogue for TDA on the SDSS DR17
galaxy sample fetched from SkyServer (sdss_dr17_galaxies_ra140_220_dec0_50.csv).

Method (documented, per GROUND RULES -- construction lives here AND in
MANIFEST.json, not only in code):
  1. Read the real catalogue's RA, Dec, z, zErr columns.
  2. Keep RA, Dec exactly as observed (this preserves the survey's real
     angular footprint / mask -- edges, holes, fiber-collision structure --
     without needing a separate mask file).
  3. Independently shuffle the z column (Fisher-Yates via numpy permutation)
     so each galaxy's redshift is reassigned to a different sky position.
     zErr is shuffled with the same permutation so it stays attached to its z.
  4. This destroys real 3-D clustering (the cosmic web) while leaving the
     angular selection function and the redshift distribution n(z) intact,
     which is the standard "randoms from shuffled-z" construction used when
     a survey random catalogue is not otherwise available.

Seed: 20260919 (fixed, stated here, in the output JSON note, and in
data/real2/MANIFEST.json).
"""
import csv
import hashlib
import json
import sys

import numpy as np

SEED = 20260919
IN_PATH = "sdss_dr17_galaxies_ra140_220_dec0_50.csv"
OUT_PATH = "sdss_dr17_random_shuffled_z_ra140_220_dec0_50.csv"


def main():
    rows = []
    with open(IN_PATH, newline="") as f:
        r = csv.reader(f)
        header_marker = next(r)  # '#Table1'
        header = next(r)  # ra,dec,z,zErr
        assert header == ["ra", "dec", "z", "zErr"], header
        for row in r:
            rows.append(row)

    ra = np.array([float(x[0]) for x in rows])
    dec = np.array([float(x[1]) for x in rows])
    z = np.array([float(x[2]) for x in rows])
    zerr = np.array([float(x[3]) for x in rows])

    rng = np.random.default_rng(SEED)
    perm = rng.permutation(len(z))
    z_shuf = z[perm]
    zerr_shuf = zerr[perm]

    with open(OUT_PATH, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["#Table1"])
        w.writerow(["ra", "dec", "z", "zErr"])
        for i in range(len(ra)):
            w.writerow([ra[i], dec[i], z_shuf[i], zerr_shuf[i]])

    sha = hashlib.sha256(open(OUT_PATH, "rb").read()).hexdigest()
    summary = {
        "input": IN_PATH,
        "output": OUT_PATH,
        "n_rows": int(len(ra)),
        "seed": SEED,
        "method": "RA,Dec kept as observed; z (and paired zErr) independently "
        "permuted via numpy.random.default_rng(seed).permutation over all rows",
        "sha256_output": sha,
    }
    print(json.dumps(summary, indent=2))
    with open("random_catalogue_summary.json", "w") as f:
        json.dump(summary, f, indent=2)


if __name__ == "__main__":
    sys.exit(main())
