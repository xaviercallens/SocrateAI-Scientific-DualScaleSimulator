# Exact fetch commands (data/real2/)

Run 2026-09-19, worktree `/mnt/disks/disk-socrateai-local-1/dualscale-wt-reverse`, branch `loop/reverse-zero`.
Every download was checked with `file`/`head -c`/`wc -c` before hashing (check-then-hash).

## 1. DESI DR2 BAO

```
curl -s --max-time 30 "https://api.github.com/repos/CobayaSampler/bao_data/contents/desi_bao_dr2" -o dr2_listing.json
curl -s --max-time 30 "https://api.github.com/repos/CobayaSampler/bao_data/contents/" -o root_listing.json
# then, for each entry's download_url in dr2_listing.json (16 files: ALL_GCcomb + 7 per-tracer mean+cov pairs):
curl -s --max-time 30 "<download_url>" -o "desi_dr2/<name>"
```

## 2. Pantheon+ full STAT+SYS covariance

```
curl -sL --max-time 180 \
  "https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/Pantheon%2BSH0ES_STAT%2BSYS.cov" \
  -o "Pantheon+SH0ES_STAT+SYS.cov"
```

## 3. NANOGrav 15yr HD angular correlation — ABSENT (see MANIFEST.json for full URL trail)

## 4. SDSS DR17 cosmic web

```
# row-count check
curl -s --max-time 60 "https://skyserver.sdss.org/dr17/SkyServerWS/SearchTools/SqlSearch?cmd=SELECT%20COUNT%28%2A%29%20as%20n%20FROM%20SpecObj%20WHERE%20class%3D%27GALAXY%27%20AND%20zWarning%3D0%20AND%20z%20BETWEEN%200.02%20AND%200.12%20AND%20ra%20BETWEEN%20140%20AND%20220%20AND%20dec%20BETWEEN%200%20AND%2050&format=csv" -o count.csv
# -> 193536

# full pull (no paging needed: row count matched exactly, no cap hit)
curl -s --max-time 180 "https://skyserver.sdss.org/dr17/SkyServerWS/SearchTools/SqlSearch?cmd=SELECT%20ra%2Cdec%2Cz%2CzErr%20FROM%20SpecObj%20WHERE%20class%3D%27GALAXY%27%20AND%20zWarning%3D0%20AND%20z%20BETWEEN%200.02%20AND%200.12%20AND%20ra%20BETWEEN%20140%20AND%20220%20AND%20dec%20BETWEEN%200%20AND%2050&format=csv" -o sdss_dr17_galaxies_ra140_220_dec0_50.csv

# random comparison catalogue (seed=20260919), see cosmic_web/make_random_catalogue.py:
/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python make_random_catalogue.py
```

## 5. WMAP 9yr ILC + KQ85 mask

```
curl -s --max-time 180 "https://lambda.gsfc.nasa.gov/data/map/dr5/dfp/ilc/wmap_ilc_9yr_v5.fits" -o wmap_ilc_9yr_v5.fits
curl -s --max-time 120 "https://lambda.gsfc.nasa.gov/data/map/dr5/ancillary/masks/wmap_temperature_kq85_analysis_mask_r9_9yr_v5.fits" -o wmap_temperature_kq85_analysis_mask_r9_9yr_v5.fits
```

Both are HEALPix FITS, 25174080 bytes each; both confirmed with `file` as "FITS image data" before hashing.
404 paths tried and rejected for the mask before the one that worked:
`.../data/map/dr5/ancillary/masks/wmap_temperature_analysis_mask_r9_9yr_v5.fits` (404),
`.../data/map/dr5/masks/wmap_temperature_analysis_mask_r9_9yr_v5.fits` (404).

## 6. Fifth-force arXiv:2002.11761 source

```
curl -sL --max-time 60 "https://arxiv.org/e-print/2002.11761" -o 2002.11761.tar.gz
tar xzf 2002.11761.tar.gz -C extracted
find extracted -type f
grep -n "tabular\|pgfplotstable\|begin{table}" extracted/FB_ISL_pdf.tex
```
Result: no data tables in the .tex, only 9 figure PDFs -> ABSENT for machine-readable exclusion table (see MANIFEST.json). Tarball itself kept on disk as provenance.

Note: the first `curl -s` (no `-L`) attempt returned an HTML redirect stub (217 bytes, `Redirecting...` to `/src/2002.11761`), caught by the check-then-hash gate before any hash was taken; the retry with `-L` produced the real gzip tarball.
