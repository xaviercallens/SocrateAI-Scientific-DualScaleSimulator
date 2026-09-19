"""X4 step 0: fetch the CMB distance-prior source (Chen, Huang & Wang, arXiv:1808.05724).
Downloads the arXiv abstract page and the TeX source (e-print), records sha256, extracts only the
.tex file (the PDF figures and tarball are not kept in git), and writes data/fetch_log.json.
Up to 3 attempts per URL. Run:
  cd audit/reverse_zero_r2/X4-bao-cmb-sn && <venv python> x4_fetch_prior.py
"""
import hashlib, json, pathlib, re, tarfile, time, urllib.request, io
HERE = pathlib.Path(__file__).resolve().parent
DATA = HERE / "data"; DATA.mkdir(exist_ok=True)
URLS = {"abs": "https://arxiv.org/abs/1808.05724",
        "pdf": "https://arxiv.org/pdf/1808.05724",
        "eprint": "https://arxiv.org/e-print/1808.05724"}
log = {"attempts": {}, "sha256": {}, "bytes": {}}
blobs = {}
for k, u in URLS.items():
    for attempt in range(1, 4):
        try:
            req = urllib.request.Request(u, headers={"User-Agent": "Mozilla/5.0 (X4 audit fetch)"})
            b = urllib.request.urlopen(req, timeout=90).read()
            blobs[k] = b; log["attempts"][k] = attempt
            log["sha256"][k] = hashlib.sha256(b).hexdigest(); log["bytes"][k] = len(b)
            break
        except Exception as e:
            log["attempts"][k] = f"failed attempt {attempt}: {e!r}"; time.sleep(3)
    else:
        log["sha256"][k] = "ABSENT"
html = blobs["abs"].decode("utf8", "replace")
meta = {m: re.findall(r'name="citation_%s" content="([^"]*)"' % m, html) for m in ["title", "author", "doi", "date", "arxiv_id"]}
log["arxiv_metadata"] = meta
tf = tarfile.open(fileobj=io.BytesIO(blobs["eprint"]), mode="r:gz")
tex = [m for m in tf.getmembers() if m.name.endswith(".tex")]
assert len(tex) == 1, [m.name for m in tex]
(DATA / "distance-priors-2018.tex").write_bytes(tf.extractfile(tex[0]).read())
log["tex_member"] = tex[0].name
log["tex_sha256"] = hashlib.sha256((DATA / "distance-priors-2018.tex").read_bytes()).hexdigest()
(DATA / "fetch_log.json").write_text(json.dumps(log, indent=1))
print(json.dumps(log, indent=1))
