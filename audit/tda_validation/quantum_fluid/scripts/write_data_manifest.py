"""Write audit/tda_validation/quantum_fluid/data_manifest.json: sha256, size and
provenance (source URL or generating command) of every file stored under
DATA_ROOT (outside git), plus the local file assessed but not used.
Command: prlimit --as=8589934592 -- .venv-tda/bin/python write_data_manifest.py
"""
import hashlib
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qf_common as qc  # noqa: E402

PY_TDA = "/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python"
PY_MC = "/mnt/disks/disk-socrateai-local-1/venv-tdaval/bin/python"
S = "audit/tda_validation/quantum_fluid/scripts"
ZEN = "https://zenodo.org/api/records/14780459/files/"


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def provenance(rel):
    parts = rel.split("/")
    if rel.startswith("raw/zenodo_14780459/STM_data.zip"):
        return {"source_url": ZEN + "STM%20data.zip/content", "zenodo_md5": "5b97a8b60660f1f7b521f0f7a526502f",
                "license": "CC-BY-4.0", "doi": "10.5281/zenodo.14780459"}
    if rel.startswith("raw/zenodo_14780459/raw_image_addition.png"):
        return {"source_url": ZEN + "raw_image_addition.png/content", "zenodo_md5": "5033419cd4911d1024eb264891eed22d"}
    if rel.startswith("raw/zenodo_14780459/extracted/"):
        return {"derived_from": "raw/zenodo_14780459/STM_data.zip", "command": "unzip -q -o STM_data.zip -d extracted"}
    if parts[0] == "xy":
        L = parts[1][1:]
        T = parts[2][1:-4]
        return {"command": f"prlimit --as=8589934592 -- {PY_MC} {S}/xy_mc.py --L {L} --temps <grid incl. {T}> --ntherm 2000 --nmeas 400 --gap 2 --save-configs 1",
                "seed": f"20260919 + 1000*{L} + round(1000*{T})"}
    if parts[0] == "xy_tda" and parts[1].endswith("_shufonly"):
        L = parts[1][1:].split("_")[0]
        return {"command": f"prlimit --as=8589934592 -- {PY_TDA} {S}/xy_tda_features.py --L {L} --stride 4 --shuffle-only 50"}
    if parts[0] == "xy_tda":
        L = parts[1][1:]
        extra = "--stride 2 --shuffle-n 50 --alpha-n 100" if L in ("32", "64") else "--stride 4 --shuffle-n 0 --alpha-n 0"
        return {"command": f"prlimit --as=8589934592 -- {PY_TDA} {S}/xy_tda_features.py --L {L} {extra}"}
    if parts[0] == "gpe" and parts[1].startswith("psi_") and parts[1].endswith("_ext.npz"):
        om = parts[1][len("psi_Omega"):-len("_ext.npz")]
        return {"command": f"prlimit --as=8589934592 -- {PY_TDA} {S}/gpe_rotating.py --omega {om} --resume {qc.DATA_ROOT}/gpe/psi_Omega{om}.npz --max1 0 --max2 40000 --tag _ext"}
    if parts[0] == "gpe" and parts[1].startswith("psi_"):
        om = parts[1][len("psi_Omega"):-4]
        return {"command": f"prlimit --as=8589934592 -- {PY_TDA} {S}/gpe_rotating.py --omega {om}", "seed": 7}
    if parts[0] == "gpe":
        return {"command": f"prlimit --as=8589934592 -- {PY_TDA} {S}/gpe_tda.py"}
    if parts[0] == "stm":
        return {"command": f"prlimit --as=8589934592 -- {PY_TDA} {S}/stm_vortex_tda.py"}
    return {"note": "cache / other"}


def main():
    files = []
    for root, dirs, fs in os.walk(qc.DATA_ROOT):
        if "numba_cache" in root:
            continue
        for f in sorted(fs):
            p = os.path.join(root, f)
            rel = os.path.relpath(p, qc.DATA_ROOT)
            files.append({"path": p, "size_bytes": os.path.getsize(p), "sha256": sha(p), **provenance(rel)})
    files.sort(key=lambda d: d["path"])
    loc = "/mnt/disks/disk-socrateai-local-1/dual_scale_datasets/domain06_supercon_abrikosov_vortices.npz"
    man = {"data_root": qc.DATA_ROOT, "n_files": len(files), "files": files,
           "assessed_not_used": {"path": loc, "sha256": sha(loc), "size_bytes": os.path.getsize(loc),
                                 "note": "content assessed in results/domain06_assessment.json; not used as validation data"},
           "generated_by": f"prlimit --as=8589934592 -- {PY_TDA} {S}/write_data_manifest.py"}
    out = os.path.abspath(os.path.join(qc.HERE, "..", "data_manifest.json"))
    with open(out, "w") as fh:
        json.dump(man, fh, indent=1)
    print(out, len(files))


if __name__ == "__main__":
    main()
