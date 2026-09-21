#!/usr/bin/env python3
"""
publish_to_zenodo.py
Zenodo automated deposit and publication pipeline for LeanFlow & K3 x T2 String Theory manuscript.

Publishes:
  - T_duality_Alone.pdf (34-page compiled manuscript)
  - T_duality_Alone.tex (LaTeX source)
  - lean4_formal_proofs.tar.gz (14 sorry-free Lean 4 proof certificates)
  - rust_stiff_dualscale_simulator.tar.gz (Rust stiff solver crate)
  - tda_mapper_skeleton.json (187 nodes, 557 edges, beta_1=376)
  - kummer_langevin_summary.json (48x48 Kummer Langevin telemetry)
  - vacuum_decay_cdl_summary.json (Coleman-De Luccia bounce action telemetry)
  - tachyon_condensation_summary.json (Boundary tachyon condensation telemetry)
  - zenodo_deposit_bundle.zip (Unified bundle)
"""

import os
import sys
import json
import argparse
import urllib.request
import urllib.error
from pathlib import Path

DEFAULT_TOKEN_FILE = os.path.expanduser("~/.config/zenodo/token")
WORKSPACE_DIR = Path(__file__).resolve().parent.parent
BUNDLE_DIR = WORKSPACE_DIR / "zenodo_bundle"
ZENODO_BASE_URL = "https://zenodo.org/api"
SANDBOX_BASE_URL = "https://sandbox.zenodo.org/api"

def get_token():
    token = os.environ.get("ZENODO_TOKEN")
    if token:
        return token.strip()
    if os.path.exists(DEFAULT_TOKEN_FILE):
        with open(DEFAULT_TOKEN_FILE, "r") as f:
            return f.read().strip()
    return None

def load_metadata():
    meta_path = BUNDLE_DIR / ".zenodo.json"
    if not meta_path.exists():
        meta_path = WORKSPACE_DIR / ".zenodo.json"
    with open(meta_path, "r") as f:
        return json.load(f)

def make_request(url, method="GET", data=None, headers=None, token=None, proxy=None):
    req_headers = {
        "Accept": "application/json",
        "User-Agent": "ZenodoPublicationClient/1.0 (LeanFlow-Academic; SocrateAI)"
    }
    if token:
        req_headers["Authorization"] = f"Bearer {token}"
    if headers:
        req_headers.update(headers)

    body = None
    if data is not None:
        if isinstance(data, (dict, list)):
            body = json.dumps(data).encode("utf-8")
            req_headers["Content-Type"] = "application/json"
        elif isinstance(data, (bytes, bytearray)):
            body = data
        else:
            body = str(data).encode("utf-8")

    req = urllib.request.Request(url, data=body, headers=req_headers, method=method)
    
    opener_handlers = []
    if proxy:
        opener_handlers.append(urllib.request.ProxyHandler({"http": proxy, "https": proxy}))
    
    opener = urllib.request.build_opener(*opener_handlers)
    try:
        with opener.open(req, timeout=60) as resp:
            content_type = resp.headers.get("Content-Type", "")
            raw = resp.read()
            if "application/json" in content_type:
                return json.loads(raw.decode("utf-8")), resp.status
            return raw, resp.status
    except urllib.error.HTTPError as e:
        err_body = e.read().decode("utf-8", errors="replace")
        print(f"[-] HTTP Error {e.code} for {url}: {err_body[:300]}", file=sys.stderr)
        raise
    except urllib.error.URLError as e:
        print(f"[-] Network URL Error for {url}: {e.reason}", file=sys.stderr)
        raise

def create_deposit(base_url, token, proxy=None):
    url = f"{base_url}/deposit/depositions"
    print(f"[*] Initializing new deposition at {url}...")
    res, status = make_request(url, method="POST", data={}, token=token, proxy=proxy)
    dep_id = res["id"]
    bucket_url = res.get("links", {}).get("bucket")
    print(f"[+] Deposition created successfully. ID: {dep_id}")
    if bucket_url:
        print(f"[+] Direct file bucket URL: {bucket_url}")
    return dep_id, bucket_url, res

def upload_file(base_url, dep_id, bucket_url, filepath, token, proxy=None):
    filename = os.path.basename(filepath)
    size_bytes = os.path.getsize(filepath)
    size_mb = size_bytes / (1024 * 1024)
    print(f"[*] Uploading '{filename}' ({size_mb:.2f} MB)...")

    with open(filepath, "rb") as f:
        file_bytes = f.read()

    if bucket_url:
        # Preferred modern bucket upload
        upload_url = f"{bucket_url}/{filename}"
        res, status = make_request(
            upload_url,
            method="PUT",
            data=file_bytes,
            headers={"Content-Type": "application/octet-stream"},
            token=token,
            proxy=proxy
        )
        print(f"[+] Uploaded '{filename}' (Status {status})")
        return res
    else:
        # Fallback multipart / deposition files endpoint
        upload_url = f"{base_url}/deposit/depositions/{dep_id}/files"
        # Multipart form data upload
        boundary = "----WebKitFormBoundaryZenodoUploader7MA4YWxkTrZu0gW"
        body = (
            f"--{boundary}\r\n"
            f'Content-Disposition: form-data; name="file"; filename="{filename}"\r\n'
            f"Content-Type: application/octet-stream\r\n\r\n"
        ).encode("utf-8") + file_bytes + f"\r\n--{boundary}--\r\n".encode("utf-8")

        res, status = make_request(
            upload_url,
            method="POST",
            data=body,
            headers={"Content-Type": f"multipart/form-data; boundary={boundary}"},
            token=token,
            proxy=proxy
        )
        print(f"[+] Uploaded '{filename}' via fallback (Status {status})")
        return res

def set_metadata(base_url, dep_id, metadata, token, proxy=None):
    url = f"{base_url}/deposit/depositions/{dep_id}"
    print(f"[*] Updating deposition metadata for ID {dep_id}...")
    payload = {"metadata": metadata}
    res, status = make_request(url, method="PUT", data=payload, token=token, proxy=proxy)
    print(f"[+] Metadata registered successfully: '{metadata.get('title')}'")
    return res

def create_new_version(base_url, record_id, token, proxy=None):
    """Open a new-version draft of a published record so the concept DOI keeps one lineage."""
    url = f"{base_url}/deposit/depositions/{record_id}/actions/newversion"
    print(f"[*] Creating new-version draft of record {record_id}...")
    res, _ = make_request(url, method="POST", data={}, token=token, proxy=proxy)
    draft_url = res.get("links", {}).get("latest_draft")
    if not draft_url:
        raise RuntimeError("Zenodo did not return links.latest_draft for the new version")
    draft, _ = make_request(draft_url, method="GET", token=token, proxy=proxy)
    # Files are inherited from the previous version; remove them so only corrected artefacts remain.
    for f in draft.get("files", []):
        file_url = f.get("links", {}).get("self") or f"{base_url}/deposit/depositions/{draft['id']}/files/{f['id']}"
        make_request(file_url, method="DELETE", token=token, proxy=proxy)
        print(f"[+] Removed inherited file '{f.get('filename')}'")
    print(f"[+] New-version draft ID: {draft['id']}")
    return draft["id"], draft.get("links", {}).get("bucket"), draft

def publish_deposit(base_url, dep_id, token, proxy=None):
    url = f"{base_url}/deposit/depositions/{dep_id}/actions/publish"
    print(f"[*] Publishing deposition {dep_id}...")
    res, status = make_request(url, method="POST", data={}, token=token, proxy=proxy)
    doi = res.get("doi")
    doi_url = res.get("doi_url")
    record_id = res.get("record_id")
    print("\n" + "=" * 70)
    print(f"🎉 ZENODO PUBLICATION SUCCESSFUL!")
    print(f"   Record ID:  {record_id}")
    print(f"   DOI:        {doi}")
    print(f"   DOI URL:    {doi_url}")
    print(f"   Record URL: https://zenodo.org/records/{record_id}")
    print("=" * 70 + "\n")
    return res

def main():
    parser = argparse.ArgumentParser(description="Publish LeanFlow K3xT2 artifacts to Zenodo.")
    parser.add_argument("--sandbox", action="store_true", help="Publish to Zenodo Sandbox instead of Production.")
    parser.add_argument("--dry-run", action="store_true", help="Verify artifacts and metadata without uploading.")
    parser.add_argument("--proxy", type=str, default=os.environ.get("HTTPS_PROXY"), help="HTTP/HTTPS proxy URL (e.g. socks5://localhost:1080 or http://proxy:8080)")
    parser.add_argument("--no-publish", action="store_true", help="Create draft and upload files, but do not trigger final publish action.")
    parser.add_argument("--new-version-of", type=int, metavar="RECORD_ID", help="Deposit as a new version of this published record (keeps the concept DOI lineage) instead of a new record.")
    args = parser.parse_args()

    token = get_token()
    if not token:
        print("[-] ERROR: Zenodo token not found in ZENODO_TOKEN or ~/.config/zenodo/token", file=sys.stderr)
        sys.exit(1)

    base_url = SANDBOX_BASE_URL if args.sandbox else ZENODO_BASE_URL
    target_name = "Zenodo Sandbox" if args.sandbox else "Zenodo Production (zenodo.org)"
    print(f"[*] Target Repository: {target_name}")

    metadata = load_metadata()
    print(f"[*] Title: {metadata.get('title')}")
    print(f"[*] Author: {metadata.get('creators', [{}])[0].get('name')}")

    # Files to upload
    files_to_upload = [
        BUNDLE_DIR / "T_duality_Alone.pdf",
        BUNDLE_DIR / "T_duality_Alone.tex",
        BUNDLE_DIR / "lean4_formal_proofs.tar.gz",
        BUNDLE_DIR / "rust_stiff_dualscale_simulator.tar.gz",
        BUNDLE_DIR / "tda_mapper_skeleton.json",
        BUNDLE_DIR / "kummer_langevin_summary.json",
        BUNDLE_DIR / "vacuum_decay_cdl_summary.json",
        BUNDLE_DIR / "tachyon_condensation_summary.json",
        # The evidence record for the 2026-09-21 corrections. A deposit that
        # corrects claims should ship the document that establishes them.
        BUNDLE_DIR / "STREAM1_BRIDGE.md",
        BUNDLE_DIR / "PRE_REGISTRATION.md",
        BUNDLE_DIR / "TELEMETRY_PROVENANCE.json",
        WORKSPACE_DIR / "zenodo_deposit_bundle.zip"
    ]

    print("\n[*] Validating artifact files:")
    for fp in files_to_upload:
        if not fp.exists():
            print(f"[-] Missing artifact file: {fp}", file=sys.stderr)
            sys.exit(1)
        size_kb = fp.stat().st_size / 1024
        print(f"  ✓ {fp.name:35s} ({size_kb:8.1f} KB)")

    if args.dry_run:
        print("\n[+] Dry run complete. All files and metadata verified.")
        sys.exit(0)

    # 1. Create Deposition (or a new version of an existing record)
    if args.new_version_of:
        dep_id, bucket_url, dep_data = create_new_version(base_url, args.new_version_of, token, proxy=args.proxy)
    else:
        dep_id, bucket_url, dep_data = create_deposit(base_url, token, proxy=args.proxy)

    # 2. Upload Files
    for fp in files_to_upload:
        upload_file(base_url, dep_id, bucket_url, str(fp), token, proxy=args.proxy)

    # 3. Update Metadata
    set_metadata(base_url, dep_id, metadata, token, proxy=args.proxy)

    # 4. Publish
    if args.no_publish:
        print(f"[+] Draft deposit ready at: https://zenodo.org/deposit/{dep_id}")
    else:
        publish_deposit(base_url, dep_id, token, proxy=args.proxy)

if __name__ == "__main__":
    main()
