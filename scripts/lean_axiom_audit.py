#!/usr/bin/env python3
"""
scripts/lean_axiom_audit.py — gate G3 enforcement for DualScaleSimulator.

Adapted from /home/callensxavier_gmail_com/SocrateAI-Scientific-Agora-LeanMaster/tools/axiom_audit.py
and includes inlined strip_comments from /home/callensxavier_gmail_com/SocrateAI-Scientific-Agora-LeanMaster/tools/statement_lock.py

For every `theorem`/`lemma` declared in the project's libraries (via lakefile.lean), run
`#print axioms` against the *built* library and classify by axiom set. A theorem passes if it
depends only on Lean's three standard axioms: propext, Classical.choice, Quot.sound.

Exit codes:
  0: All libraries verified; all theorems are STANDARD (propext/Classical.choice/Quot.sound only)
  1: Some library or theorem failed verification (CUSTOM_AXIOM, SORRY, NATIVE, or MISSING found)
  2: Audit infrastructure error (lake/Lean missing, library not built, zero libraries discovered)

Usage:
  lake build                               # Build all libraries
  python3 scripts/lean_axiom_audit.py      # Audit all libraries
  python3 scripts/lean_axiom_audit.py MathieuVertexOperators     # Single library
  python3 scripts/lean_axiom_audit.py --scratch <file.lean>      # Test mode on scratch file

Output: audit/lean_axiom_report.json (per-theorem classification)
"""

import os
import re
import subprocess
import sys
import tempfile
import json
from pathlib import Path

# ============================================================================
# CONFIGURATION
# ============================================================================

ROOT = Path(os.environ.get("LEAN_PROJECT_ROOT") or Path(__file__).resolve().parent.parent).resolve()
STANDARD = {"propext", "Classical.choice", "Quot.sound"}
DECL = re.compile(r"^(?:theorem|lemma)\s+([A-Za-z0-9_'.]+)", re.MULTILINE)
NAMESPACE = re.compile(r"^(namespace|end)\s+([A-Za-z0-9_.]+)\s*$", re.MULTILINE)
LAKEFILE = ROOT / "lakefile.lean"


# ============================================================================
# UTILITY: strip_comments (inlined from statement_lock.py)
# ============================================================================

def strip_comments(text: str) -> str:
    """Remove Lean comments, keeping the text length-independent of documentation.
    Block comments nest in Lean, so this is a small scanner rather than a regex.

    Inlined from /home/callensxavier_gmail_com/SocrateAI-Scientific-Agora-LeanMaster/tools/statement_lock.py
    lines 46-66.
    """
    out, i, depth, n = [], 0, 0, len(text)
    while i < n:
        two = text[i:i + 2]
        if two == "/-":
            depth += 1
            i += 2
        elif two == "-/" and depth:
            depth -= 1
            i += 2
        elif depth:
            i += 1
        elif two == "--":
            j = text.find("\n", i)
            i = n if j < 0 else j
        else:
            out.append(text[i])
            i += 1
    return "".join(out)


# ============================================================================
# LIBRARY DISCOVERY FROM LAKEFILE
# ============================================================================

def discover_libraries() -> dict[str, str]:
    """Parse lakefile.lean and extract (library_name -> source_file_path) mappings.

    Each lean_lib entry with srcDir="proofs" generates a mapping:
      - MathieuVertexOperators -> proofs/MathieuVertexOperators.lean
      - LeanscratchDB.HoloAlg -> proofs/LeanscratchDB/HoloAlg.lean

    Returns dict{lib_name: str -> proofs_file_path: str}
    Raises ValueError if zero libraries found (catastrophic infrastructure failure).
    """
    if not LAKEFILE.exists():
        raise ValueError(f"lakefile.lean not found at {LAKEFILE}")

    lakefile_text = LAKEFILE.read_text(encoding="utf-8")
    lib_pattern = re.compile(
        r'lean_lib\s+([A-Za-z0-9_.]+)\s+\{\s*srcDir\s*:=\s*"([^"]+)"\s*\}',
        re.MULTILINE
    )

    libraries = {}
    for match in lib_pattern.finditer(lakefile_text):
        lib_name = match.group(1)
        src_dir = match.group(2)

        # Map dotted name to file path: LeanscratchDB.HoloAlg -> LeanscratchDB/HoloAlg.lean
        parts = lib_name.split(".")

        # For nested namespaces, build the full path: LeanscratchDB.HoloAlg -> proofs/LeanscratchDB/HoloAlg.lean
        if len(parts) > 1:
            rel_path = "/".join(parts) + ".lean"
            file_path = Path(src_dir) / rel_path
        else:
            file_path = Path(src_dir) / (parts[-1] + ".lean")

        libraries[lib_name] = str(file_path)

    if not libraries:
        raise ValueError(
            f"No lean_lib entries with srcDir found in {LAKEFILE}. "
            "Catastrophic: cannot proceed."
        )

    return libraries


# ============================================================================
# THEOREM/LEMMA EXTRACTION
# ============================================================================

def qualified_names(path: Path) -> list[str]:
    """Extract all theorem/lemma names in qualified form (with namespaces).

    Removes comments first to avoid mistaking docstring lines for declarations.
    """
    text = strip_comments(path.read_text(encoding="utf-8"))
    events = sorted(
        [(m.start(), "ns", m.group(1), m.group(2)) for m in NAMESPACE.finditer(text)]
        + [(m.start(), "decl", m.group(1), None) for m in DECL.finditer(text)]
    )
    stack: list[str] = []
    names = []
    for _, kind, a, b in events:
        if kind == "ns":
            if a == "namespace":
                stack.append(b)
            elif stack and stack[-1] == b:
                stack.pop()
        else:
            names.append(".".join(stack + [a]))
    return names


# ============================================================================
# AXIOM CLASSIFICATION
# ============================================================================

def classify_axioms(axiom_set: set[str]) -> str:
    """Classify a set of axioms into one of five categories.

    Returns:
      "STANDARD"      - Only propext, Classical.choice, Quot.sound
      "CUSTOM_AXIOM"  - Contains user-declared axioms (not in STANDARD)
      "SORRY"         - Contains sorryAx (a remaining sorry)
      "NATIVE"        - Contains Lean.ofReduceBool (native_decide)
      "MISSING"       - No axiom info returned (name resolution failed)
    """
    if not axiom_set:
        return "STANDARD"

    extra = axiom_set - STANDARD

    if "sorryAx" in extra:
        return "SORRY"
    if "Lean.ofReduceBool" in extra:
        return "NATIVE"
    if extra:
        return "CUSTOM_AXIOM"

    return "STANDARD"


# ============================================================================
# AUDIT: RUN #print axioms FOR SINGLE LIBRARY
# ============================================================================

def audit_library(lib_name: str, source_file: str, scratch_mode: bool = False) -> dict[str, any]:
    """Audit one library by running #print axioms for each theorem/lemma.

    Returns:
      {
        "library": str,
        "file": str,
        "status": "OK" | "FAIL" | "ERROR",
        "error": str | None,
        "theorems": {
          "name": {
            "axioms": [str],
            "classification": "STANDARD" | "CUSTOM_AXIOM" | "SORRY" | "NATIVE" | "MISSING"
          },
          ...
        },
        "summary": {
          "total": int,
          "standard": int,
          "custom_axiom": int,
          "sorry": int,
          "native": int,
          "missing": int
        }
      }
    """
    result = {
        "library": lib_name,
        "file": source_file,
        "status": None,
        "error": None,
        "theorems": {},
        "summary": {
            "total": 0,
            "standard": 0,
            "custom_axiom": 0,
            "sorry": 0,
            "native": 0,
            "missing": 0
        }
    }

    source_path = ROOT / source_file
    if not source_path.exists():
        result["status"] = "ERROR"
        result["error"] = f"Source file not found: {source_file}"
        return result

    names = qualified_names(source_path)
    if not names:
        # Empty library is not an error, but is suspicious — report it
        result["status"] = "OK"
        result["theorems"] = {}
        return result

    # Build probe file
    if scratch_mode:
        # For scratch tests: include the original source plus #print axioms lines
        source_text = source_path.read_text(encoding="utf-8")
        probe = source_text + "\n" + "".join(f"#print axioms {n}\n" for n in names)
    else:
        # For library audit: import the library, then #print axioms
        imports = f"import {lib_name}\n"
        probe = imports + "".join(f"#print axioms {n}\n" for n in names)

    with tempfile.NamedTemporaryFile("w", suffix=".lean", delete=False) as fh:
        fh.write(probe)
        tmp = fh.name

    try:
        proc = subprocess.run(
            ["lake", "env", "lean", "-DmaxHeartbeats=1000000", tmp],
            cwd=ROOT,
            capture_output=True,
            text=True,
            timeout=120
        )

        out = proc.stdout
        combined = proc.stdout + proc.stderr

        # Check for elaboration errors
        import_errors = [ln for ln in combined.splitlines() if ": error" in ln]
        if import_errors:
            result["status"] = "ERROR"
            result["error"] = f"Probe elaboration failed (library not built or import failed): {import_errors[0]}"
            return result

        # Parse #print axioms output
        blocks = re.split(r"(?='[^']+' (?:depends on axioms|does not depend on any axioms))", out)
        results = {}
        for blk in blocks:
            m = re.match(
                r"'([^']+)' (depends on axioms: \[([^\]]*)\]|does not depend on any axioms)",
                blk.strip(),
                re.S
            )
            if m:
                axiom_str = m.group(3) or ""
                axioms = {a.strip() for a in axiom_str.replace("\n", " ").split(",") if a.strip()}
                results[m.group(1)] = axioms

        # Classify each theorem
        failed_count = 0
        for n in names:
            if n not in results:
                classification = "MISSING"
                axioms_list = []
                failed_count += 1
            else:
                axioms_list = sorted(results[n])
                classification = classify_axioms(results[n])
                if classification != "STANDARD":
                    failed_count += 1

            result["theorems"][n] = {
                "axioms": axioms_list,
                "classification": classification
            }
            result["summary"][classification.lower()] += 1

        result["summary"]["total"] = len(names)
        result["status"] = "FAIL" if failed_count else "OK"
        return result

    finally:
        Path(tmp).unlink(missing_ok=True)


# ============================================================================
# SCRATCH FILE MODE (for testing)
# ============================================================================

def audit_scratch_file(scratch_path: str) -> dict[str, any]:
    """Audit a scratch Lean file directly (for testing).

    The file should contain theorems/lemmas but no imports (stays outside project).
    """
    scratch = Path(scratch_path)
    if not scratch.exists():
        return {
            "status": "ERROR",
            "error": f"Scratch file not found: {scratch_path}",
            "theorems": {},
            "summary": {}
        }

    # This is a scratch file audit — don't import anything
    return audit_library(scratch.stem, scratch_path, scratch_mode=True)


# ============================================================================
# MAIN: ORCHESTRATE AUDITS FOR ALL LIBRARIES
# ============================================================================

def main() -> int:
    """Main entry point: audit all or specified libraries."""

    args = sys.argv[1:]

    # Check for scratch mode
    if args and args[0] == "--scratch":
        if len(args) < 2:
            print("Usage: --scratch <file.lean>")
            return 2
        result = audit_scratch_file(args[1])
        print(json.dumps(result, indent=2))
        status = result.get("status")
        return 1 if status == "FAIL" else (2 if status == "ERROR" else 0)

    # Discover libraries
    try:
        libraries = discover_libraries()
    except ValueError as e:
        print(f"AUDIT ERROR — {e}")
        return 2

    if args:
        # Audit specific libraries
        target_libs = {name: libraries[name] for name in args if name in libraries}
        if not target_libs:
            print(f"AUDIT ERROR — None of the specified libraries found: {args}")
            return 2
    else:
        # Audit all
        target_libs = libraries

    # Run audit on each library
    all_results = {}
    all_failed = False

    for lib_name in sorted(target_libs.keys()):
        source_file = target_libs[lib_name]
        result = audit_library(lib_name, source_file)
        all_results[lib_name] = result

        if result["status"] == "ERROR":
            all_failed = True
            print(f"ERROR   {lib_name}  {result['error']}")
        elif result["status"] == "FAIL":
            all_failed = True
            # Print summary for failed library
            summary = result["summary"]
            standard_ok = summary.get("standard", 0) == summary.get("total", 1)
            if not standard_ok:
                print(f"FAIL    {lib_name}  {summary['standard']}/{summary['total']} theorems standard")
        else:
            summary = result["summary"]
            print(f"OK      {lib_name}  {summary.get('total', 0)} theorems verified")

    # Print detailed failures
    print()
    for lib_name in sorted(all_results.keys()):
        result = all_results[lib_name]
        if result["status"] == "FAIL":
            print(f"\n{lib_name}:")
            for theo_name, theo_data in result["theorems"].items():
                classification = theo_data["classification"]
                if classification != "STANDARD":
                    axioms_str = ", ".join(theo_data["axioms"]) if theo_data["axioms"] else "[]"
                    print(f"  {classification:15} {theo_name:40} {axioms_str}")

    # Write JSON report
    report_path = ROOT / "audit" / "lean_axiom_report.json"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(all_results, indent=2) + "\n")
    print(f"\nReport written to: {report_path}")

    # Summary
    total_count = sum(r["summary"].get("total", 0) for r in all_results.values())
    failed_count = sum(
        r["summary"].get("total", 0) - r["summary"].get("standard", 0)
        for r in all_results.values()
        if r["status"] != "ERROR"
    )
    error_count = sum(1 for r in all_results.values() if r["status"] == "ERROR")

    print(f"\nAudit summary: {len(all_results)} libraries, {total_count} theorems audited")
    if error_count:
        print(f"  {error_count} library(ies) errored (not built or missing)")
    if failed_count:
        print(f"  {failed_count} theorem(s) with non-standard axioms")
    else:
        print(f"  All theorems are STANDARD ✓")

    # Determine exit code: FAIL wins over ERROR
    # 0 = all standard, 1 = some failure (CUSTOM_AXIOM/SORRY/NATIVE/MISSING), 2 = only infrastructure errors
    if failed_count:
        return 1
    if error_count:
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
