#!/usr/bin/env python3
"""
Rule-8 Cross-Consistency Checker for DualScaleSimulator
Verifies that quantities in simulation results, JSON summaries, and LaTeX papers
are consistent across all sources.
"""

import json
import re
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Any

def load_json(path: str) -> Any:
    """Load JSON file."""
    with open(path, 'r') as f:
        return json.load(f)

def get_canonical_value(params: Dict[str, Any], quantity: str) -> Tuple[Any, str]:
    """Extract canonical value from canonical source."""
    q = params["quantities"][quantity]
    source = q["canonical_source"]

    if source["type"] == "json_file":
        data = load_json(source["path"])
        # Navigate nested keys
        keys = source["key_path"].split(".")
        value = data
        for key in keys:
            value = value[key]
        return value, "json"
    elif source["type"] == "latex_constant":
        # For constants defined in LaTeX, use the provided value
        return q["canonical_value"], "latex"
    elif source["type"] == "computed":
        # For computed values, use the expected value
        return q["canonical_value"], "computed"
    else:
        raise ValueError(f"Unknown source type: {source['type']}")

def grep_occurrences(file_path: str, regex: str) -> List[str]:
    """Find all occurrences of a pattern in a file."""
    try:
        with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
    except FileNotFoundError:
        return []

    matches = re.findall(regex, content)
    return matches

def parse_value(val_str: str, comparison_type: str) -> Any:
    """Parse a string value to the appropriate type."""
    val_str = val_str.strip()

    # Try to parse as scientific notation (e.g., "1.24", "1e-4")
    if 'e' in val_str.lower() or '.' in val_str:
        try:
            return float(val_str)
        except ValueError:
            pass

    # Try to parse as integer
    try:
        return int(val_str)
    except ValueError:
        pass

    # Return as string
    return val_str

def compare_values(canonical: Any, found: Any, comparison_type: str, tolerance: float = 1e-9) -> Tuple[bool, str]:
    """
    Compare two values according to the comparison type.
    Returns (match, reason).
    """
    if comparison_type == "exact_integer":
        match = int(canonical) == int(found)
        reason = f"{canonical} vs {found}"
        return match, reason
    elif comparison_type == "relative_tolerance":
        try:
            c_float = float(canonical)
            f_float = float(found)
            if c_float == 0:
                match = abs(f_float - c_float) < tolerance
            else:
                rel_diff = abs(f_float - c_float) / abs(c_float)
                match = rel_diff < tolerance
            reason = f"{c_float} vs {f_float} (rel_diff={abs(f_float - c_float) / abs(c_float) if c_float != 0 else 'N/A'})"
            return match, reason
        except (ValueError, TypeError):
            return False, f"Could not parse as float: {canonical} vs {found}"
    else:
        return False, f"Unknown comparison type: {comparison_type}"

def main():
    # Load parameters
    params = load_json("audit/parameters.json")

    results = []
    mismatches = []

    # Check each quantity
    for quantity, q_data in params["quantities"].items():
        canonical_value, source_type = get_canonical_value(params, quantity)
        comparison_type = q_data["comparison"]
        tolerance = q_data.get("tolerance", 1e-9)

        # Build result row
        row = {
            "quantity": quantity,
            "canonical": canonical_value,
            "source_type": source_type,
            "occurrences": [],
            "mismatches": []
        }

        # Check occurrences
        for occurrence in q_data.get("occurrences", []):
            file_path = occurrence["file"]
            regex = occurrence["regex"]

            matches = grep_occurrences(file_path, regex)

            if not matches:
                row["occurrences"].append({
                    "file": file_path,
                    "status": "NOT_FOUND",
                    "values": []
                })
            else:
                for match_value in matches:
                    parsed = parse_value(match_value, comparison_type)
                    match, reason = compare_values(canonical_value, parsed, comparison_type, tolerance)

                    if match:
                        status = "OK"
                    else:
                        status = "MISMATCH"
                        row["mismatches"].append({
                            "file": file_path,
                            "found": parsed,
                            "reason": reason
                        })
                        mismatches.append({
                            "quantity": quantity,
                            "canonical": canonical_value,
                            "file": file_path,
                            "found": parsed
                        })

                    row["occurrences"].append({
                        "file": file_path,
                        "status": status,
                        "value": parsed,
                        "reason": reason
                    })

        results.append(row)

    # Generate markdown report
    report_lines = [
        "# Cross-Consistency Check Report",
        "",
        "## Summary",
        f"- Total mismatches: {len(mismatches)}",
        ""
    ]

    if mismatches:
        report_lines.append("## Mismatches Detected")
        report_lines.append("")
        report_lines.append("| Quantity | Canonical Value | Location | Found Value |")
        report_lines.append("|----------|-----------------|----------|-------------|")

        for mismatch in mismatches:
            report_lines.append(
                f"| {mismatch['quantity']} | {mismatch['canonical']} | "
                f"{mismatch['file']} | {mismatch['found']} |"
            )

    # Full table
    report_lines.append("")
    report_lines.append("## Detailed Results")
    report_lines.append("")
    report_lines.append("| Quantity | Canonical | Type | Status |")
    report_lines.append("|----------|-----------|------|--------|")

    for row in results:
        if not row["occurrences"]:
            # No occurrences to check
            status_str = "NO_CHECKS"
        else:
            statuses = [occ["status"] for occ in row["occurrences"]]
            if any(s == "MISMATCH" for s in statuses):
                status_str = "MISMATCH"
            elif any(s == "NOT_FOUND" for s in statuses):
                status_str = "NOT_FOUND"
            else:
                status_str = "OK"

        report_lines.append(
            f"| {row['quantity']} | {row['canonical']} | {row['source_type']} | {status_str} |"
        )

    # Write report
    report_content = "\n".join(report_lines)
    with open("audit/cross_consistency.md", "w") as f:
        f.write(report_content)

    print(report_content)

    # Exit code: 1 if any mismatches found, 0 otherwise
    sys.exit(1 if mismatches else 0)

if __name__ == "__main__":
    main()
