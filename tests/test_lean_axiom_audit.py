#!/usr/bin/env python3
"""
tests/test_lean_axiom_audit.py — Testing lean_axiom_audit.py gate (G3).

Three core controls:
  1. Positive: Scratch theorem using `rfl` → must be STANDARD, exit 0
  2. Negative: Scratch file declaring an axiom and using it → must be CUSTOM_AXIOM, exit 1
  3. Malformed: Bad Lean syntax → must exit 2 (elaboration error)

Skip if lake is not available (required for audit).
"""

import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

# Determine the project root
TEST_DIR = Path(__file__).parent
ROOT = TEST_DIR.parent
AUDIT_SCRIPT = ROOT / "scripts" / "lean_axiom_audit.py"


@pytest.mark.skipif(shutil.which("lake") is None, reason="lake not found in PATH")
class TestLeanAxiomAudit:
    """Gate G3: axiom audit tests."""

    def test_positive_control_standard_theorem(self):
        """Positive control: theorem using rfl should be STANDARD, exit 0."""
        scratch_content = """theorem test_one_plus_one : 1 + 1 = 2 := rfl
"""
        returncode, result = run_scratch_audit(scratch_content)
        assert returncode == 0, f"Expected exit 0, got {returncode}"
        assert result["status"] == "OK", f"Expected OK, got {result['status']}: {result.get('error')}"
        assert "test_one_plus_one" in result["theorems"], "Theorem not found in results"
        theorem_data = result["theorems"]["test_one_plus_one"]
        assert theorem_data["classification"] == "STANDARD", \
            f"Expected STANDARD, got {theorem_data['classification']}: {theorem_data['axioms']}"
        print(f"✓ Positive control passed: exit {returncode}, {theorem_data['classification']}")

    def test_negative_control_custom_axiom(self):
        """Negative control: theorem depending on custom axiom should be CUSTOM_AXIOM, exit 1."""
        scratch_content = """axiom my_custom_axiom : True

theorem uses_custom_axiom : True := my_custom_axiom
"""
        returncode, result = run_scratch_audit(scratch_content)
        assert returncode == 1, f"Expected exit 1 for CUSTOM_AXIOM, got {returncode}"
        # Status should be FAIL when custom axioms are found
        assert result["status"] == "FAIL", \
            f"Expected FAIL, got {result['status']}: {result.get('error')}"
        assert "uses_custom_axiom" in result["theorems"], "Theorem not found in results"
        theorem_data = result["theorems"]["uses_custom_axiom"]
        assert theorem_data["classification"] == "CUSTOM_AXIOM", \
            f"Expected CUSTOM_AXIOM, got {theorem_data['classification']}: {theorem_data['axioms']}"
        print(f"✓ Negative control passed: exit {returncode}, {theorem_data['classification']}")

    def test_malformed_lean_elaboration_error(self):
        """Malformed Lean in a theorem should result in elaboration error, exit 2."""
        scratch_content = """theorem test_bad : NonexistentType := sorry
"""
        returncode, result = run_scratch_audit(scratch_content)
        assert returncode == 2, f"Expected exit 2 for elaboration error, got {returncode}"
        # Elaboration errors should give status ERROR
        assert result["status"] == "ERROR", \
            f"Expected ERROR for malformed Lean, got {result['status']}: {result.get('error')}"
        assert result.get("error") is not None, "Expected error message"
        print(f"✓ Malformed Lean test passed: exit {returncode}, {result['error'][:60]}")


def run_scratch_audit(lean_code: str) -> tuple[int, dict]:
    """
    Write lean_code to a temporary scratch file and audit it.
    Returns (returncode, result_dict) where result_dict is the parsed JSON.
    """
    with tempfile.NamedTemporaryFile("w", suffix=".lean", delete=False) as fh:
        fh.write(lean_code)
        scratch_path = fh.name

    try:
        proc = subprocess.run(
            [sys.executable, str(AUDIT_SCRIPT), "--scratch", scratch_path],
            cwd=ROOT,
            capture_output=True,
            text=True,
            timeout=30
        )

        try:
            result = json.loads(proc.stdout)
        except json.JSONDecodeError:
            result = {
                "status": "ERROR",
                "error": f"Failed to parse audit output: {proc.stdout[:200]}",
                "theorems": {}
            }

        return proc.returncode, result
    finally:
        Path(scratch_path).unlink(missing_ok=True)


if __name__ == "__main__":
    # Run tests manually if executed directly
    if shutil.which("lake") is None:
        print("⚠️  lake not found in PATH; skipping tests")
        sys.exit(0)

    test = TestLeanAxiomAudit()

    try:
        print("Running positive control (theorem with rfl)...")
        test.test_positive_control_standard_theorem()

        print("\nRunning negative control (custom axiom)...")
        test.test_negative_control_custom_axiom()

        print("\nRunning malformed Lean test...")
        test.test_malformed_lean_elaboration_error()

        print("\n✅ All tests passed!")
        sys.exit(0)
    except AssertionError as e:
        print(f"\n❌ Test failed: {e}")
        sys.exit(1)
