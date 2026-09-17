"""
pytest configuration and fixtures for WorkshopCosmo test suite.
"""

import os
import pytest


@pytest.fixture(autouse=True)
def setup_output_dir(tmp_path, monkeypatch):
    """
    Autouse fixture that sets OUTPUT_DIR env var to tmp_path for all tests.
    This ensures all test outputs go to temp directory and don't modify committed files.
    Also seeds the TDA Mapper for deterministic output.
    """
    monkeypatch.setenv("OUTPUT_DIR", str(tmp_path))
    monkeypatch.setenv("MAPPER_SEED", "42")
    yield
