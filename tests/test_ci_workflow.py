"""
=============================================================================
CI Workflow Validation Test Suite
=============================================================================
Validates that the .github/workflows/ci.yml file conforms to the specification:
- Four jobs: rust, python, lean, consistency
- Correct job configurations and step structures
- Proper error handling and artifact uploads

Tests include positive controls (valid YAML passes) and negative controls
(mutated YAML correctly fails validation).
=============================================================================
"""

import os
import json
import pytest
import yaml
import copy


REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
CI_WORKFLOW_PATH = os.path.join(REPO_ROOT, ".github", "workflows", "ci.yml")


def load_ci_workflow():
    """Load and parse the CI workflow YAML."""
    with open(CI_WORKFLOW_PATH, 'r') as f:
        return yaml.safe_load(f)


def validate_ci_workflow(workflow_dict):
    """
    Validate CI workflow structure against spec.
    Raises AssertionError if validation fails.
    """
    # Check required top-level structure
    assert isinstance(workflow_dict, dict), "Workflow must be a dict"
    assert 'jobs' in workflow_dict, "Workflow missing 'jobs' key"
    jobs = workflow_dict['jobs']

    # Check required job names
    required_jobs = ['rust', 'python', 'lean', 'consistency']
    for job_name in required_jobs:
        assert job_name in jobs, f"Missing required job: {job_name}"

    # Validate rust job
    rust_job = jobs['rust']
    assert 'steps' in rust_job, "Rust job missing 'steps'"
    rust_steps = rust_job['steps']
    test_step = None
    for step in rust_steps:
        if 'run' in step and 'cargo test --release' in step['run']:
            test_step = step
            break
    assert test_step is not None, "Rust job missing 'cargo test --release' step"
    assert test_step.get('working-directory') == 'rust_simulator', \
        "Rust test step must have working-directory: rust_simulator"

    # Validate python job
    python_job = jobs['python']
    assert 'steps' in python_job, "Python job missing 'steps'"
    python_steps = python_job['steps']
    pytest_step = None
    for step in python_steps:
        if 'run' in step and 'pytest -q' in step['run']:
            pytest_step = step
            break
    assert pytest_step is not None, "Python job missing 'pytest -q' step"

    # Validate lean job
    lean_job = jobs['lean']
    assert 'steps' in lean_job, "Lean job missing 'steps'"
    lean_steps = lean_job['steps']

    # Check elan PATH configuration
    elan_step = None
    for step in lean_steps:
        if 'Install elan' in step.get('name', ''):
            elan_step = step
            break
    assert elan_step is not None, "Lean job missing 'Install elan' step"
    assert '"$HOME/.elan/bin" >> $GITHUB_PATH' in elan_step.get('run', ''), \
        "Lean job must add $HOME/.elan/bin to PATH"

    # Check that no step sources a non-existent activate script
    for step in lean_steps:
        run_cmd = step.get('run', '')
        assert 'source "$HOME/.elan/toolchains/' not in run_cmd, \
            "Lean job must not source activate scripts from toolchains directory"

    # Verify lake build step exists and is simple
    lake_step = None
    for step in lean_steps:
        if 'run' in step and 'lake build' in step['run']:
            lake_step = step
            break
    assert lake_step is not None, "Lean job missing 'lake build' step"

    # Check lean axiom audit step
    audit_step = None
    for step in lean_steps:
        if 'Run lean axiom audit' in step.get('name', ''):
            audit_step = step
            break
    assert audit_step is not None, "Lean job missing 'Run lean axiom audit' step"
    assert audit_step.get('continue-on-error') is True, \
        "Lean axiom audit step must have continue-on-error: true"

    # Check lean audit artifact upload
    audit_upload_step = None
    for step in lean_steps:
        if step.get('uses', '').startswith('actions/upload-artifact'):
            if step.get('with', {}).get('name') == 'lean_axiom_report':
                audit_upload_step = step
                break
    assert audit_upload_step is not None, \
        "Lean job missing upload for lean_axiom_report artifact"
    assert audit_upload_step.get('with', {}).get('path') == 'audit/lean_axiom_report.json', \
        "Lean audit upload must have path: audit/lean_axiom_report.json"

    # Validate consistency job
    consistency_job = jobs['consistency']
    assert 'steps' in consistency_job, "Consistency job missing 'steps'"
    consistency_steps = consistency_job['steps']

    # Check consistency check step
    check_step = None
    for step in consistency_steps:
        if 'Run cross consistency check' in step.get('name', ''):
            check_step = step
            break
    assert check_step is not None, "Consistency job missing 'Run cross consistency check' step"
    assert check_step.get('continue-on-error') is True, \
        "Consistency check step must have continue-on-error: true"

    # Check consistency artifact upload
    consistency_upload_step = None
    for step in consistency_steps:
        if step.get('uses', '').startswith('actions/upload-artifact'):
            if step.get('with', {}).get('name') == 'cross_consistency_report':
                consistency_upload_step = step
                break
    assert consistency_upload_step is not None, \
        "Consistency job missing upload for cross_consistency artifact"
    assert consistency_upload_step.get('with', {}).get('path') == 'audit/cross_consistency.md', \
        "Consistency upload must have path: audit/cross_consistency.md"


class TestCIWorkflowPositiveControl:
    """Positive control: the actual CI workflow should pass validation."""

    def test_ci_workflow_file_exists(self):
        """Verify CI workflow file exists."""
        assert os.path.exists(CI_WORKFLOW_PATH), \
            f"CI workflow file not found at {CI_WORKFLOW_PATH}"

    def test_ci_workflow_is_valid_yaml(self):
        """Verify CI workflow is syntactically valid YAML."""
        workflow = load_ci_workflow()
        assert workflow is not None, "CI workflow failed to parse as YAML"

    def test_ci_workflow_has_required_jobs(self):
        """Verify CI workflow has all required jobs."""
        workflow = load_ci_workflow()
        jobs = workflow.get('jobs', {})
        required_jobs = ['rust', 'python', 'lean', 'consistency']
        for job_name in required_jobs:
            assert job_name in jobs, f"Missing required job: {job_name}"

    def test_ci_workflow_passes_full_validation(self):
        """Positive control: valid CI workflow passes all validation rules."""
        workflow = load_ci_workflow()
        validate_ci_workflow(workflow)  # Should not raise AssertionError


class TestCIWorkflowNegativeControl:
    """Negative controls: mutations of the CI workflow should fail validation."""

    def test_missing_rust_job_fails_validation(self):
        """Negative control: CI workflow without rust job should fail."""
        workflow = load_ci_workflow()
        workflow_mutated = copy.deepcopy(workflow)
        del workflow_mutated['jobs']['rust']

        with pytest.raises(AssertionError, match="Missing required job: rust"):
            validate_ci_workflow(workflow_mutated)

    def test_missing_python_job_fails_validation(self):
        """Negative control: CI workflow without python job should fail."""
        workflow = load_ci_workflow()
        workflow_mutated = copy.deepcopy(workflow)
        del workflow_mutated['jobs']['python']

        with pytest.raises(AssertionError, match="Missing required job: python"):
            validate_ci_workflow(workflow_mutated)

    def test_missing_lean_job_fails_validation(self):
        """Negative control: CI workflow without lean job should fail."""
        workflow = load_ci_workflow()
        workflow_mutated = copy.deepcopy(workflow)
        del workflow_mutated['jobs']['lean']

        with pytest.raises(AssertionError, match="Missing required job: lean"):
            validate_ci_workflow(workflow_mutated)

    def test_missing_consistency_job_fails_validation(self):
        """Negative control: CI workflow without consistency job should fail."""
        workflow = load_ci_workflow()
        workflow_mutated = copy.deepcopy(workflow)
        del workflow_mutated['jobs']['consistency']

        with pytest.raises(AssertionError, match="Missing required job: consistency"):
            validate_ci_workflow(workflow_mutated)

    def test_rust_missing_cargo_test_fails_validation(self):
        """Negative control: rust job without cargo test --release should fail."""
        workflow = load_ci_workflow()
        workflow_mutated = copy.deepcopy(workflow)
        # Remove cargo test step
        rust_steps = workflow_mutated['jobs']['rust']['steps']
        workflow_mutated['jobs']['rust']['steps'] = [
            step for step in rust_steps
            if not ('run' in step and 'cargo test' in step['run'])
        ]

        with pytest.raises(AssertionError, match="cargo test --release"):
            validate_ci_workflow(workflow_mutated)

    def test_rust_wrong_working_directory_fails_validation(self):
        """Negative control: rust test without rust_simulator working-directory should fail."""
        workflow = load_ci_workflow()
        workflow_mutated = copy.deepcopy(workflow)
        # Find and modify the cargo test step
        for step in workflow_mutated['jobs']['rust']['steps']:
            if 'run' in step and 'cargo test --release' in step['run']:
                step['working-directory'] = 'wrong_dir'
                break

        with pytest.raises(AssertionError, match="rust_simulator"):
            validate_ci_workflow(workflow_mutated)

    def test_python_missing_pytest_fails_validation(self):
        """Negative control: python job without pytest -q should fail."""
        workflow = load_ci_workflow()
        workflow_mutated = copy.deepcopy(workflow)
        # Remove pytest step
        python_steps = workflow_mutated['jobs']['python']['steps']
        workflow_mutated['jobs']['python']['steps'] = [
            step for step in python_steps
            if not ('run' in step and 'pytest -q' in step['run'])
        ]

        with pytest.raises(AssertionError, match="pytest -q"):
            validate_ci_workflow(workflow_mutated)

    def test_lean_wrong_elan_path_fails_validation(self):
        """Negative control: lean job with wrong elan PATH should fail."""
        workflow = load_ci_workflow()
        workflow_mutated = copy.deepcopy(workflow)
        # Find and modify the elan install step
        for step in workflow_mutated['jobs']['lean']['steps']:
            if 'Install elan' in step.get('name', ''):
                step['run'] = 'curl https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh -sSf | sh -s -- -y --default-toolchain none\necho "wrong_path" >> $GITHUB_PATH'
                break

        with pytest.raises(AssertionError, match="\\$HOME/.elan/bin"):
            validate_ci_workflow(workflow_mutated)

    def test_lean_with_activate_source_fails_validation(self):
        """Negative control: lean job with activate script sourcing should fail."""
        workflow = load_ci_workflow()
        workflow_mutated = copy.deepcopy(workflow)
        # Find and modify the lake build step to add activate sourcing
        for step in workflow_mutated['jobs']['lean']['steps']:
            if 'run' in step and 'lake build' in step['run']:
                step['run'] = 'source "$HOME/.elan/toolchains/lean4/activate"\nlake build'
                break

        with pytest.raises(AssertionError, match="activate scripts"):
            validate_ci_workflow(workflow_mutated)

    def test_lean_audit_missing_continue_on_error_fails_validation(self):
        """Negative control: lean axiom audit without continue-on-error should fail."""
        workflow = load_ci_workflow()
        workflow_mutated = copy.deepcopy(workflow)
        # Find and modify the audit step
        for step in workflow_mutated['jobs']['lean']['steps']:
            if 'Run lean axiom audit' in step.get('name', ''):
                step['continue-on-error'] = False
                break

        with pytest.raises(AssertionError, match="continue-on-error: true"):
            validate_ci_workflow(workflow_mutated)

    def test_lean_audit_upload_wrong_path_fails_validation(self):
        """Negative control: lean audit upload with wrong path should fail."""
        workflow = load_ci_workflow()
        workflow_mutated = copy.deepcopy(workflow)
        # Find and modify the upload step
        for step in workflow_mutated['jobs']['lean']['steps']:
            if step.get('uses', '').startswith('actions/upload-artifact'):
                if step.get('with', {}).get('name') == 'lean_axiom_report':
                    step['with']['path'] = 'wrong/path.json'
                    break

        with pytest.raises(AssertionError, match="audit/lean_axiom_report.json"):
            validate_ci_workflow(workflow_mutated)

    def test_consistency_missing_continue_on_error_fails_validation(self):
        """Negative control: consistency check without continue-on-error should fail."""
        workflow = load_ci_workflow()
        workflow_mutated = copy.deepcopy(workflow)
        # Find and modify the consistency check step
        for step in workflow_mutated['jobs']['consistency']['steps']:
            if 'Run cross consistency check' in step.get('name', ''):
                step['continue-on-error'] = False
                break

        with pytest.raises(AssertionError, match="continue-on-error: true"):
            validate_ci_workflow(workflow_mutated)

    def test_consistency_upload_wrong_path_fails_validation(self):
        """Negative control: consistency upload with wrong path should fail."""
        workflow = load_ci_workflow()
        workflow_mutated = copy.deepcopy(workflow)
        # Find and modify the upload step
        for step in workflow_mutated['jobs']['consistency']['steps']:
            if step.get('uses', '').startswith('actions/upload-artifact'):
                if step.get('with', {}).get('name') == 'cross_consistency_report':
                    step['with']['path'] = 'wrong/path.md'
                    break

        with pytest.raises(AssertionError, match="audit/cross_consistency.md"):
            validate_ci_workflow(workflow_mutated)
