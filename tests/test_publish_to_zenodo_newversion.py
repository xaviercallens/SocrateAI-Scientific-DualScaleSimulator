import importlib.util
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parent.parent / "scripts" / "publish_to_zenodo.py"


@pytest.fixture
def zen():
    spec = importlib.util.spec_from_file_location("publish_to_zenodo", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_new_version_uses_newversion_action_and_clears_inherited_files(zen, monkeypatch):
    calls = []
    draft = {"id": 999, "links": {"bucket": "https://z/bucket/999"},
             "files": [{"id": "f1", "filename": "old.pdf", "links": {"self": "https://z/files/f1"}}]}

    def fake(url, method="GET", data=None, headers=None, token=None, proxy=None):
        calls.append((method, url))
        if url.endswith("/actions/newversion"):
            return {"links": {"latest_draft": "https://z/api/deposit/depositions/999"}}, 201
        if method == "GET":
            return draft, 200
        return {}, 204

    monkeypatch.setattr(zen, "make_request", fake)
    dep_id, bucket, _ = zen.create_new_version("https://z/api", 22683565, "tok")
    assert (dep_id, bucket) == (999, "https://z/bucket/999")
    assert ("POST", "https://z/api/deposit/depositions/22683565/actions/newversion") in calls
    assert ("DELETE", "https://z/files/f1") in calls
    assert not any(u.endswith("/deposit/depositions") and m == "POST" for m, u in calls)


def test_new_version_fails_without_latest_draft(zen, monkeypatch):
    monkeypatch.setattr(zen, "make_request", lambda *a, **k: ({"links": {}}, 201))
    with pytest.raises(RuntimeError):
        zen.create_new_version("https://z/api", 1, "tok")
