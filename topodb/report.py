"""Machine-written JSON summary of TopoDB.

Nothing in here is an interpretation: every field is a count or a verbatim
row from the database. The one editorial choice is the grouping, and it is
stated in the output under `_what_this_is`.
"""
from __future__ import annotations

import datetime as _dt
import json

from .api import TopoDB


def _rows(db: TopoDB, sql: str, args: tuple = ()) -> list[dict]:
    return [dict(r) for r in db.con.execute(sql, args).fetchall()]


def _counts(db: TopoDB, sql: str, args: tuple = ()) -> dict:
    return {str(k): int(v) for k, v in db.con.execute(sql, args).fetchall()}


def build(db: TopoDB | None = None) -> dict:
    db = db or TopoDB()
    out: dict = {
        "_what_this_is": (
            "A machine-written census of the TopoDB tables. Counts only; no claim is made "
            "here about any physical or mathematical result. Tiers are as recorded at "
            "ingestion: A/B/L/C/X. A run counted under 'controls_passed' passed the control "
            "its own source declared, nothing wider."
        ),
        "generated": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
        "db_path": db.path,
        "totals": {
            "datasets": db.con.execute("SELECT count(*) FROM dataset").fetchone()[0],
            "runs": db.con.execute("SELECT count(*) FROM run").fetchone()[0],
            "betti_rows": db.con.execute("SELECT count(*) FROM betti").fetchone()[0],
            "bars": db.con.execute("SELECT count(*) FROM bar").fetchone()[0],
            "statistics": db.con.execute("SELECT count(*) FROM statistic").fetchone()[0],
            "statistics_with_p_value": db.con.execute(
                "SELECT count(*) FROM statistic WHERE p_value IS NOT NULL").fetchone()[0],
            "statistics_without_p_value": db.con.execute(
                "SELECT count(*) FROM statistic WHERE p_value IS NULL").fetchone()[0],
            "controls": db.con.execute("SELECT count(*) FROM control").fetchone()[0],
            "findings": db.con.execute("SELECT count(*) FROM finding").fetchone()[0],
        },
    }
    out["datasets_by_domain"] = _counts(db, "SELECT domain, count(*) FROM dataset GROUP BY domain")
    out["datasets_by_provenance"] = _counts(
        db, "SELECT provenance, count(*) FROM dataset GROUP BY provenance")
    out["runs_by_domain"] = _counts(
        db, "SELECT d.domain, count(*) FROM run r JOIN dataset d ON d.id=r.dataset_id GROUP BY d.domain")
    out["runs_by_method"] = _counts(db, "SELECT method, count(*) FROM run GROUP BY method")
    out["runs_by_tier"] = _counts(db, "SELECT tier, count(*) FROM run GROUP BY tier")
    out["runs_by_domain_and_tier"] = {
        f"{r['domain']}/{r['tier']}": r["n"] for r in _rows(
            db, "SELECT d.domain, r.tier, count(*) AS n FROM run r JOIN dataset d ON d.id=r.dataset_id"
                " GROUP BY d.domain, r.tier ORDER BY d.domain, r.tier")}
    out["runs_by_coeff_field"] = _counts(
        db, "SELECT coeff_field, count(*) FROM run GROUP BY coeff_field")
    out["findings_by_verdict"] = _counts(db, "SELECT verdict, count(*) FROM finding GROUP BY verdict")
    out["findings_by_tier"] = _counts(db, "SELECT tier, count(*) FROM finding GROUP BY tier")
    out["controls_by_kind_and_outcome"] = {
        f"{r['kind']}/{r['outcome']}": r["n"] for r in _rows(
            db, "SELECT kind, CASE passed WHEN 1 THEN 'passed' WHEN 0 THEN 'failed' ELSE 'not_run' END"
                " AS outcome, count(*) AS n FROM control GROUP BY kind, outcome ORDER BY kind, outcome")}

    # --- which datasets have a KNOWN-ANSWER control that passed -------------
    out["known_answer_controls"] = {
        "_meaning": ("a dataset is listed under `all_passed` when it carries at least one control of "
                     "kind 'known_answer' and none of its known-answer controls failed"),
        "all_passed": _rows(db, """
            SELECT d.id AS dataset, d.domain, count(*) AS n_known_answer_controls
            FROM control c JOIN run r ON r.id=c.run_id JOIN dataset d ON d.id=r.dataset_id
            WHERE c.kind='known_answer'
            GROUP BY d.id
            HAVING sum(CASE WHEN c.passed=0 THEN 1 ELSE 0 END)=0
               AND sum(CASE WHEN c.passed=1 THEN 1 ELSE 0 END)>0
            ORDER BY d.domain, d.id"""),
        "some_failed": _rows(db, """
            SELECT d.id AS dataset, d.domain,
                   sum(CASE WHEN c.passed=0 THEN 1 ELSE 0 END) AS n_failed,
                   sum(CASE WHEN c.passed=1 THEN 1 ELSE 0 END) AS n_passed
            FROM control c JOIN run r ON r.id=c.run_id JOIN dataset d ON d.id=r.dataset_id
            WHERE c.kind='known_answer' GROUP BY d.id
            HAVING n_failed > 0 ORDER BY d.domain, d.id"""),
        "datasets_with_no_known_answer_control": _rows(db, """
            SELECT d.id AS dataset, d.domain FROM dataset d
            WHERE d.id NOT IN (SELECT r.dataset_id FROM control c JOIN run r ON r.id=c.run_id
                               WHERE c.kind='known_answer')
            ORDER BY d.domain, d.id"""),
    }

    # --- betti agreement where an expectation was pre-stated ----------------
    out["betti_vs_prestated_expectation"] = {
        "_meaning": ("only rows whose source pre-stated an expectation carry `expected`; rows with "
                     "expected IS NULL were never gated against an expectation and are not counted here"),
        "matched": db.con.execute("SELECT count(*) FROM betti WHERE matches=1").fetchone()[0],
        "mismatched": db.con.execute("SELECT count(*) FROM betti WHERE matches=0").fetchone()[0],
        "no_expectation_recorded": db.con.execute(
            "SELECT count(*) FROM betti WHERE expected IS NULL").fetchone()[0],
        "mismatches": _rows(db, """
            SELECT b.run_id, d.id AS dataset, r.method, b.dim, b.value AS observed, b.expected
            FROM betti b JOIN run r ON r.id=b.run_id JOIN dataset d ON d.id=r.dataset_id
            WHERE b.matches=0 ORDER BY d.id, b.run_id, b.dim"""),
    }

    # --- the honest bit: what came out null or inconclusive -----------------
    out["findings_null_or_inconclusive"] = _rows(db, """
        SELECT f.id, f.run_id, coalesce(f.dataset_id, r.dataset_id) AS dataset, f.verdict, f.tier,
               f.claim, f.caveat, f.reference
        FROM finding f LEFT JOIN run r ON r.id=f.run_id
        WHERE f.verdict IN ('null','inconclusive') ORDER BY f.verdict, f.id""")
    out["findings_failed_or_artefact"] = _rows(db, """
        SELECT f.id, f.run_id, coalesce(f.dataset_id, r.dataset_id) AS dataset, f.verdict, f.tier,
               f.claim, f.caveat, f.reference
        FROM finding f LEFT JOIN run r ON r.id=f.run_id
        WHERE f.verdict IN ('failed','artefact') ORDER BY f.verdict, f.id""")
    out["findings_recovered"] = _rows(db, """
        SELECT f.id, f.run_id, coalesce(f.dataset_id, r.dataset_id) AS dataset, f.tier, f.claim, f.caveat
        FROM finding f LEFT JOIN run r ON r.id=f.run_id WHERE f.verdict='recovered' ORDER BY f.id""")

    out["statistics_carrying_a_p_value"] = _rows(db, """
        SELECT s.run_id, d.id AS dataset, s.name, s.value, s.p_value, s.p_method, s.null_model,
               s.n_null, s.multiplicity
        FROM statistic s JOIN run r ON r.id=s.run_id JOIN dataset d ON d.id=r.dataset_id
        WHERE s.p_value IS NOT NULL ORDER BY s.p_value""")
    out["runs_without_any_bars"] = {
        "_meaning": "exact chain-level / Mayer-Vietoris runs carry Betti numbers but no filtration",
        "count": db.con.execute(
            "SELECT count(*) FROM run WHERE id NOT IN (SELECT DISTINCT run_id FROM bar)").fetchone()[0],
        "by_method": _counts(db, "SELECT method, count(*) FROM run"
                                 " WHERE id NOT IN (SELECT DISTINCT run_id FROM bar) GROUP BY method"),
    }
    return out


def write(path: str, db: TopoDB | None = None) -> dict:
    data = build(db)
    with open(path, "w") as fh:
        json.dump(data, fh, indent=2, sort_keys=False)
        fh.write("\n")
    return data


if __name__ == "__main__":  # pragma: no cover
    print(json.dumps(build(), indent=2))
