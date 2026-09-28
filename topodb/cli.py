"""Command-line query tool for TopoDB.

    python -m topodb.cli summary
    python -m topodb.cli search <query>
    python -m topodb.cli datasets [--domain DOMAIN]
    python -m topodb.cli runs [--dataset DATASET]
    python -m topodb.cli betti [--domain DOMAIN]
    python -m topodb.cli show <run_id>
    python -m topodb.cli export --json <path>
    python -m topodb.cli nearest <run_id> [--dim D] [--k K] [--metric bottleneck|wasserstein]
    python -m topodb.cli distance <run_a> <run_b> [--dim D] [--metric ...]

`--db PATH` overrides the database (default: $TOPODB_PATH, else the on-disk default).
"""
from __future__ import annotations

import argparse
import json
import sys

from .api import TopoDB


# ------------------------------------------------------------------ printing
def table(rows, columns=None, widths=None, out=None) -> None:
    """Print a list of dicts/sqlite3.Row as a plain aligned table."""
    out = out or sys.stdout
    rows = [dict(r) for r in rows]
    if not rows:
        print("(no rows)", file=out)
        return
    columns = columns or list(rows[0].keys())
    widths = widths or {}
    cells = [[_fmt(r.get(c), widths.get(c)) for c in columns] for r in rows]
    w = [max(len(str(c)), *(len(row[i]) for row in cells)) for i, c in enumerate(columns)]
    print("  ".join(str(c).ljust(w[i]) for i, c in enumerate(columns)), file=out)
    print("  ".join("-" * w[i] for i in range(len(columns))), file=out)
    for row in cells:
        print("  ".join(row[i].ljust(w[i]) for i in range(len(columns))), file=out)


def _fmt(v, maxw=None) -> str:
    if v is None:
        return "-"
    if isinstance(v, float):
        s = f"{v:.6g}"
    else:
        s = str(v)
    s = s.replace("\n", " ")
    if maxw and len(s) > maxw:
        s = s[: maxw - 1] + "…"
    return s


def kv(d: dict, indent: str = "  ", out=None) -> None:
    out = out or sys.stdout
    if not d:
        print(f"{indent}(none)", file=out)
        return
    w = max(len(str(k)) for k in d)
    for k, v in d.items():
        print(f"{indent}{str(k).ljust(w)} : {v}", file=out)


# ------------------------------------------------------------- subcommands
def cmd_summary(db: TopoDB, args) -> int:
    s = db.summary()
    print("TopoDB summary  (%s)" % db.path)
    for section in ("datasets_by_domain", "runs_by_method", "runs_by_tier",
                    "controls_passed", "findings_by_verdict"):
        print(f"\n{section}:")
        kv(s[section])
    print("\ntotals:")
    kv({"runs": db.con.execute("SELECT count(*) FROM run").fetchone()[0],
        "datasets": db.con.execute("SELECT count(*) FROM dataset").fetchone()[0],
        "betti_rows": db.con.execute("SELECT count(*) FROM betti").fetchone()[0],
        "bars": s["bars"],
        "statistics": db.con.execute("SELECT count(*) FROM statistic").fetchone()[0],
        "statistics_with_p_value": s["statistics_with_p"],
        "controls": db.con.execute("SELECT count(*) FROM control").fetchone()[0],
        "findings": db.con.execute("SELECT count(*) FROM finding").fetchone()[0]})
    print("\ntiers: A=machine-checked proof  B=exact arithmetic  L=literature  "
          "C=convention/definition  X=numerics")
    return 0


def cmd_search(db: TopoDB, args) -> int:
    rows = db.search(args.query, limit=args.limit)
    if not rows:
        print(f"no match for {args.query!r}")
        return 1
    table(rows, ["kind", "ref_id", "title", "hit"], {"title": 70, "hit": 60})
    return 0


def cmd_datasets(db: TopoDB, args) -> int:
    q = ("SELECT d.id, d.domain, d.provenance, d.n_objects, d.ambient_dim,"
         " (SELECT count(*) FROM run r WHERE r.dataset_id=d.id) AS runs, d.title"
         " FROM dataset d")
    a: tuple = ()
    if args.domain:
        q += " WHERE d.domain=?"
        a = (args.domain,)
    q += " ORDER BY d.domain, d.id"
    table(db.con.execute(q, a).fetchall(),
          ["id", "domain", "provenance", "n_objects", "ambient_dim", "runs", "title"],
          {"id": 48, "title": 60})
    return 0


def cmd_runs(db: TopoDB, args) -> int:
    q = ("SELECT r.id, r.dataset_id, r.method, r.coeff_field AS field, r.max_dim, r.tier,"
         " (SELECT count(*) FROM bar b WHERE b.run_id=r.id) AS bars,"
         " (SELECT count(*) FROM statistic s WHERE s.run_id=r.id) AS stats,"
         " (SELECT count(*) FROM control c WHERE c.run_id=r.id) AS ctrls,"
         " (SELECT group_concat(b.dim||':'||b.value,' ') FROM betti b WHERE b.run_id=r.id) AS betti"
         " FROM run r")
    a: tuple = ()
    where = []
    if args.dataset:
        where.append("r.dataset_id LIKE ?")
        a += (f"%{args.dataset}%",)
    if args.tier:
        where.append("r.tier=?")
        a += (args.tier,)
    if where:
        q += " WHERE " + " AND ".join(where)
    q += " ORDER BY r.id LIMIT ?"
    table(db.con.execute(q, (*a, args.limit)).fetchall(),
          ["id", "dataset_id", "method", "field", "max_dim", "tier", "betti", "bars", "stats", "ctrls"],
          {"dataset_id": 46, "betti": 34})
    return 0


def cmd_betti(db: TopoDB, args) -> int:
    q = ("SELECT d.domain, r.dataset_id AS dataset, r.id AS run, r.method, r.tier,"
         " (SELECT group_concat(b.dim||':'||b.value,' ') FROM betti b WHERE b.run_id=r.id) AS observed,"
         " (SELECT group_concat(b.dim||':'||coalesce(b.expected,'?'),' ') FROM betti b"
         "   WHERE b.run_id=r.id AND b.expected IS NOT NULL) AS expected,"
         " (SELECT CASE WHEN min(b.matches) IS NULL THEN '-'"
         "              WHEN min(b.matches)=1 THEN 'match' ELSE 'MISMATCH' END"
         "   FROM betti b WHERE b.run_id=r.id) AS agree"
         " FROM run r JOIN dataset d ON d.id=r.dataset_id"
         " WHERE EXISTS (SELECT 1 FROM betti b WHERE b.run_id=r.id)")
    a: tuple = ()
    if args.domain:
        q += " AND d.domain=?"
        a = (args.domain,)
    q += " ORDER BY d.domain, r.dataset_id, r.id LIMIT ?"
    rows = db.con.execute(q, (*a, args.limit)).fetchall()
    table(rows, ["domain", "dataset", "run", "method", "tier", "observed", "expected", "agree"],
          {"dataset": 44, "observed": 30, "expected": 30})
    print("\n'expected' is present only where the source pre-stated the expectation; "
          "'-' under agree means no expectation was recorded.")
    return 0


def cmd_show(db: TopoDB, args) -> int:
    r = db.con.execute("SELECT * FROM run WHERE id=?", (args.run_id,)).fetchone()
    if r is None:
        print(f"no run {args.run_id}", file=sys.stderr)
        return 1
    d = db.con.execute("SELECT * FROM dataset WHERE id=?", (r["dataset_id"],)).fetchone()
    print(f"run {r['id']}  [tier {r['tier']}]")
    print("\ndataset:")
    kv({k: d[k] for k in d.keys() if d[k] is not None})
    print("\nrun:")
    kv({k: r[k] for k in r.keys() if k != "params_json" and r[k] is not None})
    print("\nparams:")
    try:
        kv(json.loads(r["params_json"]))
    except Exception:
        print("  " + str(r["params_json"]))

    print("\nbetti:")
    table(db.con.execute(
        "SELECT dim, value AS observed, expected,"
        " CASE matches WHEN 1 THEN 'match' WHEN 0 THEN 'MISMATCH' ELSE '-' END AS agree"
        " FROM betti WHERE run_id=? ORDER BY dim", (args.run_id,)).fetchall())

    print("\ncontrols:")
    table(db.con.execute(
        "SELECT kind, CASE passed WHEN 1 THEN 'passed' WHEN 0 THEN 'FAILED' ELSE 'not_run' END"
        " AS outcome, description, detail FROM control WHERE run_id=? ORDER BY kind", (args.run_id,)
    ).fetchall(), ["kind", "outcome", "description", "detail"], {"description": 60, "detail": 60})

    print("\nstatistics:")
    table(db.con.execute(
        "SELECT name, value, p_value, p_method, null_model, n_null, multiplicity"
        " FROM statistic WHERE run_id=? ORDER BY name", (args.run_id,)).fetchall(),
        ["name", "value", "p_value", "p_method", "null_model", "n_null", "multiplicity"],
        {"name": 40, "null_model": 46, "multiplicity": 26})
    print("  (a statistic with no p_value is one whose source reported no null model; "
          "no p-value was invented for it)")

    print(f"\ntop bars (longest {args.bars} per dimension; the DB stores at most the top 50):")
    dims = [row[0] for row in db.con.execute(
        "SELECT DISTINCT dim FROM bar WHERE run_id=? ORDER BY dim", (args.run_id,)).fetchall()]
    if not dims:
        print("  (no bars stored for this run — exact chain-level runs have no filtration)")
    for dim in dims:
        n = db.con.execute("SELECT count(*) FROM bar WHERE run_id=? AND dim=?",
                           (args.run_id, dim)).fetchone()[0]
        print(f"  H{dim}  ({n} bars stored)")
        table(db.con.execute(
            "SELECT rank, birth, death, CASE WHEN death IS NULL THEN NULL ELSE death-birth END"
            " AS persistence FROM bar WHERE run_id=? AND dim=? ORDER BY rank LIMIT ?",
            (args.run_id, dim, args.bars)).fetchall())

    print("\nfindings:")
    table(db.con.execute(
        "SELECT id, verdict, tier, claim, caveat, reference FROM finding WHERE run_id=? ORDER BY id",
        (args.run_id,)).fetchall(), ["id", "verdict", "tier", "claim", "caveat", "reference"],
        {"claim": 68, "caveat": 60, "reference": 44})
    return 0


def cmd_export(db: TopoDB, args) -> int:
    from . import report as _report
    data = _report.write(args.json, db)
    print(f"wrote {args.json}")
    print(json.dumps(data["totals"], indent=2))
    return 0


def cmd_nearest(db: TopoDB, args) -> int:
    from . import compare as _compare
    res = _compare.nearest(args.run_id, k=args.k, dim=args.dim, db=db, metric=args.metric,
                           same_method=not args.any_method, infinite=args.infinite)
    ref = res["reference"]
    print(f"run {args.run_id}  ({ref['dataset']}, {ref['method']}, tier {ref['tier']})  "
          f"H{args.dim}: {ref['n_points_used']} finite points of {ref['n_bars_stored']} stored, "
          f"{ref['n_infinite_dropped']} infinite dropped")
    print(f"metric: {args.metric}")
    table(res["neighbours"], ["run_id", "distance", "dataset", "method", "tier",
                              "n_points_used", "n_bars_stored"], {"dataset": 46})
    print(f"\ncompared against {res['n_compared']} run(s); "
          f"{res['excluded_no_bars']} run(s) excluded for having no stored bars in H{args.dim} "
          f"(exact chain-level runs have no filtration); "
          f"{res['excluded_other_method']} excluded as a different filtration method.")
    print("APPROXIMATE: " + res["caveat"])
    return 0


def cmd_distance(db: TopoDB, args) -> int:
    from . import compare as _compare
    fn = _compare.bottleneck if args.metric == "bottleneck" else _compare.wasserstein
    d, info = fn(args.run_a, args.run_b, args.dim, db=db, infinite=args.infinite)
    print(f"{info['metric']}  H{args.dim}  run {args.run_a} vs run {args.run_b}  =  {d:.6g}")
    kv({"run_a": info["run_a"], "run_b": info["run_b"],
        "methods_comparable": info["methods_comparable"]})
    if not info["methods_comparable"]:
        print("\nWARNING: the two runs use different filtration methods; their filtration values "
              "are in different units and this number is not interpretable as a similarity.")
    print("APPROXIMATE: " + info["caveat"])
    return 0


# ----------------------------------------------------------------- argparse
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="topodb", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--db", default=None, help="path to topodb.sqlite (default: $TOPODB_PATH)")
    sub = p.add_subparsers(dest="cmd", required=True)

    sub.add_parser("summary", help="counts by domain, method, tier, verdict").set_defaults(fn=cmd_summary)

    s = sub.add_parser("search", help="full-text search over dataset titles and finding claims")
    s.add_argument("query")
    s.add_argument("--limit", type=int, default=20)
    s.set_defaults(fn=cmd_search)

    s = sub.add_parser("datasets", help="list datasets")
    s.add_argument("--domain", choices=["astro", "quantum_fluid", "biology", "mathematics", "synthetic"])
    s.set_defaults(fn=cmd_datasets)

    s = sub.add_parser("runs", help="list runs")
    s.add_argument("--dataset", help="substring of the dataset id")
    s.add_argument("--tier", choices=list("ABLCX"))
    s.add_argument("--limit", type=int, default=500)
    s.set_defaults(fn=cmd_runs)

    s = sub.add_parser("betti", help="Betti numbers, with the pre-stated expectation where there was one")
    s.add_argument("--domain", choices=["astro", "quantum_fluid", "biology", "mathematics", "synthetic"])
    s.add_argument("--limit", type=int, default=500)
    s.set_defaults(fn=cmd_betti)

    s = sub.add_parser("show", help="full record of one run")
    s.add_argument("run_id", type=int)
    s.add_argument("--bars", type=int, default=8, help="top bars per dimension to print")
    s.set_defaults(fn=cmd_show)

    s = sub.add_parser("export", help="write the machine-readable DB report")
    s.add_argument("--json", required=True, metavar="PATH")
    s.set_defaults(fn=cmd_export)

    s = sub.add_parser("nearest", help="k topologically nearest runs (APPROXIMATE, see the output)")
    s.add_argument("run_id", type=int)
    s.add_argument("--k", type=int, default=5)
    s.add_argument("--dim", type=int, default=1)
    s.add_argument("--metric", choices=["bottleneck", "wasserstein"], default="bottleneck")
    s.add_argument("--any-method", action="store_true",
                   help="also rank runs from a different filtration (units are incommensurable)")
    s.add_argument("--infinite", choices=["drop", "cap", "keep"], default="drop")
    s.set_defaults(fn=cmd_nearest)

    s = sub.add_parser("distance", help="bottleneck/Wasserstein between two runs (APPROXIMATE)")
    s.add_argument("run_a", type=int)
    s.add_argument("run_b", type=int)
    s.add_argument("--dim", type=int, default=1)
    s.add_argument("--metric", choices=["bottleneck", "wasserstein"], default="bottleneck")
    s.add_argument("--infinite", choices=["drop", "cap", "keep"], default="drop")
    s.set_defaults(fn=cmd_distance)
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    db = TopoDB(args.db) if args.db else TopoDB()
    return args.fn(db, args)


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
