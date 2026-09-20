#!/usr/bin/env python3
"""Four corrections to rows already in TopoDB, after a review of the record.

Each is a row edit, not a re-computation; the numbers are unchanged.

F1  Five point-cloud `null_calibration passed=1` rows are the same vacuous
    control that was removed from the sky maps: 20 DISJOINT draws from one
    random catalogue are exchangeable BY CONSTRUCTION, so the check cannot
    fail.  "A check that cannot fail is not a check" (LeanFlow CLAUDE.md
    rule 2).  Set them to passed=NULL with the honest description.

F2  Run 146 (LRG, absolute cap) carries a `negative passed=0` control saying
    its H2 statistics are saturated and must not be read, but its FINDING says
    only "recovered, 4 of 6 outside the null range" with the generic DESI
    caveat.  A reader pulling `finding` never sees the control.  Append the
    warning, and name which two statistics are the missing 2 of 6.

F3  COBE-DMR's rejection is an ARTEFACT, measured:
      - 1 to 4 DMR quad-cube pixels land in each nside-16 HEALPix pixel
        (702 / 1758 / 522 / 90 pixels get 1 / 2 / 3 / 4), a coverage ratio of 4;
      - SERROR varies by 1.92x across the map and N_OBS by 3.68x.
    An isotropic Gaussian null contains neither the repixelisation unevenness
    nor DMR's anisotropic scan coverage, so it would reject on those alone,
    before any sky signal.  f_sky = 1.0 only says every pixel got at least one
    sample; it says nothing about how many.  Verdict recovered -> artefact.

F4  Planck SMICA and WMAP ILC are the SAME SKY and disagree: SMICA has 0 of 5
    statistics outside the null range (max separation 2.1 sigma) while the ILC
    rejects at the coarse-curve floor on both b0 and b1.  Two maps of one sky
    cannot both be reporting the sky, so the statistic is reading each map's
    PROCESSING.  Neither finding referenced the other.  Cross-reference both,
    and downgrade the ILC verdict accordingly.

Tier X throughout.  Nothing here is "proved".
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import topodb  # noqa: E402

VACUOUS = ("DIAGNOSTIC, passed=None by design: a held-out null realisation ranked against the "
           "remaining n-1 tests exchangeability AMONG THE NULL DRAWS, which is TRUE BY CONSTRUCTION "
           "for disjoint subsamples of one random catalogue. It has no power to detect a null "
           "mismatched to the data (every draw shares the same mismatch) and, being one Uniform(0,1) "
           "draw, it fails ~5% of the time on a perfectly good null. It never gates a verdict.")

LRG_WARN = (" SATURATED, DO NOT READ THE H2 STATISTICS OF THIS RUN. The pre-declared absolute alpha "
            "cap of 30 Mpc/h is only 2.1x this sample's own H0 death median (14.12 Mpc/h), so the "
            "complex is truncated before voids of the sample's characteristic size can be born. The "
            "data AND all 20 nulls contain exactly zero H2 bars over 5 Mpc/h - a statistic with no "
            "variance in the null, which is saturation, not a null result. That is precisely why this "
            "run shows 4 of 6 outside the null range and not 5 of 6: the two that are not are "
            "n_h2_bars_over_5mpc and max_persistence_h2, the saturated pair. See AMENDMENT A3 and the "
            "scale-free re-run (run 155), where lifting the cap moves the longest H2 bar from 3.81 to "
            "34.33 Mpc/h.")

COBE_ART = (" THIS REJECTION IS AN ARTEFACT OF THE REPIXELISATION AND OF DMR'S OWN NOISE, not a "
            "measurement of the sky. Measured: 1 to 4 DMR quad-cube pixels land in each nside-16 "
            "HEALPix pixel (702/1758/522/90 HEALPix pixels receive 1/2/3/4), a coverage ratio of 4; "
            "SERROR varies by 1.92x across the map and N_OBS by 3.68x. COBE-DMR at 53 GHz is "
            "noise-dominated at its 7 deg beam. An isotropic Gaussian null contains NEITHER the "
            "uneven repixelisation NOR the anisotropic scan coverage, so it must reject on those "
            "alone, before any sky signal is considered. f_sky = 1.0000 only certifies that every "
            "HEALPix pixel received at least one DMR sample; it says nothing about how many. The null "
            "is inapplicable to this map and no departure of the SKY from Gaussianity is claimed.")

SAMESKY = (" SAME-SKY DISCREPANCY, and it is the most informative thing in this pair: Planck SMICA "
           "(run 152) has 0 of 5 statistics outside the null range with a maximum separation of 2.1 "
           "sigma, while WMAP 9-yr ILC (run 201) rejects at the coarse-curve rank floor 1/101 on both "
           "b0 and b1. These are two maps of ONE sky, so they cannot both be reporting the sky: the "
           "statistic is reading each map's PROCESSING, not the CMB. Candidate causes, neither of "
           "which an isotropic Gaussian null contains: the WMAP ILC's spatially varying 12-region "
           "weighting, and WMAP's anisotropic scan noise (Planck SMICA is a different component "
           "separation on a different instrument with a different noise geometry). No cosmological "
           "departure is claimed from either map.")

SPACING = (" One clause on h0_death_iqr_over_median, because its very large separation will otherwise "
           "be read as revising the prior campaign's lesson: that lesson ('in 3-D use the COUNT, not "
           "the spacing spread') was measured against a CLUSTERING-MATCHED null and an injected "
           "signal. What is separated here is a clustered field from an UNCLUSTERED one, which is a "
           "far easier discrimination. This number is not evidence that the spacing spread works as a "
           "defect detector in 3-D, and it does not revise that lesson.")

PEAKMEM = ("NOTE on peak_mb for multi-dataset scripts: ru_maxrss is a PROCESS high-water mark, so the "
           "2nd and later runs of one script inherit the peak of the first and their peak_mb is an "
           "upper bound, not an attribution.")


def main():
    with topodb() as db:
        c = db.con

        # ---- F1
        rows = c.execute("SELECT r.id FROM control x JOIN run r ON r.id=x.run_id "
                         "WHERE x.kind='null_calibration' AND x.passed=1 AND r.method='alpha'").fetchall()
        ids = [r[0] for r in rows]
        n1 = c.execute("UPDATE control SET passed=NULL, description=? WHERE kind='null_calibration' "
                       "AND passed=1 AND run_id IN (SELECT id FROM run WHERE method='alpha')",
                       (VACUOUS,)).rowcount
        c.commit()
        print(f"F1: {n1} vacuous null_calibration rows set to passed=NULL (runs {ids})")

        # ---- F2
        n2 = c.execute("UPDATE finding SET caveat = caveat || ? WHERE run_id=146", (LRG_WARN,)).rowcount
        c.commit()
        print(f"F2: {n2} LRG finding(s) now carry the saturation warning")

        # ---- F3
        n3a = c.execute("UPDATE finding SET verdict='artefact', caveat = caveat || ? "
                        "WHERE dataset_id='astro/cobe_dmr'", (COBE_ART,)).rowcount
        c.execute("INSERT OR REPLACE INTO control(run_id,kind,description,passed,detail) "
                  "SELECT id,'negative','the Gaussian null is INAPPLICABLE to this map: the "
                  "repixelisation is uneven (1-4 DMR pixels per HEALPix pixel, ratio 4) and DMR is "
                  "noise-dominated with anisotropic coverage (SERROR varies 1.92x, N_OBS 3.68x), "
                  "neither of which an isotropic Gaussian null contains',0,"
                  "'702/1758/522/90 nside-16 pixels receive 1/2/3/4 DMR samples; SERROR 0.0572-0.1097; "
                  "N_OBS 44469-163713' FROM run WHERE dataset_id='astro/cobe_dmr'")
        c.commit()
        print(f"F3: {n3a} COBE finding(s) -> verdict 'artefact', with the measured diagnostic")

        # ---- F4
        n4 = c.execute("UPDATE finding SET caveat = caveat || ? WHERE dataset_id IN "
                       "('astro/planck_smica','astro/wmap_ilc')", (SAMESKY,)).rowcount
        c.execute("UPDATE finding SET verdict='inconclusive' WHERE dataset_id='astro/wmap_ilc'")
        c.commit()
        print(f"F4: {n4} CMB findings cross-referenced; wmap_ilc verdict 'recovered' -> 'inconclusive'")

        # ---- minor
        n5 = c.execute("UPDATE finding SET caveat = caveat || ? WHERE dataset_id LIKE 'astro/desi%' "
                       "OR dataset_id='astro/sdss_spec_bulk'", (SPACING,)).rowcount
        c.execute("UPDATE run SET preprocessing = preprocessing || ? WHERE dataset_id LIKE 'astro/%' "
                  "AND preprocessing IS NOT NULL", (" || " + PEAKMEM,))
        c.commit()
        print(f"minor: {n5} point-cloud findings carry the spacing-statistic clause; peak_mb note on all astro runs")

        print("\n--- astro verdicts after the fixes ---")
        for r in c.execute("SELECT dataset_id, verdict FROM finding WHERE dataset_id LIKE 'astro/%' "
                           "ORDER BY dataset_id").fetchall():
            print(f"  [{r[1]:12s}] {r[0]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
