# Errata: reverse-to-zero (2026-09-19)

These are corrections from the statistics and framing skeptics (runs wf_32d1ef7e-667 and wf_c11901b6-1ed). Script-generated JSON and code comments are left as committed, so reruns still reproduce them; this table gives the corrected reading. The hand-written Markdown files were corrected inline. Where anything disagrees, `REPORT.md` is authoritative.

| Where | As committed | Correct reading |
|---|---|---|
| `_prior_results_wf_32d1ef7e.json` ground E4; `E4-cosmic-web-tda/cosmic_web_tda.py:244,268`; `e4_cosmic_web_tda_report.json:58-59`; `e2_pta_and_tda_extension/IMPROVEMENT_PROPOSALS.md:80-82,99-102`; `IMPROVEMENT_PROPOSAL.md:16` | Quotes LeanMaster: "a TDA prediction would be ... This is frozen as Tier C", and "(H') ... to be frozen before comparison" as a LeanMaster E4 note | **No such text exists in LeanMaster.** A fixed-string grep over all tags finds nothing. LeanMaster STREAM8 says "Observables: none" (lines 112, 153, 374 at eb791e7). |
| ground E4 quote | Octad stabiliser in M24 = (Z2)^4 x\| A7 (order 40320) | It is (Z2)^4 x\| A8. (Z2)^4 x\| A7, order 40320, is TW's overarching group, a maximal subgroup of M23 (STREAM8:134-136). |
| ground P6.2 quote; `E3-M0/E3_RESULT.md:56` | kappa = 1, quoted without its qualifier | kappa = 1 holds only "with alpha' = s^2 (the Tier C dual-scale identification of Stream 3)". |
| ground P1 consequence | LeanMaster Stream 6 P1 and PRE_REGISTRATION P1 treated as the same | They differ: LeanMaster P1 is the extra-dimension radius; PRE_REGISTRATION P1 is the cosmological constant. |
| `exp:E3` parameter_effect | M0 "building on Tier A (LeanMaster kappa=1)" | M0 uses no kappa=1 input. Its Omega_Lambda = 0.68885 is an imported Planck18 value (tier L). |
| `e2_pta_result.json:107`; `e2_pta_sensitivity_scan.py:171` | "LeanMaster Streams 6-8 confirm mu_sym and c4_pta_product are unchanged" | LeanMaster never mentions either parameter. "Untested" is this project's statement. |
| `E3-M0/m0_model.py:295`; `E3_RESULT.md:40` | 9.00 called the "pre-registered rule" for M0 vs fitted LCDM | PRE_REGISTRATION.md:34-36 is only a threshold table. The M0 comparison was not pre-registered. |
| `E3_RESULT.md:40` | "~0.1 sigma" | 0.86 sigma: chi2.sf(0.746, 1) = 0.388. |
| `e1_desi_dr2.py:9`; exp:E1 verdict | "P1 NOT rejected; 3-sigma and 5-sigma bars missed" | P1's falsification rule names DR3/final BAO + CMB + SN, Euclid or Rubin. DR2 BAO+SN without CMB is not one of them, so no P1 verdict follows. The pipeline action was executed (PRE_REGISTRATION addendum A1). |
| `E2_RESULT.md:78` | slope 6.171 | 1.000 (`skeptic_statistics/check_e2_slope.json`). |
| `e5b_lognormal_and_poisson_null_report.json:1062` | "the actual M0 test" | Not an M0 test: the excess is confined to r < 20 Mpc/h, a range the gate and mocks do not validate. Verdict INCONCLUSIVE. The xi(r) gate has no power: xi = 0 passes it. |
| e5b "preregistered_rule" | Gate threshold called pre-registered | It was fixed in script source only. PRE_REGISTRATION.md has no TDA threshold (addendum A5). |
| `E5-cmb-tda/e5_cmb_tda_report.json` verdict | "M0's Gaussian prediction survives" | Not established. The null's pseudo-C_ell does not match the data (z = -6.6 / -11.6 / +11.0 in three ell bands), and the recalibration failed its gate. With about 2 effective tests, the corrected p for superlevel b1 is 0.022. |
| `e5_cmb_tda_report.json:2067` | "K3xT2 makes no ... prediction beyond the same Gaussian statistics" | K3xT2 makes no CMB-topology prediction. |
| E1 vs E3 | 1590 vs 1580 SNe | Different cuts. DR2 Delta chi2 is 4.547 and 4.735; the verdict is unchanged. |
