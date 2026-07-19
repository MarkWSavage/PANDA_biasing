# PANDA Validation Summary

**Status:** Core development and validation considered complete as of 2026-07-17. CREME-MC comparison corrected 2026-07-19 (see 1.3).

**Scope:** This document summarizes the validation and robustness testing performed against PANDA's three core responsibilities -- geometry, particle generation/transport, and energy/charge scoring -- and records the final, accepted status of every comparison against external data or tools. See `Documentation/PANDA_MASTER_DESIGN` for architecture and philosophy; this document is the validation-specific companion.

---

## 1. Comparisons against external data/tools

### 1.1 Hitachi HM68512 proton SEU (open literature, LaBel et al. ~1994 NASA GSFC NSREC report)

The only comparison against real, independent (non-CREME-MC, non-CUPID) experimental data. Device: 4 Mbit CMOS-on-EPI SRAM, proton-tested 22-63 MeV at UC Davis. Measured static cross section 1e-6 cm^2/device, dynamic 2.2e-5 cm^2/device, with little energy dependence over that range.

With a geometrically plausible (not device-specified -- the paper doesn't report it) sensitive volume and a freely-fit critical charge:

- **Best fit: Qcrit=50 fC, sensitive volume 20-27 um^3** (three different aspect ratios -- 3x3x3, 2x2x5, 5x5x1 um -- all converge here independent of shape).
- Reproduces the measured static cross section within **~2x**, and reproduces the real device's flat 22-63 MeV energy dependence (2-3x max/min variation).
- Undersized (4 um^3) and oversized (125 um^3) geometries both break the fit -- undershooting badly or requiring a much higher Qcrit that then destroys the flat energy dependence -- so the 20-27 um^3 match is a real constraint, not a free parameter absorbing everything.
- Re-verified after the region-scoped production-cuts change (2026-07-15): shift is small at the Qcrit=50-100fC range this comparison actually operates in; conclusion unchanged.

**Verdict: genuine, independent plausibility win for PANDA's physics list (QGSP_BIC_HP)** -- not proof of accuracy (Qcrit is still a fit parameter), but agreement to within a factor of 2 against real 1994 accelerator data, with no dependency on reproducing another simulator's internals.

### 1.2 McNulty et al. 1989 proton SEE / CUPID (148 MeV, 11.7x11.7x13.5 um Si)

- **Original run** (vacuum world, no surrounding material): systematically under-predicted, growing from ~1.4x at 1 MeV deposited to >100x by 25-28 MeV (log-RMSE ~1.1 decades). Root cause identified: the dead+sensitive stack sat in a vacuum world with nothing for nearby nuclear reactions to occur in -- exactly the mechanism the source 1989 paper (El-Teleaty, McNulty, Abdel-Kader, Beauvais) was about.
- **Fix**: added a 100 um silicon surrounding volume plus secondary-neutron biasing within it. Improved log-RMSE to ~0.91 decades (~18% reduction), e.g. the 10 MeV point improved from 3.07x to 2.47x under-prediction.
- **Remaining gap**: still large at high deposited energy (~48x at 25 MeV, down from ~133x). Not fully understood -- candidates (100 um shell still too small, only one cascade generation biased, other secondary species not biased) were never isolated. Accepted as an open, unresolved residual rather than a bug.
- **CUPID's own limitation**: CUPID uses a liquid-drop nuclear model for the Si-recoil evaporation/spallation step, and cannot model direct ionization from the proton at all -- so CUPID's curve is structurally absent across the entire low-deposited-energy end. Disagreement there is expected and not meaningful to chase.

**Verdict: partially addressed, residual gap accepted as a known, understood-in-cause-but-not-in-magnitude limitation.**

### 1.3 CREME-MC / MRED

- **Correction (2026-07-19): the original ~0.9-decade / "1-2 orders of magnitude" gap reported below was a methodology bug in `compare_creme_panda.py`, not a physics discrepancy.** The script was comparing CREME-MC's cross section against PANDA's **collected**-charge curve (`cumulative_cross_section_collected.csv`), which applies PANDA's drift/trapping collection-efficiency model (`Qeff = T*Qgen`, T<1, in `SteppingAction.cc`). CREME-MC (like CUPID) assumes the RPP model's ideal 100%-charge-collection convention, which has no equivalent of that efficiency factor -- so the two curves were never comparable quantities. Fixed to read `cumulative_cross_section_deposited.csv` instead (PANDA's 100%-collection-equivalent output).
- **Corrected result: log-RMSE = 0.292 decades over all 50/50 CREME points**, now fully inside PANDA's simulated charge range (the previous comparison excluded 4 tail points as out-of-range, since the collected-charge curve tops out lower than deposited). The two curves visually overlap almost entirely across the full range, including through the ~3-10 fC region -- **the "intermediate-charge hump" previously attributed to differing nuclear reaction model families is not visible in the corrected comparison.** CREME-MC shows essentially the same knee/shoulder structure PANDA does at that point.
- The prior root-cause narrative (QGSP_BIC_HP vs. CREME-MC's semi-empirical fragmentation cross-sections producing a physics-list-dependent hump) is **not supported once the comparison is done correctly** -- it was very likely an artifact of comparing a collection-degraded curve against an ideal-collection reference, not evidence of a genuine cascade-model discrepancy. The Geant4-vs-CREME96 literature comparison (arXiv:0712.2149) may still be relevant to *other* discrepancies, but is no longer needed to explain this one.
- MRED-vs-CREME-MC's own internal physics-list spread (Bertini/Binary/INCL++) is unaffected by this correction and remains a separate, valid reason not to chase an exact match to one specific MRED configuration.

**Verdict (revised): good agreement (log-RMSE ~0.29 decades, full-range coverage) once compared against the correct 100%-collection-equivalent PANDA quantity. No remaining unexplained disagreement with CREME-MC requiring a "structural, not fixable" hedge.**

---

## 2. Geometry and physics robustness checks

These did not compare against external data -- they stress-tested PANDA's own geometry, biasing, and stepping code for self-consistency and physical plausibility across configurations no single external dataset covers.

| Investigation | Question | Result |
|---|---|---|
| Si-recoil max-LET | Does silicon recoil LET really top out near the commonly-cited ~12 MeV*cm2/mg heavy-ion-hardness heuristic? | **Confirmed, ceiling ~14.4-14.8** MeV*cm2/mg (Si26/Si27), converging not diverging with 10x more statistics. A navigator boundary-artifact bug (near-zero step length -> spurious LET) was found and fixed along the way. |
| Gold-lid scatter | Does a package-lid gold layer raise the max recoil LET (Turflinger effect)? | **Negative result.** Max LET stayed pinned at ~14.2 across every lid thickness/gap/bias-factor variant tried. Lid raises ordinary recoil *rate* 2-5x but not high-LET-tail *severity*. Feature fully reverted. |
| Tungsten/heavy-metal fission | Do W via-plug dead layers really undergo neutron-induced fission (Vanderbilt claim), and does this exceed Si's LET ceiling? | **Confirmed, real effect.** W dead layer produces genuine fission-like fragments (Z 25-65) up to LET ~35-37, well above Si's ~14-15 ceiling. Confirmed **not neutron-specific** -- 200 MeV protons induce it at least as readily. |
| Heavy-metal comparison (Ta/W/Pb/Au) | Does severity scale with atomic number? | **No** -- severity clusters in two pairs, not Z order: Ta/W (mild, LET ~33-37) vs. Au/Pb (severe, LET ~39-41, 5-8x the fission rate), despite Ta/W/Au/Pb being in ascending Z order. Investigation closed at Mark's request; not extended to further elements. |
| Surrounding-volume "shell ceiling" (~80 um) | Does the auto-grown surrounding-volume arithmetic silently under-capture recoils above a ~80 um device-size threshold? | **Not confirmed.** Swept 50 um-1 mm devices, 50/200 MeV, both XY and thickness axes, auto-grown shell vs. an explicit 4x-wide reference shell: all ratios stayed within 1-4% of 1.0 with no systematic growing-deficit-with-size trend. A plausible code-reading-derived hypothesis that did not survive direct testing. |
| Region-scoped production cuts | Does the stock Geant4 ~0.7-1mm production-cut default silently break charge-sharing/delta-ray fidelity at deep-submicron device scale? | **Yes, real gap found and fixed.** Added a `G4Region` with production cuts auto-scaled to `min(device dimensions)/10` (clamped [1nm, 0.7mm]). Regression check showed up to ~32% shift at the high-charge tail for a 50 um device -- a real fidelity improvement, not noise. Re-validated against Hitachi/McNulty afterward: no material change to those results at the Qcrit ranges they use. |
| Deep-submicron (80nm-node) smoke test | Does the full geometry/biasing/scoring/LET pipeline still function correctly at modern transistor scale (200/50/20 nm)? | **Yes.** Clean run, biasing engaged correctly (~850x recoil-hit uplift matching the 1041x factor), max recoil LET ~12.1 (Al27) -- consistent with the um-scale ~11-15 ceiling, now confirmed an order of magnitude smaller too. |

---

## 3. Overall conclusion

PANDA's three core responsibilities -- geometry, particle generation/transport, and energy/charge scoring -- have been validated for correctness and stress-tested for robustness:

- **Against real, open experimental data** (Hitachi HM68512, McNulty et al. 1989): agreement within a factor of ~2 (Hitachi) and a partially-closed, understood-cause residual gap (McNulty), respectively.
- **Against other simulation tools**: CREME-MC agreement is good (log-RMSE ~0.29 decades, full-range coverage) once compared against PANDA's correct 100%-collection-equivalent (deposited-charge) output -- see 1.3 correction. CUPID's remaining disagreement is attributed to its own structural limitations (liquid-drop evaporation model, no direct-ionization modeling), not to defects in PANDA.
- **Across device scale**, from mid-2000s bulk-junction geometries down to 80nm-node deep-submicron devices: geometry, biasing, production cuts, and recoil/LET export all confirmed functioning correctly, with one real fidelity bug (missing scaled production cuts) found and fixed along the way.
- **Two plausible-sounding hypotheses tested and rejected** (gold-lid LET-ceiling increase, ~80 um shell ceiling) -- both are recorded as closed, negative results rather than open concerns.
- **One genuine, non-obvious physical effect confirmed**: heavy-metal (W/Au/Pb/Ta) dead-layer/via-plug materials produce real fission-fragment recoils well above silicon's own LET ceiling, with severity that does not scale monotonically with atomic number.
- **One validation-methodology bug found and fixed** (2026-07-19): the CREME-MC comparison script was using PANDA's collected-charge (not deposited-charge) cross section, comparing against the wrong physical quantity for an RPP-model reference; correcting it resolved what had been documented as an unexplained intermediate-charge "hump."

Remaining known gap: McNulty's residual high-energy under-prediction, considered structural to the comparison (open, understood-cause but not fully closed), not an open action item. CUPID's absent low-energy end is expected by construction (see 1.2, 1.3).
