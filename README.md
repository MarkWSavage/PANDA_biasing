# PANDA

**Current release: [v1.1.0](https://github.com/MarkWSavage/PANDA/releases/tag/v1.1.0) — Material Selection.** Sensitive volume and dead layer/electrode are now independently selectable materials (see "Materials" below), documented as a scoped exception to the [v1.0.0](https://github.com/MarkWSavage/PANDA/releases/tag/v1.0.0) Architecture Freeze; see `Documentation/PANDA_MASTER_DESIGN` (Section 5, Version History) for what's frozen and what's in scope.

**PANDA** — *Protons And Neutron charge Deposition in mAterials* — is a Geant4 Monte Carlo simulation for calculating charge deposition from energetic particles (protons, neutrons, ions) in semiconductor structures, for Single Event Effect (SEE) analysis.

The core physical quantity computed is **collected charge** per event; upset probability, cross-section spectra, and other SEE metrics are derived from it in post-processing (`PANDA_Analyze.py`), keeping particle transport and interpretation cleanly separated.

Compared against CREME-MC (see `compare_creme_panda.py`) and against McNulty et al.'s published proton-induced experimental/CUPID data (see `compare_mcnulty_panda.py`). Includes optional cross-section biasing (`SEEBiasingOperator`) to efficiently sample the rare nuclear-recoil tail, with correct raw-charge/event-weight separation so biased and unbiased runs produce statistically consistent spectra.

## Materials

The sensitive volume and dead layer/electrode each have an independently selectable material, set via `/sim/sensitiveMaterial` and `/sim/deadMaterial` in a macro, or the matching dropdowns in `PANDA_GUI.py` (next to the Sensitive/Dead Thickness fields):

- **Sensitive volume**: Si, GaAs, Ge, SiC, GaN
- **Dead layer/electrode**: SiO2, Al2O3, TiO2, Si

The sensitive volume's material also selects the pair-creation energy and carrier mobility/saturation velocity used to convert deposited energy to charge -- see `DetectorConstruction::GetSensitivePairCreationEnergy()` and neighbors for the per-material constants and sources. Defaults are Si/Si, matching PANDA's original silicon-only behavior.

## Known limitations

- **CREME-MC comparison shoulder**: PANDA's cross-section curve sits ~1-2 orders of magnitude above CREME-MC's in the ~2-20 fC charge range (log-RMSE ~0.9 decades over the full curve), reproducible at both low and high statistics (i.e. not a sampling artifact). This is expected: PANDA (Geant4, QGSP_BIC_HP) and CREME-MC use different underlying nuclear reaction model families -- Geant4's cascade/pre-compound physics vs. CREME-MC's semi-empirical fragmentation cross-sections -- and published comparisons of Geant4 against CREME96 report exactly this kind of physics-list-dependent discrepancy (see [arXiv:0712.2149](https://arxiv.org/pdf/0712.2149)). Not considered fixable in PANDA's code; treat sub-decade-scale disagreement in this charge range as expected rather than a bug.

- **McNulty et al. comparison under-prediction**: comparing PANDA's integral cross-section-vs-deposited-energy curve against McNulty et al.'s 148 MeV proton data for an 11.7x11.7x13.5 um silicon sensitive volume (`compare_mcnulty_panda.py`) shows PANDA systematically under-predicting, growing smoothly from ~1.4x at 1 MeV deposited to >100x by 25-28 MeV (log-RMSE ~1.1 decades vs. the experimental curve). Root cause: PANDA's sensitive+dead-layer stack sits in a vacuum world, with no material outside it for a nearby nuclear reaction to occur in. McNulty's own CUPID model explicitly included a silicon "surrounding volume" for this reason -- it's the subject of the 1989 paper this comparison was originally aimed at (El-Teleaty, McNulty, Abdel-Kader, Beauvais, Nucl. Instr. and Meth. B40/41, 1300-1305). The discrepancy growing with deposited energy (where surrounding-volume-sourced recoils matter increasingly more) rather than being flat is consistent with this being the missing mechanism. Not fixable without adding a surrounding-material volume to `DetectorConstruction` -- tracked as a known geometry limitation, not a nuclear-physics-list defect like the CREME-MC discrepancy above.

## Build

```
mkdir -p build && cd build
cmake ..
make -j$(nproc)
```

## Run

```
./PANDA run.mac
python3 PANDA_Analyze.py
```

See `Documentation/PANDA_MASTER_DESIGN` for the full design philosophy.
