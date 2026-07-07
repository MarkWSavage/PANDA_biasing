# PANDA

**Current release: [v1.4.0](https://github.com/MarkWSavage/PANDA/releases/tag/v1.4.0) — Per-Species Recoil LET Analysis.** `PANDAEX_Analyze.py` now analyzes the per-hit recoil/LET output added in [v1.3.0](https://github.com/MarkWSavage/PANDA/releases/tag/v1.3.0) (see "Per-hit recoil/LET export" below), producing a per-species summary table and a differential LET spectrum plot. See `Documentation/PANDA_MASTER_DESIGN` (Section 5, Version History) for what's frozen and what's in scope.

**PANDA** — *Protons And Neutron charge Deposition in mAterials* — is a Geant4 Monte Carlo simulation for calculating charge deposition from energetic particles (protons, neutrons, ions) in semiconductor structures, for Single Event Effect (SEE) analysis.

The core physical quantity computed is **collected charge** per event; upset probability, cross-section spectra, and other SEE metrics are derived from it in post-processing (`PANDA_Analyze.py`), keeping particle transport and interpretation cleanly separated.

Compared against CREME-MC (see `compare_creme_panda.py`) and against McNulty et al.'s published proton-induced experimental/CUPID data (see `compare_mcnulty_panda.py`). Includes optional cross-section biasing (`SEEBiasingOperator`) to efficiently sample the rare nuclear-recoil tail, with correct raw-charge/event-weight separation so biased and unbiased runs produce statistically consistent spectra. The dead+sensitive stack sits inside a configurable bulk "surrounding volume" (matching the sensitive material, `/sim/surroundingXY`/`/sim/surroundingThickness`), with optional biasing of secondary neutrons produced there too (`/sim/secondaryNeutronBiasFactor`), so nearby nuclear reactions can contribute recoils into the sensitive volume.

## Materials

The sensitive volume and dead layer/electrode each have an independently selectable material, set via `/sim/sensitiveMaterial` and `/sim/deadMaterial` in a macro, or the matching dropdowns in `PANDA_GUI.py` (next to the Sensitive/Dead Thickness fields):

- **Sensitive volume**: Si, GaAs, Ge, SiC, GaN
- **Dead layer/electrode**: SiO2, Al2O3, TiO2, Si

The sensitive volume's material also selects the pair-creation energy and carrier mobility/saturation velocity used to convert deposited energy to charge -- see `DetectorConstruction::GetSensitivePairCreationEnergy()` and neighbors for the per-material constants and sources. Defaults are Si/Si, matching PANDA's original silicon-only behavior.

## Per-hit recoil/LET export

Set `/sim/logRecoilHits true` (default `false`) to additionally write `Results/Current/recoil_hits.csv`, one row per energy-depositing hit in the sensitive volume (species, Z, A, LET in MeV·cm²/mg, position, EventWeight), filtered to recoils only (excludes `proton`/`e-`). This is a finer-grained companion to `events.csv`'s per-event `Proton_keV`/`Electron_keV`/`Recoil_keV` sums -- useful for studying the recoil-species/LET spectrum directly, e.g. checking where the highest LET a given recoil species reaches actually lands (a common reference point: silicon recoils are often cited as topping out around LET~12, informing the assumption that heavy-ion hardness above LET=20 implies proton/neutron-SEE immunity). Off by default since per-step file writes add real overhead most runs don't need.

Run `python3 PANDAEX_Analyze.py` to analyze it: a per-species summary table (count, weighted count, Z, A, max/mean LET, saved as `PANDAEX_recoil_species_summary.csv`) and a differential LET spectrum plot (`PANDAEX_LET_spectrum.png`, overall plus the top species by hit count). Falls back gracefully if `recoil_hits.csv` doesn't exist.

## Known limitations

- **CREME-MC comparison shoulder**: PANDA's cross-section curve sits ~1-2 orders of magnitude above CREME-MC's in the ~2-20 fC charge range (log-RMSE ~0.9 decades over the full curve), reproducible at both low and high statistics (i.e. not a sampling artifact). This is expected: PANDA (Geant4, QGSP_BIC_HP) and CREME-MC use different underlying nuclear reaction model families -- Geant4's cascade/pre-compound physics vs. CREME-MC's semi-empirical fragmentation cross-sections -- and published comparisons of Geant4 against CREME96 report exactly this kind of physics-list-dependent discrepancy (see [arXiv:0712.2149](https://arxiv.org/pdf/0712.2149)). Not considered fixable in PANDA's code; treat sub-decade-scale disagreement in this charge range as expected rather than a bug.

- **McNulty et al. comparison under-prediction (partially addressed)**: comparing PANDA's integral cross-section-vs-deposited-energy curve against McNulty et al.'s 148 MeV proton data for an 11.7x11.7x13.5 um silicon sensitive volume (`compare_mcnulty_panda.py`) originally showed PANDA systematically under-predicting, growing smoothly from ~1.4x at 1 MeV deposited to >100x by 25-28 MeV (log-RMSE ~1.1 decades), traced to the dead+sensitive stack sitting in a vacuum world with no material for nearby nuclear reactions to occur in -- the actual subject of the 1989 paper this comparison was aimed at (El-Teleaty, McNulty, Abdel-Kader, Beauvais, Nucl. Instr. and Meth. B40/41, 1300-1305). Adding a 100 um silicon surrounding volume plus biasing secondary neutrons within it (silicon's cm-scale neutron mean free path meant unbiased secondary neutrons almost never reacted again in any practical volume) improved log-RMSE to ~0.91 decades -- a real ~18% reduction, now comparable to the CREME-MC discrepancy above rather than distinctly worse -- but a large gap remains at high deposited energy (e.g. still ~48x at 25 MeV, down from ~133x). The remaining gap isn't yet understood: candidates include the 100 um surrounding volume still being too small, only one cascade generation being biased (not a second-generation secondary neutron), and other secondary species (protons, alphas) not being biased. See `Documentation/PANDA_MASTER_DESIGN` Known Limitations for the full before/after.

## Build

```
mkdir -p build && cd build
cmake ..
make -j$(nproc)
```

## Run

From the project root (not from inside `build/`):

```
./build/PANDA run.mac
python3 PANDA_Analyze.py
```

See `Documentation/PANDA_MASTER_DESIGN` for the full design philosophy.
