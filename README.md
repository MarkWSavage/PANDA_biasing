# PANDA

**Current release: [v1.7.14](https://github.com/MarkWSavage/PANDA/releases/tag/v1.7.14) — CREME-MC Thin-Film Direct-Ionization Shoulder Documented.** A dedicated thin-film check (200 MeV proton, 0.5um sensitive Si, 0.5um SiO2 dead layer) found PANDA's cross section running up to ~75x higher than CREME-MC's in a narrow 0.1-0.5 fC window, tapering back to agreement by ~1-2 fC. Root-caused inside PANDA via a per-event `events.csv` breakdown (`Electron_keV`/`Total_keV`): 83% of events in that window carry an explicit delta-ray (knock-on electron) hit, averaging 91% of the event's deposited energy -- a real, correctly-functioning secondary-transport effect, not primary-track straggling or a production-cut artifact (the region-scoped cut here is 50nm, a 10x margin under the 500nm sensitive layer). The likely explanation for the magnitude gap -- CREME-MC's engine (MRED) hands electron transport to PENELOPE below some threshold, vs. PANDA's native Geant4 EM physics throughout -- could not be confirmed against MRED's actual configuration, which isn't publicly documented (the primary source, Reed et al. 2015 IEEE TNS, is paywalled). Documentation only -- no functional or behavioral code changes. See `Documentation/PANDA_VALIDATION_SUMMARY.md` Section 1.5 for the full writeup.

**PANDA** — *Protons And Neutron charge Deposition in mAterials* — is a Geant4 Monte Carlo simulation for calculating charge deposition from energetic particles (protons, neutrons, ions) in semiconductor structures, for Single Event Effect (SEE) analysis.

The core physical quantity computed is **collected charge** per event; upset probability, cross-section spectra, and other SEE metrics are derived from it in post-processing (`PANDA_Analyze.py`), keeping particle transport and interpretation cleanly separated.

Compared against CREME-MC (see `compare_creme_panda.py`) and against McNulty et al.'s published proton-induced experimental/CUPID data (see `compare_mcnulty_panda.py`). Includes optional cross-section biasing (`SEEBiasingOperator`) to efficiently sample the rare nuclear-recoil tail, with correct raw-charge/event-weight separation so biased and unbiased runs produce statistically consistent spectra. The dead+sensitive stack sits inside a configurable bulk "surrounding volume" (matching the sensitive material, `/sim/surroundingXY`/`/sim/surroundingThickness`), with optional biasing of secondary neutrons produced there too (`/sim/secondaryNeutronBiasFactor`), so nearby nuclear reactions can contribute recoils into the sensitive volume.

## Primary particles

Set via `/sim/particle` in a macro, or the Particle dropdown in `PANDA_GUI.py`: `proton`, `neutron`, `alpha`, `deuteron`, `triton`, `He3`, `e-`, or one of six heavy-ion primaries -- `C12`, `F19`, `Cl35`, `Ni58`, `I127`, `Au197` (ground-state, most-abundant-stable-isotope-per-element, matching common heavy-ion SEE test cocktails). Cross-section biasing (`/sim/biasCrossSectionFactor`) is validated for every species except `e-` (see `PANDA.cc`'s wrapped-process list). Unlike SRIM -- which computes a bare ion's stopping power/range but has no device or charge-collection model -- PANDA reports cross-section-vs-LET through its own geometry, charge-collection-efficiency, and angle-of-incidence models for these ions directly, covering the direct-ionization-dominated regime of SEE testing alongside its existing proton/neutron nuclear-recoil focus. See `Documentation/PANDA_MASTER_DESIGN`'s Cross-section biasing exception (Section 5) for the ion-creation ordering constraint this needed and how it was resolved.

**Nuclear-recoil tail beyond SRIM:** SRIM only computes a bare ion's stopping power along its own track, so it has no way to predict the population of secondary recoil fragments produced when a heavy ion undergoes a nuclear reaction (`ionInelastic`) inside the device -- fragments that can carry higher LET than the primary ion's own direct ionization. PANDA's differential charge spectrum shows this directly: a 200 MeV Cl35 beam through a 1 um dead layer + 1 um Si sensitive volume, with `sensitiveXY` set to 0.1 um (representative of a 65 nm SRAM cell's sensitive volume), produces a sharp direct-ionization peak matching the SRIM-predicted LET (2.694E+00 MeV*cm^2/mg -> ~28 fC expected from LET x thickness, matching the observed peak just above 20 fC), plus a separate, much rarer tail extending out past 100 fC from nuclear-recoil events that a SRIM-only hazard assessment would never predict. Despite the initial concern that such a small lateral sensitive volume might be an edge case for the geometry/charge-collection model, the resulting spectrum was well-behaved and physically sensible at this deep-submicron scale. This is the same underlying mechanism (nuclear-reaction-produced recoils, not direct ionization) behind PANDA's proton/neutron SEE-susceptibility work -- heavy ions just add a high-LET direct-ionization peak on top of it.

## Materials

The sensitive volume and dead layer/electrode each have an independently selectable material, set via `/sim/sensitiveMaterial` and `/sim/deadMaterial` in a macro, or the matching dropdowns in `PANDA_GUI.py` (next to the Sensitive/Dead Thickness fields):

- **Sensitive volume**: Si, GaAs, Ge, SiC, GaN
- **Dead layer/electrode**: SiO2, Al2O3, TiO2, Si

The sensitive volume's material also selects the pair-creation energy and carrier mobility/saturation velocity used to convert deposited energy to charge -- see `DetectorConstruction::GetSensitivePairCreationEnergy()` and neighbors for the per-material constants and sources. Defaults are Si/Si, matching PANDA's original silicon-only behavior.

## Per-hit recoil/LET export

Set `/sim/logRecoilHits true` (default `false`) to additionally write `Results/Current/recoil_hits.csv`, one row per energy-depositing hit in the sensitive volume (species, Z, A, LET in MeV·cm²/mg, position, EventWeight), filtered to recoils only (excludes `proton`/`e-`). This is a finer-grained companion to `events.csv`'s per-event `Proton_keV`/`Electron_keV`/`PrimaryIon_keV`/`Recoil_keV` sums -- useful for studying the recoil-species/LET spectrum directly, e.g. checking where the highest LET a given recoil species reaches actually lands (a common reference point: silicon recoils are often cited as topping out around LET~12, informing the assumption that heavy-ion hardness above LET=20 implies proton/neutron-SEE immunity). Off by default since per-step file writes add real overhead most runs don't need. Note this file's per-hit export applies a 100nm minimum step-length floor to avoid LET-inflation artifacts from very short steps (see `SteppingAction.cc`), so for deep-submicron geometries a heavy-ion primary's own continuous track rarely clears it -- `events.csv`'s `PrimaryIon_keV` column (not gated by that floor) is where that contribution is actually visible.

`events.csv`'s `PrimaryIon_keV` column (added alongside the existing `Proton_keV`/`Electron_keV`/`Recoil_keV`) separates the primary beam particle's own track (any species other than proton/e-, i.e. alpha and every heavier ion) from genuine secondary nuclear-recoil deposits, which previously shared one `Recoil_keV` bucket purely because both are "not proton, not e-". For a 200 MeV Au197 primary this is the majority of `Total_keV` (~88%) -- that's the primary ion's own expected electronic stopping (SRIM's "Ions" curve), not a nuclear-recoil effect, and it's why `Recoil_keV` alone should not be read as a recoil-population indicator for non-proton/e- primaries.

Run `python3 PANDAEX_Analyze.py` to analyze it: a per-species summary table (count, weighted count, Z, A, max/mean LET, saved as `PANDAEX_recoil_species_summary.csv`) and a differential LET spectrum plot (`PANDAEX_LET_spectrum.png`, overall plus the top species by hit count). Falls back gracefully if `recoil_hits.csv` doesn't exist.

## Angle of incidence (approximate)

Set `/sim/incidentAngle` (degrees, default `0`) to approximate a beam tilted off the sensitive volume's normal. Rather than rotating the beam or geometry (which would also need to reshape the dead layer, surrounding volume, and secondary-neutron-biasing region), this divides the sensitive volume's effective thickness by `cos(incidentAngle)` -- the same chord-length-elongation model used for tilt-angle corrections in heavy-ion SEE testing (`LET_eff = LET(0 deg)/cos(theta)`), valid for a broad-uniform-beam test condition (beam footprint much larger than the device). The lateral (XY) footprint is intentionally left unchanged -- no physically-grounded correction for it was found; see `Documentation/PANDA_MASTER_DESIGN` for why. Valid range `[0, 90)` deg; anything outside that raises a fatal exception rather than silently producing an infinite or negative thickness.

See `Documentation/PANDA_VALIDATION_SUMMARY.md` for the full validation summary -- comparisons against McNulty et al./CUPID and CREME-MC/MRED (both against independently-specified geometry), plus every geometry/physics robustness check run against PANDA itself. (Hitachi HM68512 open proton-SEU data is also discussed there, but reclassified as a capability demonstration rather than a validation, since its device geometry and critical charge were never published.)

## Known limitations

- ~~**CREME-MC comparison shoulder**~~ (resolved 2026-07-19): the previously reported ~1-2 order of magnitude gap was a comparison-methodology bug in `compare_creme_panda.py` -- it compared CREME-MC's ideal-100%-charge-collection reference against PANDA's *collected*-charge curve rather than the *deposited*-charge curve (PANDA's own 100%-collection-equivalent quantity). Corrected, the two curves agree closely (log-RMSE ~0.29 decades) across the full range, including the region previously described as an unexplained "hump." See `Documentation/PANDA_VALIDATION_SUMMARY.md` Section 1.3.

- ~~**McNulty et al. comparison under-prediction**~~ (reframed 2026-07-19): comparing PANDA's integral cross-section-vs-deposited-energy curve against McNulty et al.'s 148 MeV proton data for an 11.7x11.7x13.5 um silicon sensitive volume (`compare_mcnulty_panda.py`) originally showed PANDA systematically under-predicting, growing smoothly from ~1.4x at 1 MeV deposited to >100x by 25-28 MeV (log-RMSE ~1.1 decades), traced to the dead+sensitive stack sitting in a vacuum world with no material for nearby nuclear reactions to occur in -- the actual subject of the 1989 paper this comparison was aimed at (El-Teleaty, McNulty, Abdel-Kader, Beauvais, Nucl. Instr. and Meth. B40/41, 1300-1305). Adding a 100 um silicon surrounding volume plus biasing secondary neutrons within it improved log-RMSE to ~0.91 decades -- a real ~18% reduction -- but a gap remained at high deposited energy (e.g. ~48x at 25 MeV, down from ~133x). Two candidate causes were tested directly and ruled out: widening the beam to cover the surrounding volume made agreement slightly *worse* (log-RMSE 0.426 -> 0.455 decades over a matched range), and growing the surrounding volume to 300/500 um showed no improving trend (0.900 -> 0.887 -> 0.945 decades). Root cause: heavy recoil fragments (the only secondaries that deposit meaningfully) have very short range, so only reactions happening close to the sensitive volume can contribute -- a structural feature of the physics, not a fixable sampling gap. Combined with CUPID's liquid-drop nuclear model being fundamentally different from PANDA's Geant4 cascade physics (not two approximations of the same model), agreement within an order of magnitude with matching curve shape is judged a genuine validation by the standard used for cross-code nuclear-model comparisons elsewhere in the field, not an open defect. See `Documentation/PANDA_MASTER_DESIGN` Known Limitations for the full before/after.

- **Heavy-ion test energy should stay realistic (>~1 GeV for Au-mass ions)** (found 2026-07-22, mechanism resolved 2026-07-22): 200 MeV Au197 (specific energy ~1 MeV/nucleon) produces an anomalously narrow, near-delta-function differential charge spectrum (relative std 2.96%) rather than genuine Vavilov/Landau straggling. Confirmed by rerunning the identical geometry at 3000 MeV (~15 MeV/nucleon, within real cyclotron heavy-ion SEE test-beam ranges): the spectrum becomes a properly-spread peak (relative std 8.30%), well clear of the max-possible-charge edge. Separate from the mean-LET model swap in v1.7.10/v1.7.11. The mechanism is ordinary Landau/Vavilov collision-counting straggling (Weeden-Wright et al., "Effects of Stopping Ions and LET Fluctuations on Soft Error Rate Prediction," IEEE TNS, 2015 -- from the MRED/CREME-MC developers themselves), not charge-state fluctuation (considered and found unsupported in that source). The specific magnitude is explained by the Vavilov shape parameter kappa ~ z^2/beta^4 (not dE/dx): a ~215x swing in kappa between 200 MeV and 3000 MeV Au197, driven by velocity alone, correctly predicts the observed ~2.8x change in relative spread -- see `Documentation/PANDA_VALIDATION_SUMMARY.md` Section 2 for the full derivation. Not a PANDA-specific bug -- just a reminder to pick physically-achievable test energies.

- **CREME-MC thin-film direct-ionization shoulder, magnitude unexplained** (found 2026-07-23/24): a 200 MeV proton through a 0.5um sensitive Si / 0.5um SiO2 dead layer stack shows PANDA's cross section running up to ~75x higher than CREME-MC's in a narrow 0.1-0.5 fC window, converging again by ~1-2 fC. Confirmed via a per-event `events.csv` breakdown that this window is genuinely delta-ray-dominated inside PANDA (83% of events carry an explicit knock-on-electron hit, averaging 91% of that event's deposited energy) -- not primary-track straggling, and not a production-cut artifact (the region-scoped cut here is 50nm, a 10x margin under the 500nm sensitive layer). CREME-MC's engine (MRED) is itself Geant4-based but hands electron transport to PENELOPE below some threshold rather than using Geant4's native EM models throughout, which is a plausible source of the gap -- but the actual handoff threshold/cutoffs MRED uses aren't publicly documented, so this couldn't be confirmed. No evidence of a PANDA defect; see `Documentation/PANDA_VALIDATION_SUMMARY.md` Section 1.5.

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

## Linting

`ruff` (Python) and `clang-tidy` (C++) checks live in `scripts/pre-commit-checks.sh`.
Git hooks aren't tracked by version control, so each clone needs a one-time setup:

```
ln -sf ../../scripts/pre-commit-checks.sh .git/hooks/pre-commit
python3 -m venv .lint-venv && .lint-venv/bin/pip install ruff
```

clang-tidy also needs `build/compile_commands.json`, generated by building with
`cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=ON ..`. The hook blocks commits on ruff
findings; clang-tidy findings are printed but never block.

See `Documentation/PANDA_MASTER_DESIGN` for the full design philosophy.
