# PANDA

**PANDA** — *Protons And Neutron charge Deposition in mAterials* — is a Geant4 Monte Carlo simulation for calculating charge deposition from energetic particles (protons, neutrons, ions) in semiconductor structures, for Single Event Effect (SEE) analysis.

The core physical quantity computed is **collected charge** per event; upset probability, cross-section spectra, and other SEE metrics are derived from it in post-processing (`PANDA_Analyze.py`), keeping particle transport and interpretation cleanly separated.

Compared against CREME-MC (see `compare_creme_panda.py`). Includes optional cross-section biasing (`SEEBiasingOperator`) to efficiently sample the rare nuclear-recoil tail, with correct raw-charge/event-weight separation so biased and unbiased runs produce statistically consistent spectra.

## Known limitations

- **CREME-MC comparison shoulder**: PANDA's cross-section curve sits ~1-2 orders of magnitude above CREME-MC's in the ~2-20 fC charge range (log-RMSE ~0.9 decades over the full curve), reproducible at both low and high statistics (i.e. not a sampling artifact). This is expected: PANDA (Geant4, QGSP_BIC_HP) and CREME-MC use different underlying nuclear reaction model families -- Geant4's cascade/pre-compound physics vs. CREME-MC's semi-empirical fragmentation cross-sections -- and published comparisons of Geant4 against CREME96 report exactly this kind of physics-list-dependent discrepancy (see [arXiv:0712.2149](https://arxiv.org/pdf/0712.2149)). Not considered fixable in PANDA's code; treat sub-decade-scale disagreement in this charge range as expected rather than a bug.

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
