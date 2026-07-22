import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker


# Resolved relative to this script's own location (not the caller's
# cwd) -- same fix applied to PANDA_GUI.py and compare_creme_panda.py
# this session.
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(SCRIPT_DIR, "Results", "Current")
os.makedirs(RESULTS_DIR, exist_ok=True)

RECOIL_HITS_CSV = os.path.join(RESULTS_DIR, "recoil_hits.csv")


# ----------------------------
# Load PANDA output
# ----------------------------
csv_file = os.path.join(RESULTS_DIR, "events.csv")

print("Loading:", csv_file)

data = pd.read_csv(csv_file)

# Weight every histogram/statistic by EventWeight, same as
# PANDA_Analyze.py -- required so a biased run's rare, heavily-
# down-weighted tail events don't get counted as full raw hits.
# Falls back to weight=1.0 for events.csv files predating this column.
has_weight_column = "EventWeight" in data.columns

if has_weight_column:
    weight = data["EventWeight"].values
    print("EventWeight column found -- building bias-corrected (weighted) histograms.")
else:
    weight = np.ones(len(data))
    print("No EventWeight column found (older events.csv) -- assuming weight=1.0 for all events.")

total_keV = data["Total_keV"].values

proton_keV   = data["Proton_keV"].values
electron_keV = data["Electron_keV"].values
recoil_keV   = data["Recoil_keV"].values

# PrimaryIon_keV splits out the primary beam particle's own track
# (TrackID==1/ParentID==0) from Recoil_keV -- for any primary species
# that isn't proton/e- (alpha and every heavier ion), this is the
# dominant, expected continuous electronic-stopping deposit (SRIM's
# "Ions" curve), not a nuclear recoil. Falls back to zero for
# events.csv files predating this column, so the Recoil curve below
# just reverts to its old (primary+recoil-conflated) meaning rather
# than erroring.
has_primary_ion_column = "PrimaryIon_keV" in data.columns

if has_primary_ion_column:
    primary_ion_keV = data["PrimaryIon_keV"].values
else:
    primary_ion_keV = np.zeros(len(data))
    print("No PrimaryIon_keV column found (older events.csv) -- Recoil "
          "still includes any primary-ion-track deposit, as before this split.")

# ----------------------------
# Filter Data
# ----------------------------
nonzero_mask = total_keV > 0
nonzero_keV    = total_keV[nonzero_mask]
nonzero_weight = weight[nonzero_mask]

proton_mask      = proton_keV > 0
electron_mask    = electron_keV > 0
recoil_mask      = recoil_keV > 0
primary_ion_mask = primary_ion_keV > 0

proton_nonzero,   proton_weight   = proton_keV[proton_mask],     weight[proton_mask]
electron_nonzero, electron_weight = electron_keV[electron_mask], weight[electron_mask]
recoil_nonzero,   recoil_weight   = recoil_keV[recoil_mask],     weight[recoil_mask]
primary_ion_nonzero, primary_ion_weight = (
    primary_ion_keV[primary_ion_mask], weight[primary_ion_mask]
)

# ----------------------------
# Basic statistics
# ----------------------------
N_total = len(total_keV)

print("\nEnergy statistics:")
print(f"Min:  {np.min(total_keV):.3f} keV")
# Weighted mean = sum(w*q) / N, N = number of primaries simulated (NOT
# sum of weights) -- standard importance-sampling estimator, matching
# PANDA_Analyze.py.
print(f"Mean: {np.sum(total_keV * weight) / N_total:.3f} keV")
print(f"Max:  {np.max(total_keV):.3f} keV")
print(f"Total events: {N_total}")
print(f"Triggered events (raw count): {len(nonzero_keV)}")
print(f"Efficiency (weighted): {np.sum(nonzero_weight) / N_total:.8e}")

if has_primary_ion_column:
    total_sum = np.sum(total_keV * weight)
    print("\nComponent share of Total_keV (weighted):")
    print(f"  Proton:      {np.sum(proton_keV * weight) / total_sum * 100:.2f}%")
    print(f"  Electron:    {np.sum(electron_keV * weight) / total_sum * 100:.2f}%")
    print(f"  Primary ion: {np.sum(primary_ion_keV * weight) / total_sum * 100:.2f}%")
    print(f"  Recoil:      {np.sum(recoil_keV * weight) / total_sum * 100:.2f}%")

# ----------------------------
# Histogram
# ----------------------------
# Log-spaced bins, not linear: post-biasing the data spans several
# orders of magnitude (common few-keV ionization deposits vs. rare
# MeV-scale nuclear-recoil tail), and linear bins from 0 to the single
# largest observed value stretch so wide that the bulk of the
# distribution collapses into the first bin or two. Range is taken
# from the smallest and largest nonzero value across ALL FIVE series
# (not just Total) since a component can be nonzero-but-small in a row
# where other components dominate the row's Total.
all_nonzero_keV = np.concatenate([
    nonzero_keV, proton_nonzero, electron_nonzero, primary_ion_nonzero, recoil_nonzero
])

bins = np.logspace(
    np.log10(np.min(all_nonzero_keV)),
    np.log10(np.max(all_nonzero_keV)),
    200
)

hist_total, edges = np.histogram(
    nonzero_keV,
    bins=bins,
    weights=nonzero_weight
)

hist_proton, _ = np.histogram(
    proton_nonzero,
    bins=bins,
    weights=proton_weight
)

hist_electron, _ = np.histogram(
    electron_nonzero,
    bins=bins,
    weights=electron_weight
)

hist_recoil, _ = np.histogram(
    recoil_nonzero,
    bins=bins,
    weights=recoil_weight
)

hist_primary_ion, _ = np.histogram(
    primary_ion_nonzero,
    bins=bins,
    weights=primary_ion_weight
)

# Geometric mean, not arithmetic, is the correct bin-center convention
# for log-spaced bins (matches PANDA_Analyze.py's bin_centers).
centers = np.sqrt(edges[:-1] * edges[1:])


# ----------------------------
# Plot
# ----------------------------
plt.figure(figsize=(10, 6))

plt.step(
    centers,
    hist_total,
    where="mid",
    label="Total"
)

plt.step(
    centers,
    hist_proton,
    where="mid",
    label="Proton"
)

plt.step(
    centers,
    hist_electron,
    where="mid",
    label="Electron"
)

plt.step(
    centers,
    hist_primary_ion,
    where="mid",
    label="Primary ion (own track)"
)

plt.step(
    centers,
    hist_recoil,
    where="mid",
    label="Recoil (genuine secondary)"
)

plt.legend()
plt.xlabel("Deposited Energy (keV)")
plt.ylabel("Weighted Counts" if has_weight_column else "Counts")
plt.title("PANDAEX Raw Deposited Energy Spectrum")
plt.xscale("log")
plt.yscale("log")

outfile = os.path.join(
    RESULTS_DIR,
    "PANDAEX_component_spectrum.png"
)

plt.savefig(outfile, dpi=300)

print("\nSaved:")
print(outfile)

# ----------------------------
# Per-hit recoil species/LET breakdown (recoil_hits.csv, opt-in --
# see /sim/logRecoilHits). The Recoil curve above lumps every non-
# proton/e- species into one bucket; this is the finer-grained view
# that actually distinguishes recoil species and reports LET, which
# events.csv's per-event aggregate can't (it only sums energy).
# ----------------------------
if not os.path.exists(RECOIL_HITS_CSV):
    print("\nrecoil_hits.csv not found -- per-species/LET breakdown unavailable.")
    print("Run with /sim/logRecoilHits true to enable it (see README")
    print("'Per-hit recoil/LET export').")
else:
    print("\nLoading:", RECOIL_HITS_CSV)
    hits = pd.read_csv(RECOIL_HITS_CSV)

    # SteppingAction/EventAction's export filter excludes proton/e- (not
    # true nuclear recoils) but not e+ -- a positron shows up from the
    # same reactions (beta-plus decay/pair production) without being a
    # nuclear recoil either. Filtered here rather than at the source so
    # the raw CSV still has it available if ever wanted for something
    # else; every view below (summary table, "All recoils" aggregate,
    # per-species LET curves) should be genuine recoils only.
    non_recoil_species = {"e+"}
    hits = hits[~hits["Particle"].isin(non_recoil_species)]

    hit_weight = hits["EventWeight"].values
    let = hits["LET_MeV_cm2_mg"].values
    species = hits["Particle"].values

    print(f"Recoil hits: {len(hits)}")

    # Per-species summary: hit count (raw and weighted), Z/A, and
    # max/mean LET -- the mean is bias-corrected (weighted by
    # EventWeight, matching every other statistic in this codebase);
    # the max is NOT weighted, since it's asking "what LET was ever
    # actually reached", not an expectation value.
    def _species_row(g):
        w = g["EventWeight"].values
        let = g["LET_MeV_cm2_mg"].values
        return pd.Series({
            "Count": len(g),
            "WeightedCount": w.sum(),
            "Z": g["Z"].iloc[0],
            "A": g["A"].iloc[0],
            "MaxLET_MeV_cm2_mg": let.max(),
            "MeanLET_MeV_cm2_mg": np.average(let, weights=w),
        })

    summary = (
        hits.groupby("Particle")
            .apply(_species_row, include_groups=False)
            .sort_values("Count", ascending=False)
    )

    summary_path = os.path.join(RESULTS_DIR, "PANDAEX_recoil_species_summary.csv")
    summary.to_csv(summary_path)

    print("\nRecoil species summary (by hit count):")
    print(summary.to_string())
    print("\nSaved:", summary_path)

    # Differential LET spectrum: overall + the top species by hit
    # count (the rest would clutter the legend without adding much,
    # since hit count falls off steeply after the first few species --
    # see the summary CSV above for the complete per-species table).
    TOP_N_SPECIES = 6
    top_species = summary.head(TOP_N_SPECIES).index.tolist()

    let_positive = let[let > 0]

    # Floor the spectrum's range at 1e-4 MeV*cm2/mg: a handful of steps
    # (mostly deuteron/light-fragment msc/Transportation steps with
    # near-zero edep) have real but physically uninteresting LET down
    # to ~1e-19, which stretches the log-scale x-axis across 20+
    # decades and squashes the actual nuclear-recoil peak into a
    # sliver. Values below the floor are outside the histogram range
    # and are simply not counted (np.histogram default behavior).
    LET_PLOT_FLOOR = 1e-4

    let_bins = np.logspace(
        np.log10(max(let_positive.min(), LET_PLOT_FLOOR)),
        np.log10(let_positive.max()),
        100
    )
    let_centers = np.sqrt(let_bins[:-1] * let_bins[1:])

    plt.figure(figsize=(10, 6))

    hist_all_let, _ = np.histogram(let, bins=let_bins, weights=hit_weight)
    plt.step(let_centers, hist_all_let, where="mid", label="All recoils", linewidth=2, color="black")

    for sp in top_species:
        mask = species == sp
        h, _ = np.histogram(let[mask], bins=let_bins, weights=hit_weight[mask])
        plt.step(let_centers, h, where="mid", label=sp)

    plt.xlabel("LET (MeV cm$^2$/mg)")
    plt.ylabel("Weighted Counts")
    plt.title("PANDAEX Recoil LET Spectrum (sensitive volume)")
    plt.xscale("log")
    plt.yscale("log")
    plt.legend()

    let_outfile = os.path.join(RESULTS_DIR, "PANDAEX_LET_spectrum.png")
    plt.savefig(let_outfile, dpi=300)

    print("\nSaved:", let_outfile)

    # Small-multiples breakdown: every recoil species gets its own panel.
    # Overlaying all of them on one axes would mean cycling through 100+
    # generated hues once past the top few -- that reads as noise, not
    # signal, so each species gets its own subplot (identified by its
    # title, not a color) instead of fighting for a spot in one legend.
    all_species = summary.index.tolist()
    n_species = len(all_species)
    n_cols = 12
    n_rows = int(np.ceil(n_species / n_cols))

    SPECIES_LINE_COLOR = "#2a78d6"  # dataviz reference palette, categorical slot 1

    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(n_cols * 1.6, n_rows * 1.3),
        sharex=True,
    )
    axes_flat = axes.flatten()

    for ax, sp in zip(axes_flat, all_species):
        mask = species == sp
        h, _ = np.histogram(let[mask], bins=let_bins, weights=hit_weight[mask])
        ax.step(let_centers, h, where="mid", linewidth=1, color=SPECIES_LINE_COLOR)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(let_bins[0], let_bins[-1])
        # Cap major ticks so sparse species (few nonzero bins spanning
        # under a decade) don't pile up overlapping log-minor-tick labels.
        ax.yaxis.set_major_locator(mticker.LogLocator(base=10.0, numticks=3))
        ax.yaxis.set_minor_locator(mticker.NullLocator())
        ax.set_title(sp, fontsize=7, pad=2)
        ax.tick_params(labelsize=5)

    # Small-multiples convention: only the outer row/column carry tick
    # labels, so 108 panels don't each repeat the same axis text.
    for ax in axes_flat[:n_species]:
        ax.label_outer()

    for ax in axes_flat[n_species:]:
        ax.axis("off")

    fig.suptitle(
        "PANDAEX Recoil LET Spectrum -- all species (small multiples)",
        fontsize=12
    )
    fig.supxlabel("LET (MeV cm$^2$/mg)", fontsize=9)
    fig.supylabel("Weighted Counts", fontsize=9)
    fig.tight_layout(rect=[0, 0, 1, 0.97])

    small_multiples_outfile = os.path.join(
        RESULTS_DIR, "PANDAEX_LET_spectrum_all_species.png"
    )
    fig.savefig(small_multiples_outfile, dpi=150)
    plt.close(fig)

    print("Saved:", small_multiples_outfile)

print("\nDone.")

plt.show()
