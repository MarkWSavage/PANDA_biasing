import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
import os

RESULTS_DIR = os.path.join("Results", "Current")
os.makedirs(RESULTS_DIR, exist_ok=True)

# -----------------------------
# Parse run.mac
# -----------------------------
particle = None
energy = None
beamXY = None
sensitiveXY = None
sensitiveThickness = None
deadThickness = None
incidentAngle = None
nParticles = None
criticalCharge = None
useCollectionModel = None

with open("run.mac", "r") as f:
    for line in f:
        if "/sim/particle" in line:
            particle = line.split()[-1]

        elif "/sim/energy" in line:
            m = re.search(r'/sim/energy\s+([\d\.]+)', line)
            if m:
                energy = float(m.group(1))

        elif "/sim/beamXY" in line:
            m = re.search(r'/sim/beamXY\s+([\d\.]+)', line)
            if m:
                beamXY = float(m.group(1))

        elif "/sim/sensitiveXY" in line:
            m = re.search(r'/sim/sensitiveXY\s+([\d\.]+)', line)
            if m:
                sensitiveXY = float(m.group(1))

        elif "/sim/sensitiveThickness" in line:
            m = re.search(r'/sim/sensitiveThickness\s+([\d\.]+)', line)
            if m:
                sensitiveThickness = float(m.group(1))

        elif "/sim/deadThickness" in line:
            m = re.search(r'/sim/deadThickness\s+([\d\.]+)', line)
            if m:
                deadThickness = float(m.group(1))

        elif "/sim/incidentAngle" in line:
            m = re.search(r'/sim/incidentAngle\s+([\d\.]+)', line)
            if m:
                incidentAngle = float(m.group(1))

        elif "/sim/criticalCharge" in line:
            m = re.search(r'/sim/criticalCharge\s+([\d\.]+)', line)
            if m:
                criticalCharge = float(m.group(1))

        elif "/sim/useCollectionModel" in line:
            useCollectionModel = line.split()[-1].strip().lower() == "true"

        elif "/run/beamOn" in line:
            m = re.search(r'/run/beamOn\s+(\d+)', line)
            if m:
                nParticles = int(m.group(1))

beam_area_cm2 = (beamXY * 1e-4)**2
fluence = nParticles / beam_area_cm2

print("\nRun conditions:")
print("Particle:", particle)
print("Energy (MeV):", energy)
print("Beam XY (um):", beamXY)
print("Sensitive XY (um):", sensitiveXY)
print("Sensitive thickness (um):", sensitiveThickness)
print("Dead layer (um):", deadThickness)
print("Incident angle (deg):", incidentAngle if incidentAngle is not None else 0.0)
print("Particles:", nParticles)
print("Critical charge (fC):", criticalCharge)
print("Upset criterion uses:",
      "CollectedCharge_fC" if useCollectionModel else "DepositedCharge_fC")

# -----------------------------
# Metrics to analyze
# -----------------------------
# Both charge quantities are always present in the CSV now, so we always
# analyze both rather than guessing which one the user "meant". This
# keeps the analysis in sync with the simulation regardless of how
# /sim/useCollectionModel was set for a given run.
METRICS = {
    "DepositedCharge_fC": {
        "label": "Deposited Charge",
        "suffix": "deposited",
    },
    "CollectedCharge_fC": {
        "label": "Collected Charge",
        "suffix": "collected",
    },
}

events_csv_path = os.path.join(RESULTS_DIR, "events.csv")
csv_header = pd.read_csv(events_csv_path, nrows=0).columns.tolist()
has_weight_column = "EventWeight" in csv_header

# -----------------------------
# First pass: find the true charge range
# -----------------------------
# A fixed 1e-3..1e3 fC bin range worked for the small default geometry
# this was tuned against, but silently breaks for a larger sensitive
# volume: np.histogram DROPS values outside [bins[0], bins[-1]]
# entirely rather than clipping them into the edge bins, so any run
# depositing more than 1000 fC (e.g. a bigger sensitiveXY/thickness)
# had those events vanish from every histogram, cumulative curve, and
# CSV output -- invisible unless you noticed the printed Max exceeded
# 1000 fC. Bins are now sized from the real observed range instead.
global_min = {key: np.inf for key in METRICS}
global_max = {key: -np.inf for key in METRICS}

for chunk in pd.read_csv(
    events_csv_path,
    usecols=list(METRICS.keys()),
    chunksize=1000000
):
    for key in METRICS:
        q = chunk[key].values
        positive = q[q > 0]
        if positive.size:
            global_min[key] = min(global_min[key], np.min(positive))
        if q.size:
            global_max[key] = max(global_max[key], np.max(q))

# Shared bins across both metrics -- keeps the Deposited/Collected
# plots and CSVs on a directly comparable charge-threshold grid, same
# as the previous single hardcoded `bins`.
range_min = min(v for v in global_min.values() if np.isfinite(v))
range_max = max(v for v in global_max.values() if np.isfinite(v))

bins = np.logspace(np.log10(range_min), np.log10(range_max), 400)
bin_centers = np.sqrt(bins[:-1] * bins[1:])

stats = {
    key: {
        "hist": np.zeros(len(bins) - 1),
        "min": np.inf,
        "max": -np.inf,
        "sum": 0.0,
        "count": 0,
    }
    for key in METRICS
}

# -----------------------------
# Second pass: stream CSV, updating all metrics together
# -----------------------------

if has_weight_column:
    usecols = list(METRICS.keys()) + ["EventWeight"]
    print("EventWeight column found -- building bias-corrected (weighted) histograms.")
else:
    usecols = list(METRICS.keys())
    print("No EventWeight column found (older events.csv) -- assuming weight=1.0 for all events.")

for chunk in pd.read_csv(
    events_csv_path,
    usecols=usecols,
    chunksize=1000000
):
    if has_weight_column:
        w = chunk["EventWeight"].values
    else:
        w = np.ones(len(chunk))

    for key in METRICS:
        q = chunk[key].values

        # Weighted histogram: each event contributes its own weight
        # (1.0 unless biasing is active) rather than a flat count of 1.
        # This is what keeps the differential spectrum and cumulative
        # cross section statistically correct when cross-section
        # biasing (SEEBiasingOperator) has been used -- otherwise the
        # artificially boosted rare events would be overcounted by
        # roughly the bias factor.
        h, _ = np.histogram(q, bins=bins, weights=w)
        s = stats[key]
        s["hist"] += h
        s["min"] = min(s["min"], np.min(q))
        s["max"] = max(s["max"], np.max(q))
        # Weighted sum for an unbiased mean estimate: mean = sum(w*q) / N,
        # where N is the actual number of primaries simulated (NOT the
        # sum of weights) -- standard importance-sampling estimator.
        s["sum"] += np.sum(q * w)
        s["count"] += len(q)

# -----------------------------
# Plots + CSV output, per metric
# -----------------------------
for key, meta in METRICS.items():
    label = meta["label"]
    suffix = meta["suffix"]
    s = stats[key]
    hist = s["hist"]
    countQ = s["count"]

    # Cumulative probability (computed here, ahead of the differential
    # spectrum plot below, so its plateau-derived x-axis cutoff -- see
    # plot_xmin below -- can be reused for that plot too).
    cum_counts = np.cumsum(hist[::-1])[::-1]
    prob = cum_counts / countQ

    # P(Q>=q) is monotonically non-increasing as the threshold rises, so
    # its max is always prob[0] -- but a single stray near-zero charge
    # value (e.g. a boundary-tolerance artifact, more likely to show up
    # at deep-submicron device scale) can drag the bin range down many
    # extra decades below where the real distribution starts, stretching
    # a flat P=1 plateau across most of the plot. Find where the curve
    # actually leaves that plateau and start the plotted (not the
    # underlying data/CSV -- both still cover the full range) x-axis
    # there instead of at the true bin minimum. Reused below for the
    # differential spectrum too, since the same stray values also
    # stretch that plot's x-axis out to the same uninformative range.
    plateau_val = prob[0]
    still_plateau = np.isclose(prob, plateau_val, rtol=1e-9, atol=0)
    plot_xmin = None if still_plateau.all() else bin_centers[np.where(still_plateau)[0][-1]]

    # Differential spectrum
    plt.figure(figsize=(8, 6))
    plt.step(bin_centers, hist, where='mid')
    plt.xscale("log")
    plt.yscale("log")
    if plot_xmin is not None:
        plt.xlim(left=plot_xmin)
    plt.xlabel(f"{label} (fC)")
    plt.ylabel("Counts")
    plt.title(f"PANDA Differential Charge Spectrum ({label})")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(
        os.path.join(RESULTS_DIR, f"differential_charge_spectrum_{suffix}.png")
    )
    plt.close()

    plt.figure(figsize=(8, 6))
    plt.plot(bin_centers, prob)
    plt.xscale("log")
    plt.yscale("log")
    if plot_xmin is not None:
        plt.xlim(left=plot_xmin)
    plt.xlabel(f"{label} Threshold (fC)")
    plt.ylabel("P(Q ≥ q)")
    plt.title(f"PANDA Cumulative Probability ({label})")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(
        os.path.join(RESULTS_DIR, f"cumulative_probability_{suffix}.png")
    )
    plt.close()

    # Cross section
    sigma = beam_area_cm2 * prob

    plt.figure(figsize=(8, 6))
    plt.plot(bin_centers, sigma)
    plt.xscale("log")
    plt.yscale("log")
    if plot_xmin is not None:
        plt.xlim(left=plot_xmin)
    plt.xlabel(f"{label} Threshold (fC)")
    plt.ylabel("Cross Section (cm²)")
    plt.title(f"PANDA Cumulative Cross Section ({label})")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(
        os.path.join(RESULTS_DIR, f"cumulative_cross_section_{suffix}.png")
    )
    plt.close()

    out = pd.DataFrame({
        "ChargeThreshold_fC": bin_centers,
        "Probability": prob,
        "CrossSection_cm2": sigma
    })

    out.to_csv(
        os.path.join(RESULTS_DIR, f"cumulative_cross_section_{suffix}.csv"),
        index=False
    )

    print(f"\n[{label}] statistics:")
    print("  Min:", s["min"], "fC")
    print("  Mean:", s["sum"] / countQ, "fC")
    print("  Max:", s["max"], "fC")

    if criticalCharge is not None:
        # Interpolate P(Q >= criticalCharge) from the cumulative curve
        p_at_qc = np.interp(criticalCharge, bin_centers, prob)
        print(f"  P(Q >= {criticalCharge} fC):", p_at_qc)

print("\nSaved to:", RESULTS_DIR)
for meta in METRICS.values():
    suffix = meta["suffix"]
    print(f"  differential_charge_spectrum_{suffix}.png")
    print(f"  cumulative_probability_{suffix}.png")
    print(f"  cumulative_cross_section_{suffix}.png")
    print(f"  cumulative_cross_section_{suffix}.csv")

print("\nDone.")
