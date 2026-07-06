import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# Resolved relative to this script's own location (not the caller's
# cwd) -- same fix applied to PANDA_GUI.py and compare_creme_panda.py
# this session.
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(SCRIPT_DIR, "Results", "Current")
os.makedirs(RESULTS_DIR, exist_ok=True)


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

# ----------------------------
# Filter Data
# ----------------------------
nonzero_mask = total_keV > 0
nonzero_keV    = total_keV[nonzero_mask]
nonzero_weight = weight[nonzero_mask]

proton_mask   = proton_keV > 0
electron_mask = electron_keV > 0
recoil_mask   = recoil_keV > 0

proton_nonzero,   proton_weight   = proton_keV[proton_mask],     weight[proton_mask]
electron_nonzero, electron_weight = electron_keV[electron_mask], weight[electron_mask]
recoil_nonzero,   recoil_weight   = recoil_keV[recoil_mask],     weight[recoil_mask]

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

# ----------------------------
# Histogram
# ----------------------------
bins = np.linspace(
    0,
    np.max(nonzero_keV),
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

centers = 0.5 * (edges[:-1] + edges[1:])


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
    hist_recoil,
    where="mid",
    label="Recoil"
)

plt.legend()
plt.xlabel("Deposited Energy (keV)")
plt.ylabel("Weighted Counts" if has_weight_column else "Counts")
plt.title("PANDAEX Raw Deposited Energy Spectrum")
plt.yscale("log")

outfile = os.path.join(
    RESULTS_DIR,
    "PANDAEX_component_spectrum.png"
)

plt.savefig(outfile, dpi=300)

print("\nSaved:")
print(outfile)

plt.show()
