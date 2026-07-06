import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# =========================
# McNulty et al. data
# =========================
# S. El-Teleaty, P.J. McNulty, W.G. Abdel-Kader, W.J. Beauvais,
# "Soft Fails in Microelectronic Circuits Due to Proton-Induced Nuclear
# Reactions in Material Surrounding the SEU-Sensitive Volume," Nucl.
# Instr. and Meth. in Physics Research B40/41, 1300-1305 (1989) is the
# paper this comparison was originally requested against. Its own data
# table/figure wasn't accessible (paywalled, no open-access copy
# found). Instead this uses a closely related, freely accessible
# result from the same research group and the same physical quantity:
# a McNulty-authored review chapter reproduced in a DTIC technical
# report (AD-A285581, Defense Nuclear Agency DNA-TR-92-163), which
# contains -- across two of its appendices -- both:
#   - Fig. 12 / Fig. 35: integral cross section (for depositing at
#     least energy E) versus E, for the IDT 6116V NMOS SRAM exposed to
#     148 MeV protons -- both the real pulse-height measurement
#     ("measured") and McNulty's own CUPID Monte Carlo simulation.
#   - The sensitive-volume geometry CUPID assumed: an 11.7 x 11.7 x
#     13.5 um silicon parallelepiped (stated explicitly in the report
#     text, confirmed self-consistent with the ~137 um^2 junction area
#     estimate given nearby).
#
# Data points were digitized directly from the scanned figure (report
# page 303): the log-scale y-axis and linear x-axis were calibrated by
# locating the actual axis-label ink pixel-by-pixel in a high-
# resolution render, then both curves were traced column-by-column,
# filtering out the horizontal SEU-cross-section reference line and
# two vertical threshold-marker lines the figure also draws. The two
# curves are coincident below ~20 MeV (hence one shared table below)
# and diverge above it. Treat this as a digitization with real but
# modest pixel-level uncertainty, same spirit as the CREME-MC table in
# compare_creme_panda.py.
#
# Columns: [Energy Deposited (MeV), Integral Cross Section (cm^2)]
mcnulty_shared = np.array([
    [1, 4.5771e-11], [2, 3.4643e-11], [3, 2.5826e-11], [4, 1.4352e-11],
    [5, 1.3575e-11], [6, 9.9670e-12], [7, 7.8957e-12], [8, 5.8861e-12],
    [9, 4.3880e-12], [10, 3.4410e-12], [11, 2.5011e-12], [12, 1.9614e-12],
    [13, 1.4996e-12], [14, 1.0736e-12], [15, 8.1669e-13], [16, 6.4697e-13],
    [19, 2.2111e-13], [20, 1.8803e-13],
])

mcnulty_measured_tail = np.array([
    [22, 1.0188e-13], [23, 7.3306e-14], [24, 5.5766e-14], [25, 5.0908e-14],
    [26, 3.5174e-14], [27, 2.8006e-14], [28, 2.1413e-14], [29, 1.5883e-14],
    [30, 1.2582e-14], [31, 1.1840e-14], [32, 1.1142e-14],
])

mcnulty_cupid_tail = np.array([
    [22, 5.8072e-14], [23, 3.9922e-14], [24, 2.8006e-14], [25, 2.0151e-14],
    [26, 1.1721e-14], [27, 6.7830e-15], [28, 4.2783e-15], [29, 1.4401e-15],
    [30, 9.7017e-16], [31, 9.7017e-16], [32, 9.7017e-16],
])

mcnulty_measured = np.vstack([mcnulty_shared, mcnulty_measured_tail])
mcnulty_cupid = np.vstack([mcnulty_shared, mcnulty_cupid_tail])

# =========================
# PANDA data
# =========================
# Run via Macros/run_mcnulty_148MeV.mac: 148 MeV protons, 11.7x11.7x13.5
# um silicon sensitive volume (matching CUPID's geometry above), 0.1 um
# dead layer (McNulty's model has no separate dead layer -- this is as
# close to "none" as the two-volume geometry allows), cross-section
# biasing factor 1041 (same validated approach used throughout this
# repo) to reach the rare-tail statistics needed. Archived at
# Results/McNulty_148MeV_validated/events.csv.
#
# McNulty's plot is cross section vs. RAW DEPOSITED ENERGY (not
# collected charge), so this reads Total_keV directly rather than
# either charge column -- no pair-creation-energy conversion needed,
# avoiding an extra assumption. The weighting/cross-section formula
# mirrors PANDA_Analyze.py exactly (weighted histogram, cumulative
# from the top, sigma = beam_area * P(E >= threshold)).
script_dir = os.path.dirname(os.path.abspath(__file__))
events_path = os.path.join(
    script_dir, "Results", "McNulty_148MeV_validated", "events.csv"
)

nParticles = 10_000_000
beamXY_um = 11.7
beam_area_cm2 = (beamXY_um * 1e-4) ** 2

df = pd.read_csv(events_path, usecols=["Total_keV", "EventWeight"])
E_MeV = df["Total_keV"].values / 1000.0
w = df["EventWeight"].values

positive = E_MeV[E_MeV > 0]
bins = np.logspace(np.log10(positive.min()), np.log10(positive.max()), 200)
bin_centers = np.sqrt(bins[:-1] * bins[1:])

hist, _ = np.histogram(E_MeV, bins=bins, weights=w)
cum = np.cumsum(hist[::-1])[::-1]
prob = cum / nParticles
sigma_panda = beam_area_cm2 * prob

x_panda = bin_centers
y_panda = sigma_panda

# =========================
# Interpolate PANDA onto McNulty's x-grid (log-log space)
# =========================
finite_mask = y_panda > 0
log_x_panda = np.log10(x_panda[finite_mask])
log_y_panda = np.log10(y_panda[finite_mask])

interp_panda_log = interp1d(
    log_x_panda,
    log_y_panda,
    bounds_error=False,
    fill_value=(log_y_panda[0], -300.0)
)


def panda_at(E):
    return 10 ** interp_panda_log(np.log10(E))


# =========================
# Error metrics
# =========================
# Same statistics-limited-tail treatment as compare_creme_panda.py:
# restrict RMSE to the range where PANDA has real (non-exhausted)
# statistics, since a handful of raw weighted events at the highest
# energies produce noisy, not-meaningfully-comparable tail values.
x_panda_min, x_panda_max = x_panda[finite_mask].min(), x_panda[finite_mask].max()

print("===== COMPARISON METRICS =====")
for name, table in [("measured", mcnulty_measured), ("CUPID", mcnulty_cupid)]:
    in_range = (table[:, 0] >= x_panda_min) & (table[:, 0] <= x_panda_max)
    y_p = panda_at(table[in_range, 0])
    y_m = table[in_range, 1]

    log_rmse = np.sqrt(np.mean((np.log10(y_p) - np.log10(y_m)) ** 2))
    ratio = y_m / y_p

    print(f"\nvs {name} ({in_range.sum()} points in PANDA's range "
          f"[{x_panda_min:.2f}, {x_panda_max:.2f}] MeV):")
    print(f"  Log RMSE: {log_rmse:.3f} decades")
    print(f"  McNulty/PANDA ratio: {ratio.min():.2f}x to {ratio.max():.2f}x")

# =========================
# Plot
# =========================
fig, ax = plt.subplots(figsize=(9, 7))

ax.loglog(x_panda[finite_mask], y_panda[finite_mask], label="PANDA", linewidth=2)
ax.loglog(
    mcnulty_measured[:, 0], mcnulty_measured[:, 1],
    'o-', label="McNulty et al. -- measured (IDT 6116V, 148 MeV p)",
    markersize=4, color="black"
)
ax.loglog(
    mcnulty_cupid[:, 0], mcnulty_cupid[:, 1],
    's--', label="McNulty et al. -- CUPID simulation",
    markersize=4, color="gray"
)

ax.set_xlabel("Energy Deposited (MeV)")
ax.set_ylabel("Integral Cross Section (cm$^2$)")
ax.set_title(
    "PANDA vs. McNulty et al. (1991) -- 148 MeV protons, "
    "11.7x11.7x13.5 um Si"
)
ax.grid(True, which="both", alpha=0.3)
ax.legend()

plt.tight_layout()
plt.savefig(os.path.join(script_dir, "PANDA_vs_McNulty.png"), dpi=300)
plt.show()
