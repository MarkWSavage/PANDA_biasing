import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# =========================
# CREME-MC data
# =========================
# NOTE: this table was read directly from a CREME-MC report where the
# x-axis was labeled "Energy (MeV)" -- i.e. deposited energy, NOT
# collected charge. It must be converted to charge (fC) using the same
# electron-hole pair creation energy convention used in the Geant4
# side (EventAction.cc / SteppingAction.cc: 3.6 eV per e-h pair in
# silicon) before it can be compared against PANDA's charge-threshold
# cross section curve. Comparing these raw MeV values directly against
# PANDA's fC values (as before) produces a systematic ~44.5x offset
# across the entire curve, not just the tail.
creme_data_energy_MeV = np.array([
    [0.0032696788286249545, 8.451543126112921e-7],
    [0.004147999108728953, 6.571876506768537e-7],
    [0.005251311599776524, 4.2911570045761557e-7],
    [0.00650126293121954, 2.535001231232242e-7],
    [0.00787021593701818, 1.5051507007089045e-7],
    [0.009634984707238541, 8.803910020465874e-8],
    [0.0120628742029542, 5.196700387531566e-8],
    [0.015272182366951163, 3.217342488788875e-8],
    [0.019553294583076697, 1.996742076327664e-8],
    [0.02503509969436659, 1.206804455079583e-8],
    [0.031343747446839224, 7.090814071688742e-9],
    [0.03752059677353492, 4.216449019897554e-9],
    [0.042946053644889216, 2.4282736043677004e-9],
    [0.047000745472734355, 1.3775807438855723e-9],
    [0.050018444273199736, 7.474699711406815e-10],
    [0.05216279742303146, 4.3694381400352067e-10],
    [0.05570029110405535, 2.3101932654691272e-10],
    [0.060063284312186806, 1.1191515199915118e-10],
    [0.06598315112100978, 6.047935388613379e-11],
    [0.07167880281644161, 3.2419665662664215e-11],
    [0.08644247198024499, 2.0721529447851275e-11],
    [0.11063158525419894, 1.9565188515652598e-11],
    [0.14158836926003576, 1.8637293576698272e-11],
    [0.181209638571392, 1.7519704860851085e-11],
    [0.2319197015711352, 1.636037635721731e-11],
    [0.29682176932259036, 1.521043044006078e-11],
    [0.37988803243432445, 1.4078987669434996e-11],
    [0.48620156596285213, 1.3002959437976288e-11],
    [0.6222724985551648, 1.190354743499111e-11],
    [0.7964346238795016, 1.075364600306422e-11],
    [1.0193475690533105, 9.65067365509507e-12],
    [1.304685601012328, 8.415707615968428e-12],
    [1.6699332380899023, 7.1626388162940525e-12],
    [2.137505693613958, 5.871512637925175e-12],
    [2.736134297136543, 4.554579871680787e-12],
    [3.5025207331693298, 3.4179009254943353e-12],
    [4.48399687717985, 2.312015817025102e-12],
    [5.740966619216466, 1.4317141959685357e-12],
    [7.028215016273092, 8.466751955457716e-13],
    [8.135179684653055, 4.887338512022942e-13],
    [9.52381953338351, 2.4729546138954216e-13],
    [11.022771726850628, 1.5882963856700206e-13],
    [12.337232514162446, 8.865786361339795e-14],
    [13.503316812254672, 4.5344104906146894e-14],
    [14.616137767429812, 2.1026662181649995e-14],
    [16.613634411373226, 9.97101058949741e-15],
    [19.147753322237982, 7.367774232974523e-15],
    [19.157127146352746, 4.317409216447429e-15],
    [19.16650555944337, 2.529939402708838e-15],
    [19.307776917002048, 1.290799031726123e-15]
])

# Silicon e-h pair creation energy and elementary charge (matches the
# constants used in EventAction.cc / SteppingAction.cc).
EH_PAIR_ENERGY_EV = 3.6
ELEMENTARY_CHARGE_C = 1.602176634e-19

# Convert MeV -> fC: pairs = E_eV / 3.6, charge_C = pairs * e, /1e-15 for fC
MEV_TO_FC = (1.0e6 / EH_PAIR_ENERGY_EV) * ELEMENTARY_CHARGE_C / 1.0e-15

creme_data = creme_data_energy_MeV.copy()
creme_data[:, 0] = creme_data_energy_MeV[:, 0] * MEV_TO_FC

print(f"Converted CREME-MC x-axis: Energy (MeV) -> Charge (fC), factor = {MEV_TO_FC:.4f} fC/MeV")

# =========================
# PANDA data
# =========================
# Resolved relative to this script's own location (not a hardcoded
# path to a sibling project) so this always reads whichever repo it
# actually lives in -- see PANDA_GUI.py's project_dir for the same fix.
script_dir = os.path.dirname(os.path.abspath(__file__))
panda_file = os.path.join(
    script_dir, "Results", "Current", "cumulative_cross_section_collected.csv"
)
panda = pd.read_csv(panda_file)

# Adjust column names if needed
x_panda = panda["ChargeThreshold_fC"].values
y_panda = panda["CrossSection_cm2"].values

x_creme = creme_data[:, 0]
y_creme = creme_data[:, 1]

# =========================
# Interpolate PANDA onto CREME x-grid (log-log space)
# =========================
# Cross-section-vs-threshold curves are smooth in log-log, not
# linear-linear, so interpolate/extrapolate in log space and convert
# back. This keeps any extrapolation near the data edges physically
# reasonable instead of producing linear-space artifacts.
#
# PANDA's cumulative cross section can legitimately hit exactly 0.0
# at charge thresholds above the highest value any simulated particle
# actually produced (no events, no probability). log10(0) = -inf, and
# feeding that into interp1d silently produces nan once you interpolate
# or extrapolate near it -- which then poisons the entire RMSE average
# through a single bad point. Filter those out before taking logs, and
# use explicit finite fill values for anything outside the real data
# range instead of raw log-space "extrapolate".
finite_mask = y_panda > 0
log_x_panda = np.log10(x_panda[finite_mask])
log_y_panda = np.log10(y_panda[finite_mask])

interp_panda_log = interp1d(
    log_x_panda,
    log_y_panda,
    bounds_error=False,
    # Below the lowest observed threshold: cross section is flat at
    # its maximum (matches the physical plateau seen in both curves).
    # Above the highest observed threshold: no particles produced that
    # much charge, so cross section is effectively zero (use a very
    # negative log value rather than -inf to stay finite).
    fill_value=(log_y_panda[0], -300.0)
)

y_panda_interp = 10 ** interp_panda_log(np.log10(x_creme))

# =========================
# Error metrics
# =========================
# CREME's tail extends slightly beyond the maximum charge PANDA's own
# simulation ever produced (a statistics-limited effect, not a bug --
# see the differential spectrum discussion). Comparing a real CREME
# value there against an artificially-floored PANDA extrapolation
# produces a meaningless, RMSE-dominating log-error from a handful of
# points. Restrict the RMSE to the region where PANDA has real data,
# and report the excluded tail separately instead of silently
# distorting the metric.
x_panda_min = x_panda[finite_mask].min()
x_panda_max = x_panda[finite_mask].max()
in_range = (x_creme >= x_panda_min) & (x_creme <= x_panda_max)
n_excluded = np.sum(~in_range)

rmse_linear = np.sqrt(np.mean((y_panda_interp[in_range] - y_creme[in_range])**2))

rmse_log = np.sqrt(
    np.mean(
        (np.log10(y_panda_interp[in_range]) - np.log10(y_creme[in_range]))**2
    )
)

print("\n===== COMPARISON METRICS =====")
print("(CREME-MC x-axis converted from Energy (MeV) to Charge (fC))")
print(f"RMSE computed over {np.sum(in_range)}/{len(x_creme)} CREME points "
      f"within PANDA's simulated charge range [{x_panda_min:.4g}, {x_panda_max:.4g}] fC")
if n_excluded > 0:
    print(f"  ({n_excluded} CREME point(s) beyond PANDA's max simulated charge "
          f"excluded -- statistics-limited tail, not a defect)")
print(f"Linear RMSE : {rmse_linear:.6e}")
print(f"Log RMSE    : {rmse_log:.6e} decades")

# =========================
# Plot
# =========================
fig, axs = plt.subplots(1, 2, figsize=(14, 6))

# Log-log
axs[0].loglog(
    x_panda, y_panda,
    label="PANDA",
    linewidth=2
)

axs[0].loglog(
    x_creme, y_creme,
    'o-',
    label="CREME-MC",
    markersize=4
)

axs[0].set_xlabel("Collected Charge (fC)")
axs[0].set_ylabel("Cross Section")
axs[0].set_title("Log-Log Comparison")
axs[0].grid(True, which="both")
axs[0].legend()

# Lin-log
axs[1].semilogy(
    x_panda, y_panda,
    label="PANDA",
    linewidth=2
)

axs[1].semilogy(
    x_creme, y_creme,
    'o-',
    label="CREME-MC",
    markersize=4
)

axs[1].set_xlabel("Collected Charge (fC)")
axs[1].set_ylabel("Cross Section")
axs[1].set_title("Lin-Log Comparison")
axs[1].grid(True, which="both")
axs[1].legend()

# RMSE text box on lower-right of lin-log panel
rmse_text = (
    f"Linear RMSE = {rmse_linear:.3e}\n"
    f"Log RMSE = {rmse_log:.3f} decades"
)

axs[1].text(
    0.97, 0.05,
    rmse_text,
    transform=axs[1].transAxes,
    fontsize=10,
    verticalalignment='bottom',
    horizontalalignment='right',
    bbox=dict(
        boxstyle="round",
        facecolor="white",
        alpha=0.85
    )
)

plt.tight_layout()
plt.savefig(os.path.join(script_dir, "PANDA_vs_CREME.png"), dpi=300)
plt.show()
