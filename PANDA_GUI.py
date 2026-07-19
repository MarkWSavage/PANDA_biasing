import sys
import os
import subprocess
import json

from PyQt5.QtWidgets import (
    QApplication, QWidget, QLabel, QLineEdit,
    QPushButton, QTextEdit, QVBoxLayout, QHBoxLayout,
    QGridLayout, QFileDialog, QCheckBox,
    QGroupBox, QComboBox, QDialog, QMessageBox
)

from PyQt5.QtGui import QPixmap
from PyQt5.QtCore import Qt


class PandaGUI(QWidget):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("🐼 PANDA SEE Simulator")
        self.setGeometry(100, 100, 1400, 1050)

        self.project_dir = os.path.dirname(os.path.abspath(__file__))
        self.build_dir = os.path.join(self.project_dir, "build")
        self.results_dir = os.path.join(self.project_dir, "Results", "Current")

        self.default_carrier_lifetime = "1.0 ns"
        self.default_electric_field = "1000000 volt/m"

        self.initUI()

    def style_button(self, button, bg, border):
        button.setStyleSheet(f"""
        QPushButton {{
            background-color: {bg};
            border: 2px solid {border};
            border-radius: 12px;
            font-size: 20px;
            font-weight: bold;
            padding: 20px;
        }}
        QPushButton:hover {{
            background-color: white;
        }}
        QPushButton:pressed {{
            background-color: #dddddd;
        }}
        """)

    def initUI(self):
        main_layout = QVBoxLayout()
        top_layout = QHBoxLayout()

        # ==========================
        # Simulation Parameters
        # ==========================
        param_group = QGroupBox("Simulation Parameters")
        param_layout = QGridLayout()

        self.fields = {}

        labels = [
            ("Energy (MeV)", "200"),
            ("Sensitive XY (um)", "10"),
            ("Dead XY (um)", "10"),
            ("Sensitive Thickness (um)", "8"),
            ("Dead Thickness (um)", "5"),
            ("Incident Angle (deg)", "0"),
            ("Beam XY (um)", "10"),
            ("Critical Charge (fC)", "150"),
            ("Bias Cross Section Factor", "1041"),
            ("Events", "1000000")
        ]

        param_layout.addWidget(QLabel("Particle:"), 0, 0)

        self.particle_dropdown = QComboBox()
        self.particle_dropdown.addItems([
            "proton",
            "neutron",
            "alpha",
            "deuteron",
            "triton",
            "He3",
            "e-"
        ])
        self.particle_dropdown.setToolTip(
            "Biasing (Bias Cross Section Factor below) is validated for "
            "proton, neutron, alpha, deuteron, triton, and He3 -- see "
            "PANDA.cc for the exact wrapped process per species. e- has "
            "no biasing wired up (runs unbiased regardless of the factor "
            "field)."
        )
        self.particle_dropdown.setFixedWidth(220)
        self.particle_dropdown.setFixedHeight(34)

        param_layout.addWidget(self.particle_dropdown, 0, 1)

        self.sensitive_material_dropdown = QComboBox()
        self.sensitive_material_dropdown.addItems(["Si", "GaAs", "Ge", "SiC", "GaN"])
        self.sensitive_material_dropdown.setToolTip(
            "Sensitive volume material. Also selects the pair-creation\n"
            "energy and mobility/saturation-velocity constants used to\n"
            "convert deposited energy to charge -- see\n"
            "DetectorConstruction::GetSensitivePairCreationEnergy() etc."
        )
        self.sensitive_material_dropdown.setFixedWidth(110)
        self.sensitive_material_dropdown.setFixedHeight(34)

        self.dead_material_dropdown = QComboBox()
        self.dead_material_dropdown.addItems(
            ["SiO2", "Al2O3", "TiO2", "Si", "Au", "W", "Pb", "Ta"]
        )
        self.dead_material_dropdown.setToolTip(
            "Dead-layer/electrode material. Independent of the sensitive\n"
            "volume's material. Au/W/Pb/Ta model metallization/via-plug\n"
            "materials directly against the sensitive volume -- these can\n"
            "produce neutron/proton-induced fission-fragment recoils with\n"
            "LET well above silicon's ~14-15 MeV*cm2/mg ceiling."
        )
        self.dead_material_dropdown.setFixedWidth(110)
        self.dead_material_dropdown.setFixedHeight(34)

        # Sensitive/Dead Thickness rows get a material dropdown appended
        # after the input box (column 2), keeping every row's input box
        # in the same column (1) so they all stay aligned. Every other
        # row is unaffected -- QGridLayout tolerates rows using
        # different column counts.
        material_dropdown_for_label = {
            "Sensitive Thickness (um)": self.sensitive_material_dropdown,
            "Dead Thickness (um)": self.dead_material_dropdown,
        }

        field_tooltips = {
            "Incident Angle (deg)":
                "Approximates a beam tilted off the sensitive volume's "
                "normal by dividing the effective sensitive thickness by "
                "cos(angle) -- the chord-length-elongation model used for "
                "tilt-angle corrections in heavy-ion SEE testing. Valid "
                "range [0, 90). Lateral (XY) footprint is NOT adjusted.",
        }

        for i, (label, default) in enumerate(labels, start=1):
            param_layout.addWidget(QLabel(label + ":"), i, 0)

            field = QLineEdit(default)
            field.setFixedWidth(220)
            field.setFixedHeight(34)

            tooltip = field_tooltips.get(label)
            if tooltip is not None:
                field.setToolTip(tooltip)

            self.fields[label] = field
            param_layout.addWidget(field, i, 1)

            dropdown = material_dropdown_for_label.get(label)
            if dropdown is not None:
                param_layout.addWidget(dropdown, i, 2)

        self.vis_checkbox = QCheckBox("Enable Visualization")
        param_layout.addWidget(self.vis_checkbox, len(labels)+1, 0, 1, 2)

        self.collection_checkbox = QCheckBox("Use Charge Collection Model for Upset Criterion")
        self.collection_checkbox.setChecked(True)
        self.collection_checkbox.setToolTip(
            "Both Deposited and Collected charge are always scored and analyzed.\n"
            "This only selects which one counts toward an upset in the simulation\n"
            "(the UpsetCharge_fC column in events.csv)."
        )
        param_layout.addWidget(self.collection_checkbox, len(labels)+2, 0, 1, 2)

        self.link_xy_checkbox = QCheckBox("Link Beam/Dead XY to Sensitive XY")
        self.link_xy_checkbox.setChecked(True)
        self.link_xy_checkbox.setToolTip(
            "Common case: beam spot, dead-layer footprint, and sensitive-\n"
            "volume footprint all match. Uncheck to set Beam XY and Dead XY\n"
            "independently -- useful for guard rings (dead area wider than\n"
            "the active pixel) or a beam spot that doesn't match the pixel."
        )
        self.link_xy_checkbox.toggled.connect(self.on_link_xy_toggled)
        param_layout.addWidget(self.link_xy_checkbox, len(labels)+3, 0, 1, 2)

        self.log_recoil_hits_checkbox = QCheckBox("Log Recoil Hits (recoil_hits.csv)")
        self.log_recoil_hits_checkbox.setToolTip(
            "Export one row per energy-depositing hit in the sensitive\n"
            "volume (species, Z, A, LET in MeV*cm2/mg, position,\n"
            "EventWeight), filtered to recoils only (excludes proton/e-).\n"
            "Off by default -- per-step file writes add real overhead\n"
            "most runs don't need. See README 'Per-hit recoil/LET export'."
        )
        param_layout.addWidget(self.log_recoil_hits_checkbox, len(labels)+4, 0, 1, 2)

        self.fields['Sensitive XY (um)'].textChanged.connect(self.on_sensitive_xy_changed)
        self.on_link_xy_toggled(self.link_xy_checkbox.isChecked())

        param_group.setLayout(param_layout)

        # ==========================
        # Controls
        # ==========================
        control_group = QGroupBox("Controls")
        control_layout = QGridLayout()

        plot_metric_label = QLabel("Plot Metric:")
        plot_metric_label.setFixedHeight(20)
        control_layout.addWidget(plot_metric_label, 0, 0, 1, 2)

        self.metric_dropdown = QComboBox()
        self.metric_dropdown.addItems(["Deposited", "Collected"])
        self.metric_dropdown.setFixedHeight(34)
        self.metric_dropdown.setToolTip(
            "Both charge metrics are always scored and analyzed.\n"
            "Choose which one the plot buttons below should open."
        )
        control_layout.addWidget(self.metric_dropdown, 1, 0, 1, 2)

        self.run_button = QPushButton("▶ Run PANDA")
        self.style_button(self.run_button, "#d8f5d1", "#6dbb63")
        self.run_button.clicked.connect(self.run_panda)
        control_layout.addWidget(self.run_button, 2, 0)

        self.analyze_button = QPushButton("≣ Analyze Results")
        self.style_button(self.analyze_button, "#d9ecff", "#5b9bd5")
        self.analyze_button.clicked.connect(self.run_analysis)
        control_layout.addWidget(self.analyze_button, 2, 1)

        self.plot_spec_button = QPushButton("▙ Plot Spectrum")
        self.style_button(self.plot_spec_button, "#fff0c9", "#d9a400")
        self.plot_spec_button.clicked.connect(
            lambda: self.open_file(
                f"differential_charge_spectrum_{self.metric_suffix()}.png"
            )
        )
        control_layout.addWidget(self.plot_spec_button, 3, 0)

        self.plot_sigma_button = QPushButton("Σ Plot σ(Q)")
        self.style_button(self.plot_sigma_button, "#eadbff", "#8b5cf6")
        self.plot_sigma_button.clicked.connect(
            lambda: self.open_file(
                f"cumulative_cross_section_{self.metric_suffix()}.png"
            )
        )
        control_layout.addWidget(self.plot_sigma_button, 3, 1)

        self.save_button = QPushButton("■ Save Preset")
        self.style_button(self.save_button, "#f5f5f5", "#999999")
        self.save_button.clicked.connect(self.save_preset)
        control_layout.addWidget(self.save_button, 4, 0)

        self.load_button = QPushButton("▣ Load Preset")
        self.style_button(self.load_button, "#f5f5f5", "#999999")
        self.load_button.clicked.connect(self.load_preset)
        control_layout.addWidget(self.load_button, 4, 1)

        self.open_results_button = QPushButton("▤ Open Results Folder")
        self.style_button(self.open_results_button, "#fff9c4", "#d9c400")
        self.open_results_button.clicked.connect(self.open_results_folder)
        control_layout.addWidget(self.open_results_button, 5, 0, 1, 2)

        self.run_pandaex_button = QPushButton("▤ Run PANDAEX")
        self.style_button(self.run_pandaex_button, "#d9ecff", "#5b9bd5")
        self.run_pandaex_button.clicked.connect(self.run_pandaex)
        control_layout.addWidget(self.run_pandaex_button, 6, 0, 1, 2)

        self.exit_button = QPushButton("✕ Exit PANDA")
        self.style_button(self.exit_button, "#ffd9d9", "#e35c5c")
        self.exit_button.clicked.connect(self.close)
        control_layout.addWidget(self.exit_button, 7, 0, 1, 2)

        # Send all leftover vertical space (this group box is shorter
        # than Simulation Parameters, which has more rows) into a
        # dedicated stretch row instead of letting QGridLayout spread
        # it evenly across every row above -- that's what was making
        # the "Plot Metric" label's row look absurdly tall.
        control_layout.setRowStretch(8, 1)

        control_group.setLayout(control_layout)

        top_layout.addWidget(param_group)
        top_layout.addWidget(control_group)

        main_layout.addLayout(top_layout)

        # ==========================
        # Bottom Half
        # ==========================
        bottom_layout = QHBoxLayout()

        # Console
        console_group = QGroupBox("Console / Log")
        console_layout = QVBoxLayout()

        self.console = QTextEdit()
        self.console.setReadOnly(True)
        self.console.setMouseTracking(True)
        self.console.viewport().setMouseTracking(True)

        console_layout.addWidget(self.console)
        console_group.setLayout(console_layout)

        bottom_layout.addWidget(console_group, 2)

        # Panda image
        self.panda_label = QLabel()
        panda_path = os.path.join(
            self.project_dir,
            "Assets",
            "panda.png"
        )

        if os.path.exists(panda_path):
            pixmap = QPixmap(panda_path)
            pixmap = pixmap.scaled(
                600, 500,
                Qt.KeepAspectRatio,
                Qt.SmoothTransformation
            )
            self.panda_label.setPixmap(pixmap)
        else:
            self.panda_label.setText("Place panda.png in project directory.")

        self.panda_label.setAlignment(Qt.AlignCenter)
        self.panda_label.setAttribute(Qt.WA_TransparentForMouseEvents)

        bottom_layout.addWidget(self.panda_label, 2)

        main_layout.addLayout(bottom_layout)

        self.setLayout(main_layout)

    def log(self, text):
        self.console.append(text)
        QApplication.processEvents()

    def metric_suffix(self):
        return "collected" if self.metric_dropdown.currentText() == "Collected" else "deposited"

    def on_sensitive_xy_changed(self, text):
        if self.link_xy_checkbox.isChecked():
            self.fields['Dead XY (um)'].setText(text)
            self.fields['Beam XY (um)'].setText(text)

    def on_link_xy_toggled(self, checked):
        self.fields['Dead XY (um)'].setEnabled(not checked)
        self.fields['Beam XY (um)'].setEnabled(not checked)

        if checked:
            self.on_sensitive_xy_changed(self.fields['Sensitive XY (um)'].text())

    def write_run_mac(self):
        run_mac_path = os.path.join(self.project_dir, "run.mac")

        with open(run_mac_path, "w") as f:
            f.write("/random/setSeeds 12345 67890\n")

            # DetectorConstruction owns these commands, and (unlike
            # serial mode) Geant4 MT doesn't construct PrimaryGenerator
            # Action/EventAction -- and therefore doesn't register
            # their /sim/energy, /sim/beamXY, /sim/criticalCharge,
            # /sim/verbose commands -- until /run/initialize actually
            # runs. DetectorConstruction's own commands (this block)
            # must come first since they affect geometry construction.
            f.write(f"/sim/particle {self.particle_dropdown.currentText()}\n")
            f.write(f"/sim/sensitiveXY {self.fields['Sensitive XY (um)'].text()} um\n")
            f.write(f"/sim/deadXY {self.fields['Dead XY (um)'].text()} um\n")
            f.write(f"/sim/sensitiveThickness {self.fields['Sensitive Thickness (um)'].text()} um\n")
            f.write(f"/sim/deadThickness {self.fields['Dead Thickness (um)'].text()} um\n")
            f.write(f"/sim/incidentAngle {self.fields['Incident Angle (deg)'].text()} deg\n")
            f.write(f"/sim/sensitiveMaterial {self.sensitive_material_dropdown.currentText()}\n")
            f.write(f"/sim/deadMaterial {self.dead_material_dropdown.currentText()}\n")
            f.write(f"/sim/biasCrossSectionFactor {self.fields['Bias Cross Section Factor'].text()}\n")

            if self.collection_checkbox.isChecked():
                f.write("/sim/useCollectionModel true\n")
            else:
                f.write("/sim/useCollectionModel false\n")

            f.write(f"/sim/carrierLifetime {self.default_carrier_lifetime}\n")
            f.write(f"/sim/electricField {self.default_electric_field}\n")

            f.write("/run/initialize\n")
            f.write("/run/reinitializeGeometry\n\n")

            # PrimaryGeneratorAction/EventAction-owned commands: must
            # come after /run/initialize in MT (see comment above).
            f.write(f"/sim/energy {self.fields['Energy (MeV)'].text()} MeV\n")
            f.write(f"/sim/beamXY {self.fields['Beam XY (um)'].text()} um\n")
            f.write(f"/sim/criticalCharge {self.fields['Critical Charge (fC)'].text()}\n")
            f.write("/sim/verbose false\n")

            if self.log_recoil_hits_checkbox.isChecked():
                f.write("/sim/logRecoilHits true\n")

            f.write("\n")

            if self.vis_checkbox.isChecked():
                f.write("/vis/open OGLSX\n")

                # /vis/open alone does not create a scene -- only
                # /vis/drawVolume implicitly does, later. /vis/scene/add/extent
                # below needs one to already exist, so create it explicitly
                # here; /vis/drawVolume then just adds the volume to it.
                f.write("/vis/scene/create\n")

                # Fixed camera extent scaled off total (dead + sensitive)
                # thickness, NOT sensitiveXY/deadXY. Without this,
                # /vis/drawVolume's auto-fit zoom frames the whole World
                # volume, whose XY extent is auto-grown with
                # sensitiveXY/deadXY (see DetectorConstruction::Construct())
                # while its Z extent stays pinned near a fixed
                # surroundingThickness floor -- so the same physical
                # thickness renders at a wildly different, shrinking pixel
                # scale as sensitiveXY grows (confirmed empirically: a
                # 10 um vs. 1000 um sensitiveXY run at the same 1 um/10 um
                # sensitive/dead thickness produces near-identical
                # DepositedCharge_fC stats, so the *geometry* is unaffected
                # -- only the auto-fit render scale is). Framing a window
                # sized off thickness alone keeps the stack's apparent
                # scale consistent across sensitiveXY/deadXY, at the cost
                # of not showing a large device's full lateral footprint
                # in one view.
                total_thickness_um = (
                    float(self.fields['Sensitive Thickness (um)'].text())
                    + float(self.fields['Dead Thickness (um)'].text())
                )

                # The stack isn't centered at world z=0: the sensitive
                # volume is centered on z=0, but the dead layer hangs
                # below it (see Construct()'s DeadLayer placement), so
                # the stack's true center sits deadThickness/2 below 0.
                dead_thickness_um = float(self.fields['Dead Thickness (um)'].text())
                z_center_um = -dead_thickness_um / 2.0

                half_xy_um = 3.0 * total_thickness_um
                half_z_um = total_thickness_um

                f.write(
                    "/vis/scene/add/extent "
                    f"{-half_xy_um} {half_xy_um} "
                    f"{-half_xy_um} {half_xy_um} "
                    f"{z_center_um - half_z_um} {z_center_um + half_z_um} um\n"
                )

                # Side-on viewpoint (mostly along +X, slightly elevated
                # in Z): the dead layer and sensitive volume share the
                # same XY footprint and are stacked along Z, so a
                # face-on view down Z (the default) shows only the
                # sensitive volume's front face -- the dead layer is
                # directly behind it, fully occluded. This angle shows
                # both layers as a visible cross-section stack instead.
                f.write("/vis/viewer/set/viewpointVector 1 0 0.3\n")
                f.write("/vis/drawVolume\n")
                f.write("/vis/viewer/set/style wireframe\n")

                # Without these, MT Geant4 still stores every event's
                # trajectory for the vis sub-thread (even though nothing
                # here ever asked for trajectories), and OGLSX draws
                # each one far slower than the run itself produces them
                # (~5-6x, measured). The vis queue then backs up past
                # its default cap (100 events) and beamOn blocks waiting
                # on the vis thread -- for a large event count that
                # stalls the whole run behind what looks like a static
                # geometry view, i.e. "visualization but it never runs".
                # add/trajectories + accumulate makes the viewer show a
                # sample of real tracks instead of just the wireframe;
                # actionOnEventQueueFull discard is what actually keeps
                # beamOn running at full speed once that sample is full.
                f.write("/vis/scene/add/trajectories smooth\n")
                f.write("/vis/scene/endOfEventAction accumulate 100\n")
                f.write("/vis/multithreading/actionOnEventQueueFull discard\n")
                f.write("/vis/viewer/flush\n\n")

            f.write(f"/run/beamOn {self.fields['Events'].text()}\n")

        self.log(f"Wrote {run_mac_path}")

    def run_panda(self):
        # Recoil-hit logging is opt-in and off by default (it adds real
        # per-step file-write overhead most runs don't need -- see
        # EventAction.hh), so leaving it unchecked is usually deliberate,
        # not a mistake. Warn rather than block: PANDA itself and
        # PANDA_Analyze.py's charge-spectrum analysis are unaffected
        # either way, only PANDAEX's per-species/LET breakdown needs it.
        if not self.log_recoil_hits_checkbox.isChecked():
            reply = QMessageBox.warning(
                self,
                "Recoil-hit logging disabled",
                "\"Log Recoil Hits\" is unchecked -- this run will not "
                "produce recoil_hits.csv, so PANDAEX's per-species/LET "
                "breakdown won't be available afterward (Run PANDA "
                "itself and PANDA_Analyze.py's charge-spectrum analysis "
                "are unaffected).\n\n"
                "Run anyway?",
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.Yes,
            )

            if reply != QMessageBox.Yes:
                self.log("Run PANDA cancelled -- Log Recoil Hits is unchecked.")
                return

        self.write_run_mac()
        self.log("Running PANDA...")

        # build/PANDA, not a root-level copy: pointing directly at the
        # binary CMake actually produces means a rebuild is immediately
        # picked up, with no separate manual "cp build/PANDA ." step
        # that can silently go stale (see commit history).
        panda_path = os.path.join(self.build_dir, "PANDA")

        process = subprocess.Popen(
            [panda_path, "run.mac"],
            cwd=self.project_dir,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True
        )

        for line in process.stdout:
            self.log(line.strip())

        self.log("PANDA run complete.")

    def run_analysis(self):
        self.log("Running PANDA_Analyze.py...")

        process = subprocess.Popen(
            [
                "python3",
                os.path.join(self.project_dir, "PANDA_Analyze.py"),
            ],
            cwd=self.project_dir,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True
        )

        for line in process.stdout:
            self.log(line.strip())

        self.log("Analysis complete.")

    def run_pandaex(self):
        # Ground truth for "was the run that produced the current
        # results logged with recoil hits" is the output file itself,
        # not the Log Recoil Hits checkbox -- the checkbox reflects
        # whatever the form is currently set to for the NEXT run, which
        # can easily be out of sync with what actually produced the
        # results sitting in Results/Current right now (checkbox
        # toggled after the last run, a preset loaded without
        # re-running, etc).
        recoil_hits_path = os.path.join(self.results_dir, "recoil_hits.csv")

        if not os.path.exists(recoil_hits_path):
            reply = QMessageBox.warning(
                self,
                "No recoil-hit data found",
                "Results/Current/recoil_hits.csv is missing -- the run "
                "that produced the current results was likely run with "
                "\"Log Recoil Hits\" unchecked (or hasn't been run yet).\n\n"
                "PANDAEX_Analyze.py's per-species/LET breakdown needs "
                "that file; without it, only the raw deposited-energy "
                "component spectrum will be produced.\n\n"
                "Run PANDAEX anyway?",
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.No,
            )

            if reply != QMessageBox.Yes:
                self.log("Run PANDAEX cancelled -- no recoil_hits.csv found.")
                return

        self.log("Running PANDAEX_Analyze.py...")

        # PANDAEX_Analyze.py ends with a blocking plt.show() (every plot
        # is already saved to disk via savefig() before that call) --
        # left in for standalone command-line use, where seeing the
        # figures pop up is convenient. Run here with the non-interactive
        # Agg backend instead, or that call blocks this subprocess (and
        # therefore this synchronous stdout-reading loop, and therefore
        # the whole GUI) until someone manually closes windows nobody
        # asked for.
        env = os.environ.copy()
        env["MPLBACKEND"] = "Agg"

        process = subprocess.Popen(
            [
                "python3",
                os.path.join(self.project_dir, "PANDAEX_Analyze.py"),
            ],
            cwd=self.project_dir,
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True
        )

        for line in process.stdout:
            self.log(line.strip())

        self.log("PANDAEX analysis complete.")

    def save_preset(self):
        filename, _ = QFileDialog.getSaveFileName(
            self, "Save Preset", "", "JSON Files (*.json)"
        )

        if filename:
            data = {
                "particle": self.particle_dropdown.currentText(),
                "sensitive_material": self.sensitive_material_dropdown.currentText(),
                "dead_material": self.dead_material_dropdown.currentText(),
                **{k: v.text() for k, v in self.fields.items()},
                "visualization": self.vis_checkbox.isChecked(),
                "collection_model": self.collection_checkbox.isChecked(),
                "link_xy": self.link_xy_checkbox.isChecked(),
                "log_recoil_hits": self.log_recoil_hits_checkbox.isChecked()
            }

            if not filename.endswith(".json"):
                filename += ".json"

            with open(filename, "w") as f:
                json.dump(data, f, indent=4)

            self.log(f"Saved preset: {filename}")

    def load_preset(self):
        filename, _ = QFileDialog.getOpenFileName(
            self, "Load Preset", "", "JSON Files (*.json)"
        )

        if filename:
            with open(filename, "r") as f:
                data = json.load(f)

            self.particle_dropdown.setCurrentText(
                data.get("particle", "proton")
            )

            self.sensitive_material_dropdown.setCurrentText(
                data.get("sensitive_material", "Si")
            )

            self.dead_material_dropdown.setCurrentText(
                data.get("dead_material", "SiO2")
            )

            # Set link state first: toggling it may sync Dead/Beam XY
            # from Sensitive XY, but the field loop below runs after
            # and restores the actual saved values, so this ordering
            # can't clobber an explicitly-unlinked preset's values.
            self.link_xy_checkbox.setChecked(
                data.get("link_xy", True)
            )

            for k in self.fields:
                if k in data:
                    self.fields[k].setText(data[k])

            self.vis_checkbox.setChecked(
                data.get("visualization", False)
            )

            self.collection_checkbox.setChecked(
                data.get("collection_model", True)
            )

            self.log_recoil_hits_checkbox.setChecked(
                data.get("log_recoil_hits", False)
            )

            self.log(f"Loaded preset: {filename}")

    def open_results_folder(self):
        subprocess.Popen(
            ["xdg-open", self.results_dir],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            start_new_session=True
        )

    def open_file(self, filename):
        filepath = os.path.join(self.results_dir, filename)

        if not os.path.exists(filepath):
            self.log(f"File not found: {filename}")
            return

        dialog = QDialog(self)
        dialog.setWindowTitle(filename)
        dialog.resize(1000, 800)

        layout = QVBoxLayout()

        image_label = QLabel()
        image_label.setAlignment(Qt.AlignCenter)

        pixmap = QPixmap(filepath)
        pixmap = pixmap.scaled(
            950,
            750,
            Qt.KeepAspectRatio,
            Qt.SmoothTransformation
        )

        image_label.setPixmap(pixmap)

        layout.addWidget(image_label)
        dialog.setLayout(layout)

        dialog.exec_()


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = PandaGUI()
    window.show()
    sys.exit(app.exec_())
