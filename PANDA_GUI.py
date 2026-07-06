import sys
import os
import subprocess
import json

from PyQt5.QtWidgets import (
    QApplication, QWidget, QLabel, QLineEdit,
    QPushButton, QTextEdit, QVBoxLayout, QHBoxLayout,
    QGridLayout, QFileDialog, QCheckBox,
    QGroupBox, QComboBox, QDialog
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

        for i, (label, default) in enumerate(labels, start=1):
            param_layout.addWidget(QLabel(label + ":"), i, 0)

            field = QLineEdit(default)
            field.setFixedWidth(220)
            field.setFixedHeight(34)

            self.fields[label] = field
            param_layout.addWidget(field, i, 1)

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

        self.fields['Sensitive XY (um)'].textChanged.connect(self.on_sensitive_xy_changed)
        self.on_link_xy_toggled(self.link_xy_checkbox.isChecked())

        param_group.setLayout(param_layout)

        # ==========================
        # Controls
        # ==========================
        control_group = QGroupBox("Controls")
        control_layout = QGridLayout()

        control_layout.addWidget(QLabel("Plot Metric:"), 0, 0, 1, 2)

        self.metric_dropdown = QComboBox()
        self.metric_dropdown.addItems(["Collected", "Deposited"])
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

        self.plot_spec_button = QPushButton("▥ Plot Spectrum")
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
        self.style_button(self.open_results_button, "#f5f5f5", "#999999")
        self.open_results_button.clicked.connect(self.open_results_folder)
        control_layout.addWidget(self.open_results_button, 5, 0, 1, 2)

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
            f.write(f"/sim/particle {self.particle_dropdown.currentText()}\n")
            f.write(f"/sim/energy {self.fields['Energy (MeV)'].text()} MeV\n")
            f.write(f"/sim/sensitiveXY {self.fields['Sensitive XY (um)'].text()} um\n")
            f.write(f"/sim/deadXY {self.fields['Dead XY (um)'].text()} um\n")
            f.write(f"/sim/sensitiveThickness {self.fields['Sensitive Thickness (um)'].text()} um\n")
            f.write(f"/sim/deadThickness {self.fields['Dead Thickness (um)'].text()} um\n")
            f.write(f"/sim/beamXY {self.fields['Beam XY (um)'].text()} um\n")
            f.write(f"/sim/criticalCharge {self.fields['Critical Charge (fC)'].text()}\n")
            f.write(f"/sim/biasCrossSectionFactor {self.fields['Bias Cross Section Factor'].text()}\n")

            if self.collection_checkbox.isChecked():
                f.write("/sim/useCollectionModel true\n")
            else:
                f.write("/sim/useCollectionModel false\n")

            f.write(f"/sim/carrierLifetime {self.default_carrier_lifetime}\n")
            f.write(f"/sim/electricField {self.default_electric_field}\n")

            f.write("/sim/verbose false\n\n")
            f.write("/run/initialize\n")
            f.write("/run/reinitializeGeometry\n\n")

            if self.vis_checkbox.isChecked():
                f.write("/vis/open OGLSX\n")
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
                f.write("/vis/viewer/flush\n\n")

            f.write(f"/run/beamOn {self.fields['Events'].text()}\n")

        self.log(f"Wrote {run_mac_path}")

    def run_panda(self):
        self.write_run_mac()
        self.log("Running PANDA...")

        panda_path = os.path.join(self.project_dir, "PANDA")

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

    def save_preset(self):
        filename, _ = QFileDialog.getSaveFileName(
            self, "Save Preset", "", "JSON Files (*.json)"
        )

        if filename:
            data = {
                "particle": self.particle_dropdown.currentText(),
                **{k: v.text() for k, v in self.fields.items()},
                "visualization": self.vis_checkbox.isChecked(),
                "collection_model": self.collection_checkbox.isChecked(),
                "link_xy": self.link_xy_checkbox.isChecked()
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
