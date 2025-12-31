"""
Coronary Artery Calcium Score - Guided Workflow Plugin

A 3D Slicer extension for quantification of coronary artery calcification using
the Agatston scoring method. Features per-vessel analysis (LAD, LCx, RCA, LM),
MESA/ACC risk classification, and MESA percentile calculation.

Author: Vittorio Censullo
Institution: AITeRTC (Associazione Italiana Tecnici di Radiologia Esperti in TC)
Year: 2025
Version: 1.0
License: MIT

Based on ACC/AHA guidelines and MESA Study:

STANDARD AGATSTON SCORING (Non-Contrast CT):
- Fixed threshold: 130 HU (standard for non-contrast cardiac CT)
- Per-vessel analysis: LAD, LCx, RCA, Left Main
- MESA/ACC Risk Categories: 0, 1-10, 11-100, 101-400, >400 AU
- MESA Percentiles: Age/Sex/Ethnicity-based
- 2D slice-by-slice labeling (matches commercial software)
- Minimum lesion criterion: 1 mm² per slice (Agatston standard)
- Density-based weighting (factors 1-4)
- 3D visualization and PDF reporting

Key Features:
- Per-vessel calcium scoring (LAD, LCx, RCA, LM)
- Color-coded territory visualization
- MESA percentile calculation
- Comprehensive PDF reporting with per-vessel breakdown

Repository: https://github.com/vcensullo/CoronaryCaScore
"""

import os
import json
import unittest
import vtk
import qt
import ctk
import slicer
from slicer.ScriptedLoadableModule import *
from slicer.util import VTKObservationMixin
import numpy as np
from datetime import datetime

#
# CoronaryCaScore
#

class CoronaryCaScore(ScriptedLoadableModule):
    """Coronary Artery Calcium Score Module"""

    def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        self.parent.title = "Coronary Artery Calcium Score"
        self.parent.categories = ["Cardiac", "Quantification"]
        self.parent.dependencies = []
        self.parent.contributors = ["Vittorio Censullo (AITeRTC - 2025)"]
        self.parent.helpText = """
        <h3>Coronary Artery Calcium Score - Per-Vessel Analysis (v1.0)</h3>

        <p>Quantification of coronary artery calcification using the Agatston scoring method with per-vessel analysis.</p>

        <h4>Features:</h4>
        <ul>
        <li>Guided tab-based workflow for easy use</li>
        <li><b>Per-vessel analysis:</b> LAD, LCx, RCA, Left Main</li>
        <li>Color-coded territory visualization</li>
        <li>MESA/ACC risk classification</li>
        <li>MESA percentile calculation (age/sex/ethnicity)</li>
        <li>ROI-based or manual paint segmentation</li>
        <li>3D visualization with vessel-based coloring</li>
        <li>Comprehensive PDF reporting</li>
        </ul>

        <h4>MESA/ACC Risk Categories:</h4>
        <ul>
        <li><b>Zero:</b> 0 AU - Very Low Risk</li>
        <li><b>Minimal:</b> 1-10 AU - Low Risk</li>
        <li><b>Mild:</b> 11-100 AU - Low-Moderate Risk</li>
        <li><b>Moderate:</b> 101-400 AU - Moderate-High Risk</li>
        <li><b>Severe:</b> >400 AU - High Risk</li>
        </ul>

        <h4>Author Information:</h4>
        <ul>
        <li><b>Author:</b> Vittorio Censullo</li>
        <li><b>Institution:</b> AITeRTC (Associazione Italiana Tecnici di Radiologia Esperti in TC)</li>
        <li><b>Year:</b> 2025</li>
        <li><b>Version:</b> 1.0</li>
        <li><b>License:</b> MIT</li>
        <li><b>Repository:</b> <a href="https://github.com/vcensullo/CoronaryCaScore">GitHub</a></li>
        </ul>

        <h4>References:</h4>
        <p>Based on ACC/AHA guidelines and MESA Study:</p>
        <ul>
        <li>Agatston et al. J Am Coll Cardiol 1990</li>
        <li>MESA Study - McClelland et al. 2006</li>
        <li>ACC/AHA Guidelines 2019</li>
        </ul>
        """
        self.parent.acknowledgementText = """
        <p>Developed by <b>Vittorio Censullo</b> at <b>AITeRTC</b> (Associazione Italiana Tecnici di Radiologia Esperti in TC) - 2025</p>

        <p>This plugin provides coronary artery calcium scoring according to
        ACC/AHA guidelines with MESA percentile calculation.</p>

        <p>For issues, suggestions, or contributions, please visit the
        <a href="https://github.com/vcensullo/CoronaryCaScore">GitHub repository</a>.</p>
        """

#
# CoronaryCaScoreWidget
#

class CoronaryCaScoreWidget(ScriptedLoadableModuleWidget, VTKObservationMixin):
    """Main widget with tab-based guided workflow"""

    # Territory definitions with colors
    TERRITORIES = {
        'LAD': {'name': 'Left Anterior Descending', 'color': (1.0, 0.0, 0.0), 'colorHex': '#FF0000'},
        'LCx': {'name': 'Left Circumflex', 'color': (0.0, 0.0, 1.0), 'colorHex': '#0000FF'},
        'RCA': {'name': 'Right Coronary Artery', 'color': (0.0, 0.8, 0.0), 'colorHex': '#00CC00'},
        'LM': {'name': 'Left Main', 'color': (1.0, 0.8, 0.0), 'colorHex': '#FFCC00'}
    }

    ETHNICITIES = ['Caucasian', 'African-American', 'Hispanic', 'Chinese']

    def __init__(self, parent=None):
        ScriptedLoadableModuleWidget.__init__(self, parent)
        VTKObservationMixin.__init__(self)
        self.logic = None
        self._parameterNode = None
        self._updatingGUIFromParameterNode = False

        # Current state
        self.currentTerritory = 'LAD'
        self.segmentationsByTerritory = {t: None for t in self.TERRITORIES}
        self.segmentEditorWidget = None
        self.segmentEditorNode = None

        # Patient info
        self.patientInfo = {
            'name': '',
            'id': '',
            'sex': 'M',
            'age': '',
            'ethnicity': 'Caucasian',
            'date': datetime.now().strftime("%Y-%m-%d")
        }

        # Results storage
        self.currentResults = None

    def setup(self):
        ScriptedLoadableModuleWidget.setup(self)

        # Load settings
        self.settings = qt.QSettings("AITeRTC", "CoronaryCaScore")

        # Create logic
        self.logic = CoronaryCaScoreLogic()

        # Create main layout
        self.layout.setContentsMargins(6, 6, 6, 6)
        self.layout.setSpacing(8)

        # Apply global stylesheet for cleaner look
        self.applyGlobalStylesheet()

        # Create header
        self.setupHeader()

        # Create patient info section (collapsible)
        self.setupPatientInfoSection()

        # Create tab widget
        self.tabWidget = qt.QTabWidget()
        self.tabWidget.setDocumentMode(True)
        self.layout.addWidget(self.tabWidget)

        # Setup tabs
        self.setupSettingsTab()      # Tab 0
        self.setupSetupViewTab()     # Tab 1
        self.setupSegmentationTab()  # Tab 2
        self.setupResultsTab()       # Tab 3
        self.setupReportTab()        # Tab 4

        # Initially disable tabs 2-4
        self.tabWidget.setTabEnabled(2, False)
        self.tabWidget.setTabEnabled(3, False)
        self.tabWidget.setTabEnabled(4, False)

        # Add stretch
        self.layout.addStretch(1)

        # Load saved settings
        self.loadSettings()

    def applyGlobalStylesheet(self):
        """Apply a clean, modern stylesheet to the widget"""
        stylesheet = """
            /* General styling */
            QWidget {
                font-family: 'Segoe UI', Arial, sans-serif;
            }

            /* Collapsible buttons */
            ctkCollapsibleButton {
                background-color: #f5f5f5;
                border: 1px solid #ddd;
                border-radius: 4px;
                padding: 8px;
                font-weight: bold;
                color: #333;
            }
            ctkCollapsibleButton:hover {
                background-color: #e8e8e8;
            }

            /* Primary action buttons */
            QPushButton#primaryButton {
                background-color: #2196F3;
                color: white;
                border: none;
                border-radius: 4px;
                padding: 10px 20px;
                font-weight: bold;
                font-size: 13px;
            }
            QPushButton#primaryButton:hover {
                background-color: #1976D2;
            }
            QPushButton#primaryButton:pressed {
                background-color: #0D47A1;
            }

            /* Secondary buttons */
            QPushButton#secondaryButton {
                background-color: #fff;
                color: #333;
                border: 1px solid #ccc;
                border-radius: 4px;
                padding: 8px 16px;
            }
            QPushButton#secondaryButton:hover {
                background-color: #f0f0f0;
                border-color: #999;
            }

            /* Danger buttons */
            QPushButton#dangerButton {
                background-color: #fff;
                color: #d32f2f;
                border: 1px solid #d32f2f;
                border-radius: 4px;
                padding: 8px 16px;
            }
            QPushButton#dangerButton:hover {
                background-color: #ffebee;
            }

            /* Tables */
            QTableWidget {
                border: 1px solid #ddd;
                border-radius: 4px;
                gridline-color: #eee;
                background-color: white;
            }
            QTableWidget::item {
                padding: 6px;
            }
            QHeaderView::section {
                background-color: #f5f5f5;
                border: none;
                border-bottom: 2px solid #2196F3;
                padding: 8px;
                font-weight: bold;
            }

            /* Input fields */
            QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox {
                border: 1px solid #ccc;
                border-radius: 4px;
                padding: 6px;
                background-color: white;
            }
            QLineEdit:focus, QSpinBox:focus, QDoubleSpinBox:focus, QComboBox:focus {
                border-color: #2196F3;
            }

            /* Checkboxes */
            QCheckBox {
                spacing: 8px;
            }

            /* Labels */
            QLabel#sectionTitle {
                font-size: 14px;
                font-weight: bold;
                color: #1976D2;
                padding: 4px 0;
            }
            QLabel#statusSuccess {
                color: #4CAF50;
                font-weight: bold;
            }
            QLabel#statusError {
                color: #d32f2f;
                font-weight: bold;
            }
            QLabel#statusInfo {
                color: #666;
                font-style: italic;
            }
        """
        self.parent.setStyleSheet(stylesheet)

    def setupHeader(self):
        """Setup plugin header with branding"""
        headerFrame = qt.QFrame()
        headerFrame.setStyleSheet("""
            QFrame {
                background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
                    stop:0 #1565C0, stop:1 #1976D2);
                border-radius: 6px;
                margin-bottom: 8px;
            }
        """)
        headerLayout = qt.QHBoxLayout(headerFrame)
        headerLayout.setContentsMargins(16, 12, 16, 12)

        # Title and version
        titleLabel = qt.QLabel("<span style='color: white; font-size: 18px; font-weight: bold;'>CoronaryCaScore</span>")
        headerLayout.addWidget(titleLabel)

        headerLayout.addStretch()

        # Version badge
        versionLabel = qt.QLabel("<span style='color: rgba(255,255,255,0.9); font-size: 12px;'>v1.0</span>")
        headerLayout.addWidget(versionLabel)

        self.layout.addWidget(headerFrame)

    def setupPatientInfoSection(self):
        """Setup collapsible patient information section"""
        patientCollapsible = ctk.ctkCollapsibleButton()
        patientCollapsible.text = "Patient Information"
        patientCollapsible.collapsed = False  # Expanded by default for visibility
        self.layout.addWidget(patientCollapsible)

        patientLayout = qt.QFormLayout(patientCollapsible)

        # Name
        self.patientNameEdit = qt.QLineEdit()
        self.patientNameEdit.setPlaceholderText("Patient Name (optional)")
        patientLayout.addRow("Name:", self.patientNameEdit)

        # ID
        self.patientIdEdit = qt.QLineEdit()
        self.patientIdEdit.setPlaceholderText("Patient ID (optional)")
        patientLayout.addRow("ID:", self.patientIdEdit)

        # Sex (mandatory)
        sexLayout = qt.QHBoxLayout()
        self.sexButtonGroup = qt.QButtonGroup()
        self.maleRadio = qt.QRadioButton("Male")
        self.femaleRadio = qt.QRadioButton("Female")
        self.maleRadio.setChecked(True)
        self.sexButtonGroup.addButton(self.maleRadio, 0)
        self.sexButtonGroup.addButton(self.femaleRadio, 1)
        sexLayout.addWidget(self.maleRadio)
        sexLayout.addWidget(self.femaleRadio)
        sexLayout.addStretch()
        patientLayout.addRow("Sex*:", sexLayout)

        # Age (required for MESA percentile)
        self.patientAgeSpinBox = qt.QSpinBox()
        self.patientAgeSpinBox.setRange(0, 120)
        self.patientAgeSpinBox.setValue(0)
        self.patientAgeSpinBox.setSpecialValueText("Not specified")
        patientLayout.addRow("Age*:", self.patientAgeSpinBox)

        # Ethnicity (for MESA percentiles)
        self.ethnicityComboBox = qt.QComboBox()
        self.ethnicityComboBox.addItems(self.ETHNICITIES)
        patientLayout.addRow("Ethnicity*:", self.ethnicityComboBox)

        # Mandatory note
        mandatoryLabel = qt.QLabel("<i>* Sex, Age (45+), and Ethnicity required for MESA percentile</i>")
        mandatoryLabel.setStyleSheet("color: gray;")
        patientLayout.addRow("", mandatoryLabel)

    def setupSettingsTab(self):
        """Setup Tab 0: Settings"""
        settingsTab = qt.QWidget()
        settingsLayout = qt.QVBoxLayout(settingsTab)

        # Dependencies section
        depsCollapsible = ctk.ctkCollapsibleButton()
        depsCollapsible.text = "Dependencies"
        settingsLayout.addWidget(depsCollapsible)
        depsLayout = qt.QVBoxLayout(depsCollapsible)

        depsLabel = qt.QLabel("Required libraries: reportlab, matplotlib")
        depsLayout.addWidget(depsLabel)

        self.installDepsButton = qt.QPushButton("Install/Update Dependencies")
        self.installDepsButton.setObjectName("secondaryButton")
        self.installDepsButton.setMinimumHeight(36)
        self.installDepsButton.clicked.connect(self.onInstallDependencies)
        depsLayout.addWidget(self.installDepsButton)

        self.depsStatusLabel = qt.QLabel("")
        depsLayout.addWidget(self.depsStatusLabel)

        # Threshold section
        thresholdCollapsible = ctk.ctkCollapsibleButton()
        thresholdCollapsible.text = "Agatston Threshold"
        settingsLayout.addWidget(thresholdCollapsible)
        thresholdLayout = qt.QFormLayout(thresholdCollapsible)

        self.thresholdSpinBox = qt.QSpinBox()
        self.thresholdSpinBox.setRange(100, 200)
        self.thresholdSpinBox.setValue(130)
        self.thresholdSpinBox.setSuffix(" HU")
        thresholdLayout.addRow("Threshold:", self.thresholdSpinBox)

        thresholdNote = qt.QLabel("<i>Standard Agatston threshold: 130 HU</i>")
        thresholdNote.setStyleSheet("color: gray;")
        thresholdLayout.addRow("", thresholdNote)

        # Mass calibration section
        massCollapsible = ctk.ctkCollapsibleButton()
        massCollapsible.text = "Mass Calibration"
        settingsLayout.addWidget(massCollapsible)
        massLayout = qt.QFormLayout(massCollapsible)

        self.calibrationFactorSpinBox = qt.QDoubleSpinBox()
        self.calibrationFactorSpinBox.setRange(0.50, 1.50)
        self.calibrationFactorSpinBox.setValue(0.81)
        self.calibrationFactorSpinBox.setSingleStep(0.01)
        self.calibrationFactorSpinBox.setDecimals(2)
        massLayout.addRow("Calibration Factor:", self.calibrationFactorSpinBox)

        calibrationNote = qt.QLabel(
            "<i>Scanner-specific factor from phantom calibration.<br>"
            "Default: 0.81 (for 130 HU = 114.5 mg/cm³ CaHA)</i>"
        )
        calibrationNote.setStyleSheet("color: gray;")
        massLayout.addRow("", calibrationNote)

        # Company branding section
        brandingCollapsible = ctk.ctkCollapsibleButton()
        brandingCollapsible.text = "Report Branding"
        settingsLayout.addWidget(brandingCollapsible)
        brandingLayout = qt.QFormLayout(brandingCollapsible)

        # Logo
        logoLayout = qt.QHBoxLayout()
        self.logoPathEdit = qt.QLineEdit()
        self.logoPathEdit.setPlaceholderText("Path to company logo")
        logoLayout.addWidget(self.logoPathEdit)
        self.browseLogoButton = qt.QPushButton("Browse...")
        self.browseLogoButton.setMaximumWidth(80)
        self.browseLogoButton.clicked.connect(self.onBrowseLogo)
        logoLayout.addWidget(self.browseLogoButton)
        brandingLayout.addRow("Logo:", logoLayout)

        # Company description
        self.companyDescEdit = qt.QLineEdit()
        self.companyDescEdit.setPlaceholderText("Company description for PDF header")
        brandingLayout.addRow("Description:", self.companyDescEdit)

        # Save settings button
        self.saveSettingsButton = qt.QPushButton("Save Settings")
        self.saveSettingsButton.setObjectName("secondaryButton")
        self.saveSettingsButton.setMinimumHeight(36)
        self.saveSettingsButton.clicked.connect(self.saveSettings)
        settingsLayout.addWidget(self.saveSettingsButton)

        # Plugin info
        infoCollapsible = ctk.ctkCollapsibleButton()
        infoCollapsible.text = "Plugin Information"
        infoCollapsible.collapsed = True
        settingsLayout.addWidget(infoCollapsible)
        infoLayout = qt.QVBoxLayout(infoCollapsible)

        infoText = qt.QLabel(
            "<b>Coronary Artery Calcium Score</b><br>"
            "Version: 1.0<br>"
            "Author: Vittorio Censullo<br>"
            "Institution: AITeRTC<br>"
            "Year: 2025<br><br>"
            "Purpose: Quantification of coronary artery calcification using Standard Agatston scoring<br>"
            "with per-vessel analysis and MESA percentile calculation."
        )
        infoText.setWordWrap(True)
        infoLayout.addWidget(infoText)

        settingsLayout.addStretch(1)
        self.tabWidget.addTab(settingsTab, "0. Settings")

    def setupSetupViewTab(self):
        """Setup Tab 1: Setup View"""
        setupTab = qt.QWidget()
        setupLayout = qt.QVBoxLayout(setupTab)

        # Volume selection
        volumeCollapsible = ctk.ctkCollapsibleButton()
        volumeCollapsible.text = "Volume Selection"
        setupLayout.addWidget(volumeCollapsible)
        volumeLayout = qt.QFormLayout(volumeCollapsible)

        self.volumeSelector = slicer.qMRMLNodeComboBox()
        self.volumeSelector.nodeTypes = ["vtkMRMLScalarVolumeNode"]
        self.volumeSelector.selectNodeUponCreation = True
        self.volumeSelector.addEnabled = False
        self.volumeSelector.removeEnabled = False
        self.volumeSelector.noneEnabled = True
        self.volumeSelector.showHidden = False
        self.volumeSelector.setMRMLScene(slicer.mrmlScene)
        self.volumeSelector.setToolTip("Select the CT volume for calcium scoring")
        self.volumeSelector.currentNodeChanged.connect(self.onVolumeSelected)
        volumeLayout.addRow("CT Volume:", self.volumeSelector)

        # View setup
        viewCollapsible = ctk.ctkCollapsibleButton()
        viewCollapsible.text = "View Setup"
        setupLayout.addWidget(viewCollapsible)
        viewLayout = qt.QVBoxLayout(viewCollapsible)

        self.applyLayoutButton = qt.QPushButton("Apply Optimal Layout")
        self.applyLayoutButton.setObjectName("secondaryButton")
        self.applyLayoutButton.setMinimumHeight(36)
        self.applyLayoutButton.setToolTip("Set up Four-Up view with cardiac window/level")
        self.applyLayoutButton.clicked.connect(self.onApplyLayout)
        viewLayout.addWidget(self.applyLayoutButton)

        self.enableIntersectionsCheck = qt.QCheckBox("Enable Slice Intersections")
        self.enableIntersectionsCheck.setChecked(True)
        viewLayout.addWidget(self.enableIntersectionsCheck)

        # Status
        self.setupStatusLabel = qt.QLabel("")
        setupLayout.addWidget(self.setupStatusLabel)

        setupLayout.addStretch(1)
        self.tabWidget.addTab(setupTab, "1. Setup View")

    def setupSegmentationTab(self):
        """Setup Tab 2: Calcium Segmentation with territory selection"""
        segTab = qt.QWidget()
        segLayout = qt.QVBoxLayout(segTab)

        # Territory selection
        territoryCollapsible = ctk.ctkCollapsibleButton()
        territoryCollapsible.text = "Coronary Territory Selection"
        segLayout.addWidget(territoryCollapsible)
        territoryLayout = qt.QVBoxLayout(territoryCollapsible)

        # Territory buttons with color indicators
        territoryButtonLayout = qt.QHBoxLayout()
        territoryButtonLayout.setSpacing(10)
        self.territoryButtonGroup = qt.QButtonGroup()
        self.territoryButtons = {}

        for i, (abbrev, info) in enumerate(self.TERRITORIES.items()):
            btn = qt.QPushButton(abbrev)
            btn.setCheckable(True)
            btn.setMinimumHeight(50)
            btn.setMinimumWidth(70)
            btn.setStyleSheet(f"""
                QPushButton {{
                    background-color: white;
                    border: 3px solid {info['colorHex']};
                    border-radius: 8px;
                    font-weight: bold;
                    font-size: 15px;
                    color: #333;
                }}
                QPushButton:hover {{
                    background-color: {info['colorHex']}22;
                }}
                QPushButton:checked {{
                    background-color: {info['colorHex']};
                    color: white;
                    border-color: {info['colorHex']};
                }}
            """)
            btn.setToolTip(info['name'])
            self.territoryButtonGroup.addButton(btn, i)
            self.territoryButtons[abbrev] = btn
            territoryButtonLayout.addWidget(btn)

        self.territoryButtons['LAD'].setChecked(True)
        self.territoryButtonGroup.buttonClicked.connect(self.onTerritoryChanged)
        territoryLayout.addLayout(territoryButtonLayout)

        # Current territory info
        self.currentTerritoryLabel = qt.QLabel("<b>Current Territory:</b> LAD - Left Anterior Descending")
        territoryLayout.addWidget(self.currentTerritoryLabel)

        # Segmentation tools
        toolsCollapsible = ctk.ctkCollapsibleButton()
        toolsCollapsible.text = "Segmentation Tools"
        segLayout.addWidget(toolsCollapsible)
        toolsLayout = qt.QVBoxLayout(toolsCollapsible)

        # Tool buttons - styled for better visual feedback
        toolButtonLayout = qt.QHBoxLayout()
        toolButtonLayout.setSpacing(8)

        toolButtonStyle = """
            QPushButton {
                background-color: #fff;
                border: 2px solid #ccc;
                border-radius: 6px;
                padding: 8px 12px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #f5f5f5;
                border-color: #999;
            }
            QPushButton:checked {
                background-color: #E3F2FD;
                border-color: #2196F3;
                color: #1565C0;
            }
        """

        # Click & Grow method (default, primary tool)
        self.clickGrowButton = qt.QPushButton(" Click && Grow")
        self.clickGrowButton.setCheckable(True)
        self.clickGrowButton.setMinimumHeight(44)
        self.clickGrowButton.setIcon(qt.QIcon(os.path.join(os.path.dirname(__file__), 'Resources/Icons/click_grow.png')))
        self.clickGrowButton.setStyleSheet(toolButtonStyle)
        self.clickGrowButton.clicked.connect(self.onClickGrowToggle)
        toolButtonLayout.addWidget(self.clickGrowButton)

        # Brush method (for refinements)
        self.brushMethodButton = qt.QPushButton(" Brush")
        self.brushMethodButton.setCheckable(True)
        self.brushMethodButton.setMinimumHeight(44)
        self.brushMethodButton.setIcon(qt.QIcon(os.path.join(os.path.dirname(__file__), 'Resources/Icons/brush.png')))
        self.brushMethodButton.setStyleSheet(toolButtonStyle)
        self.brushMethodButton.clicked.connect(self.onBrushToggle)
        toolButtonLayout.addWidget(self.brushMethodButton)

        # Erase button
        self.eraseButton = qt.QPushButton(" Erase")
        self.eraseButton.setCheckable(True)
        self.eraseButton.setMinimumHeight(44)
        eraseButtonStyle = """
            QPushButton {
                background-color: #fff;
                border: 2px solid #ccc;
                border-radius: 6px;
                padding: 8px 12px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #fff5f5;
                border-color: #d32f2f;
            }
            QPushButton:checked {
                background-color: #FFEBEE;
                border-color: #d32f2f;
                color: #c62828;
            }
        """
        self.eraseButton.setStyleSheet(eraseButtonStyle)
        self.eraseButton.clicked.connect(self.onEraseToggle)
        toolButtonLayout.addWidget(self.eraseButton)

        toolsLayout.addLayout(toolButtonLayout)

        # Instructions label
        self.toolInstructionsLabel = qt.QLabel("<i>Select a tool to start segmenting calcium deposits</i>")
        self.toolInstructionsLabel.setStyleSheet("color: gray; padding: 5px;")
        self.toolInstructionsLabel.setWordWrap(True)
        toolsLayout.addWidget(self.toolInstructionsLabel)

        # Track active tool
        self.activeToolName = None

        # Territory summary table
        summaryCollapsible = ctk.ctkCollapsibleButton()
        summaryCollapsible.text = "Territory Summary"
        segLayout.addWidget(summaryCollapsible)
        summaryLayout = qt.QVBoxLayout(summaryCollapsible)

        self.territorySummaryTable = qt.QTableWidget()
        self.territorySummaryTable.setColumnCount(3)
        self.territorySummaryTable.setHorizontalHeaderLabels(["Territory", "Lesions", "Status"])
        self.territorySummaryTable.setRowCount(4)
        self.territorySummaryTable.horizontalHeader().setStretchLastSection(True)
        self.territorySummaryTable.setMaximumHeight(150)
        summaryLayout.addWidget(self.territorySummaryTable)

        self.updateTerritorySummary()

        # Clear buttons
        clearLayout = qt.QHBoxLayout()
        self.clearTerritoryButton = qt.QPushButton("Clear Current Territory")
        self.clearTerritoryButton.setObjectName("dangerButton")
        self.clearTerritoryButton.setMinimumHeight(32)
        self.clearTerritoryButton.clicked.connect(self.onClearTerritory)
        clearLayout.addWidget(self.clearTerritoryButton)

        self.clearAllButton = qt.QPushButton("Clear All")
        self.clearAllButton.setObjectName("dangerButton")
        self.clearAllButton.setMinimumHeight(32)
        self.clearAllButton.clicked.connect(self.onClearAllTerritories)
        clearLayout.addWidget(self.clearAllButton)
        segLayout.addLayout(clearLayout)

        # Status
        self.segmentationStatusLabel = qt.QLabel("")
        segLayout.addWidget(self.segmentationStatusLabel)

        segLayout.addStretch(1)
        self.tabWidget.addTab(segTab, "2. Segmentation")

    def setupResultsTab(self):
        """Setup Tab 3: Results & Analysis"""
        resultsTab = qt.QWidget()
        resultsLayout = qt.QVBoxLayout(resultsTab)

        # Calculate button
        self.calculateButton = qt.QPushButton("  Calculate All Scores")
        self.calculateButton.setObjectName("primaryButton")
        self.calculateButton.setMinimumHeight(44)
        self.calculateButton.clicked.connect(self.onCalculate)
        resultsLayout.addWidget(self.calculateButton)

        # Results table by territory
        resultsCollapsible = ctk.ctkCollapsibleButton()
        resultsCollapsible.text = "Agatston Score by Territory"
        resultsLayout.addWidget(resultsCollapsible)
        tableLayout = qt.QVBoxLayout(resultsCollapsible)

        self.resultsTable = qt.QTableWidget()
        self.resultsTable.setColumnCount(5)
        self.resultsTable.setHorizontalHeaderLabels(["Territory", "Score (AU)", "Volume (mm³)", "Mass (mg)", "Lesions"])
        self.resultsTable.setRowCount(5)  # 4 territories + total
        self.resultsTable.horizontalHeader().setStretchLastSection(True)
        tableLayout.addWidget(self.resultsTable)

        # Risk classification
        riskCollapsible = ctk.ctkCollapsibleButton()
        riskCollapsible.text = "Risk Classification"
        resultsLayout.addWidget(riskCollapsible)
        riskLayout = qt.QVBoxLayout(riskCollapsible)

        # MESA/ACC Category
        self.mesaCategoryLabel = qt.QLabel("<b>MESA/ACC Category:</b> --")
        self.mesaCategoryLabel.setStyleSheet("font-size: 14px;")
        riskLayout.addWidget(self.mesaCategoryLabel)

        self.riskDescriptionLabel = qt.QLabel("")
        riskLayout.addWidget(self.riskDescriptionLabel)

        # MESA Percentile
        self.mesaPercentileLabel = qt.QLabel("<b>MESA Percentile:</b> --")
        self.mesaPercentileLabel.setStyleSheet("font-size: 14px;")
        riskLayout.addWidget(self.mesaPercentileLabel)

        self.percentileDescriptionLabel = qt.QLabel("")
        riskLayout.addWidget(self.percentileDescriptionLabel)

        # Visualization buttons
        vizLayout = qt.QHBoxLayout()
        vizLayout.setSpacing(10)

        self.show3DButton = qt.QPushButton("Show 3D")
        self.show3DButton.setObjectName("secondaryButton")
        self.show3DButton.setMinimumHeight(38)
        self.show3DButton.clicked.connect(self.onShow3D)
        vizLayout.addWidget(self.show3DButton)

        self.showMPRButton = qt.QPushButton("MPR View")
        self.showMPRButton.setObjectName("secondaryButton")
        self.showMPRButton.setMinimumHeight(38)
        self.showMPRButton.clicked.connect(self.onShowMPR)
        vizLayout.addWidget(self.showMPRButton)

        self.showChartsButton = qt.QPushButton("Show Charts")
        self.showChartsButton.setObjectName("secondaryButton")
        self.showChartsButton.setMinimumHeight(38)
        self.showChartsButton.clicked.connect(self.onShowCharts)
        vizLayout.addWidget(self.showChartsButton)

        resultsLayout.addLayout(vizLayout)

        resultsLayout.addStretch(1)
        self.tabWidget.addTab(resultsTab, "3. Results")

    def setupReportTab(self):
        """Setup Tab 4: Generate Report"""
        reportTab = qt.QWidget()
        reportLayout = qt.QVBoxLayout(reportTab)

        # Output directory
        outputCollapsible = ctk.ctkCollapsibleButton()
        outputCollapsible.text = "Output Settings"
        reportLayout.addWidget(outputCollapsible)
        outputLayout = qt.QFormLayout(outputCollapsible)

        dirLayout = qt.QHBoxLayout()
        self.outputDirEdit = qt.QLineEdit()
        self.outputDirEdit.setPlaceholderText("Select output directory")
        dirLayout.addWidget(self.outputDirEdit)
        self.browseOutputButton = qt.QPushButton("Browse...")
        self.browseOutputButton.setMaximumWidth(80)
        self.browseOutputButton.clicked.connect(self.onBrowseOutput)
        dirLayout.addWidget(self.browseOutputButton)
        outputLayout.addRow("Directory:", dirLayout)

        # Report options
        optionsCollapsible = ctk.ctkCollapsibleButton()
        optionsCollapsible.text = "Report Options"
        reportLayout.addWidget(optionsCollapsible)
        optionsLayout = qt.QVBoxLayout(optionsCollapsible)

        self.includeScreenshotsCheck = qt.QCheckBox("Include Slice Screenshots")
        self.includeScreenshotsCheck.setChecked(True)
        optionsLayout.addWidget(self.includeScreenshotsCheck)

        self.include3DCheck = qt.QCheckBox("Include 3D View")
        self.include3DCheck.setChecked(True)
        optionsLayout.addWidget(self.include3DCheck)

        self.includeChartsCheck = qt.QCheckBox("Include Charts")
        self.includeChartsCheck.setChecked(True)
        optionsLayout.addWidget(self.includeChartsCheck)

        self.includePercentileChartCheck = qt.QCheckBox("Include MESA Percentile Chart")
        self.includePercentileChartCheck.setChecked(True)
        optionsLayout.addWidget(self.includePercentileChartCheck)

        # Generate button
        self.generateReportButton = qt.QPushButton("  Generate PDF Report")
        self.generateReportButton.setObjectName("primaryButton")
        self.generateReportButton.setMinimumHeight(44)
        self.generateReportButton.clicked.connect(self.onGenerateReport)
        reportLayout.addWidget(self.generateReportButton)

        # Status
        self.reportStatusLabel = qt.QLabel("")
        reportLayout.addWidget(self.reportStatusLabel)

        reportLayout.addStretch(1)
        self.tabWidget.addTab(reportTab, "4. Report")

    # ==================== Event Handlers ====================

    def onInstallDependencies(self):
        """Install required Python packages"""
        self.depsStatusLabel.setText("Installing dependencies...")
        slicer.app.processEvents()

        try:
            import pip
            slicer.util.pip_install('reportlab')
            slicer.util.pip_install('matplotlib')
            self.depsStatusLabel.setText("<span style='color:green'>Dependencies installed successfully!</span>")
        except Exception as e:
            self.depsStatusLabel.setText(f"<span style='color:red'>Error: {str(e)}</span>")

    def onBrowseLogo(self):
        """Browse for company logo"""
        filePath = qt.QFileDialog.getOpenFileName(self.parent, "Select Logo", "", "Images (*.png *.jpg *.jpeg)")
        if filePath:
            self.logoPathEdit.setText(filePath)

    def onBrowseOutput(self):
        """Browse for output directory"""
        dirPath = qt.QFileDialog.getExistingDirectory(self.parent, "Select Output Directory")
        if dirPath:
            self.outputDirEdit.setText(dirPath)

    def saveSettings(self):
        """Save settings to QSettings"""
        self.settings.setValue("threshold", self.thresholdSpinBox.value)
        self.settings.setValue("calibrationFactor", self.calibrationFactorSpinBox.value)
        self.settings.setValue("logoPath", self.logoPathEdit.text)
        self.settings.setValue("companyDesc", self.companyDescEdit.text)
        slicer.util.infoDisplay("Settings saved successfully!")

    def loadSettings(self):
        """Load settings from QSettings"""
        if self.settings.contains("threshold"):
            self.thresholdSpinBox.setValue(int(self.settings.value("threshold")))
        if self.settings.contains("calibrationFactor"):
            self.calibrationFactorSpinBox.setValue(float(self.settings.value("calibrationFactor")))
        if self.settings.contains("logoPath"):
            self.logoPathEdit.setText(self.settings.value("logoPath"))
        if self.settings.contains("companyDesc"):
            self.companyDescEdit.setText(self.settings.value("companyDesc"))

    def onVolumeSelected(self, node):
        """Handle volume selection"""
        if node:
            self.setupStatusLabel.setText(f"<span style='color:green'>Volume selected: {node.GetName()}</span>")
            self.tabWidget.setTabEnabled(2, True)

            # Create segmentations for each territory if not exist
            self.createTerritorySegmentations(node)
        else:
            self.setupStatusLabel.setText("")
            self.tabWidget.setTabEnabled(2, False)

    def onApplyLayout(self):
        """Apply optimal cardiac layout"""
        volumeNode = self.volumeSelector.currentNode()
        if not volumeNode:
            slicer.util.warningDisplay("Please select a volume first")
            return

        self.logic.setupOptimalView(volumeNode, self.enableIntersectionsCheck.isChecked())
        self.setupStatusLabel.setText("<span style='color:green'>Layout applied successfully!</span>")

    def createTerritorySegmentations(self, volumeNode):
        """Create segmentation nodes for each territory"""
        for abbrev, info in self.TERRITORIES.items():
            # Check if node exists and is still valid in the scene
            existingNode = self.segmentationsByTerritory[abbrev]
            nodeIsValid = (existingNode is not None and
                          slicer.mrmlScene.IsNodePresent(existingNode))

            if not nodeIsValid:
                # Reset reference if node was deleted
                self.segmentationsByTerritory[abbrev] = None

                # Try to find existing node by name (in case it was recreated)
                existingNodes = slicer.util.getNodesByClass("vtkMRMLSegmentationNode")
                for node in existingNodes:
                    if node.GetName() == f"CoronaryCalcium_{abbrev}":
                        self.segmentationsByTerritory[abbrev] = node
                        break

                # If still not found, create new node
                if self.segmentationsByTerritory[abbrev] is None:
                    segNode = slicer.mrmlScene.AddNewNodeByClass(
                        "vtkMRMLSegmentationNode",
                        f"CoronaryCalcium_{abbrev}"
                    )
                    segNode.SetReferenceImageGeometryParameterFromVolumeNode(volumeNode)
                    segNode.CreateDefaultDisplayNodes()

                    # Create empty segment with territory color
                    segmentation = segNode.GetSegmentation()
                    segmentId = segmentation.AddEmptySegment(f"Calcium_{abbrev}")
                    segment = segmentation.GetSegment(segmentId)
                    segment.SetColor(*info['color'])

                    self.segmentationsByTerritory[abbrev] = segNode

        self.updateTerritorySummary()

    def onTerritoryChanged(self, button):
        """Handle territory selection change"""
        for abbrev, btn in self.territoryButtons.items():
            if btn == button:
                self.currentTerritory = abbrev
                info = self.TERRITORIES[abbrev]
                self.currentTerritoryLabel.setText(f"<b>Current Territory:</b> {abbrev} - {info['name']}")

                # If Click & Grow is active, restart it for the new territory
                if self.activeToolName == 'clickgrow':
                    self.logic.stopClickGrowMode()
                    volumeNode = self.volumeSelector.currentNode()
                    segNode = self.segmentationsByTerritory[abbrev]
                    if volumeNode and segNode:
                        threshold = self.thresholdSpinBox.value
                        self.logic.startClickGrowMode(volumeNode, segNode, threshold, abbrev,
                                                      self.onClickGrowComplete, self.segmentationsByTerritory)
                        territoryColor = info['colorHex']
                        self.toolInstructionsLabel.setText(
                            f"<b style='color:{territoryColor}'>Click & Grow active for {abbrev}</b><br>"
                            f"<i>Click on calcium deposits to segment them. Change territory above to assign to different vessel.</i>"
                        )

                # Update segment editor if active
                if self.segmentEditorWidget and self.segmentationsByTerritory[abbrev]:
                    self.segmentEditorWidget.setSegmentationNode(self.segmentationsByTerritory[abbrev])
                break

    def deactivateAllTools(self):
        """Deactivate all segmentation tools"""
        if self.activeToolName == 'clickgrow':
            self.logic.stopClickGrowMode()
        elif self.activeToolName == 'brush':
            self.logic.stopPaintMode()
        elif self.activeToolName == 'erase':
            self.logic.stopEraseMode()

        self.activeToolName = None
        self.clickGrowButton.setChecked(False)
        self.brushMethodButton.setChecked(False)
        self.eraseButton.setChecked(False)

    def onClickGrowToggle(self):
        """Toggle Click & Grow mode"""
        if self.clickGrowButton.isChecked():
            # Deactivate other tools first
            if self.activeToolName and self.activeToolName != 'clickgrow':
                self.deactivateAllTools()

            volumeNode = self.volumeSelector.currentNode()
            if not volumeNode:
                slicer.util.warningDisplay("Please select a volume first")
                self.clickGrowButton.setChecked(False)
                return

            # Ensure segmentation nodes exist (recreate if deleted)
            self.createTerritorySegmentations(volumeNode)
            segNode = self.segmentationsByTerritory[self.currentTerritory]

            if not segNode:
                slicer.util.warningDisplay("Error creating segmentation node")
                self.clickGrowButton.setChecked(False)
                return

            self.activeToolName = 'clickgrow'
            self.brushMethodButton.setChecked(False)
            self.eraseButton.setChecked(False)

            threshold = self.thresholdSpinBox.value
            self.logic.startClickGrowMode(volumeNode, segNode, threshold, self.currentTerritory,
                                          self.onClickGrowComplete, self.segmentationsByTerritory)

            territoryColor = self.TERRITORIES[self.currentTerritory]['colorHex']
            self.toolInstructionsLabel.setText(
                f"<b style='color:{territoryColor}'>Click & Grow active for {self.currentTerritory}</b><br>"
                f"<i>Click on calcium deposits to segment them. Change territory above to assign to different vessel.</i>"
            )
        else:
            self.logic.stopClickGrowMode()
            self.activeToolName = None
            self.toolInstructionsLabel.setText("<i>Select a tool to start segmenting calcium deposits</i>")

    def onBrushToggle(self):
        """Toggle Brush paint mode"""
        if self.brushMethodButton.isChecked():
            # Deactivate other tools first
            if self.activeToolName and self.activeToolName != 'brush':
                self.deactivateAllTools()

            volumeNode = self.volumeSelector.currentNode()
            if not volumeNode:
                slicer.util.warningDisplay("Please select a volume first")
                self.brushMethodButton.setChecked(False)
                return

            # Ensure segmentation nodes exist (recreate if deleted)
            self.createTerritorySegmentations(volumeNode)
            segNode = self.segmentationsByTerritory[self.currentTerritory]

            if not segNode:
                slicer.util.warningDisplay("Error creating segmentation node")
                self.brushMethodButton.setChecked(False)
                return

            self.activeToolName = 'brush'
            self.clickGrowButton.setChecked(False)
            self.eraseButton.setChecked(False)

            self.logic.startPaintMode(segNode, self.currentTerritory, volumeNode, self.thresholdSpinBox.value)

            territoryColor = self.TERRITORIES[self.currentTerritory]['colorHex']
            self.toolInstructionsLabel.setText(
                f"<b style='color:{territoryColor}'>Paint mode active for {self.currentTerritory}</b><br>"
                f"<i>Paint over calcium deposits to add them to segmentation.</i>"
            )
        else:
            self.logic.stopPaintMode()
            self.activeToolName = None
            self.toolInstructionsLabel.setText("<i>Select a tool to start segmenting calcium deposits</i>")
            self.updateTerritorySummary()
            self.checkEnableResultsTab()

    def onEraseToggle(self):
        """Toggle Erase mode"""
        if self.eraseButton.isChecked():
            # Deactivate other tools first
            if self.activeToolName and self.activeToolName != 'erase':
                self.deactivateAllTools()

            volumeNode = self.volumeSelector.currentNode()
            if not volumeNode:
                slicer.util.warningDisplay("Please select a volume first")
                self.eraseButton.setChecked(False)
                return

            # Ensure segmentation nodes exist (recreate if deleted)
            self.createTerritorySegmentations(volumeNode)
            segNode = self.segmentationsByTerritory[self.currentTerritory]

            if not segNode:
                slicer.util.warningDisplay("Error creating segmentation node")
                self.eraseButton.setChecked(False)
                return

            self.activeToolName = 'erase'
            self.clickGrowButton.setChecked(False)
            self.brushMethodButton.setChecked(False)

            self.logic.startEraseMode(segNode, self.currentTerritory, volumeNode)

            self.toolInstructionsLabel.setText(
                f"<b style='color:orange'>Erase mode active for {self.currentTerritory}</b><br>"
                f"<i>Paint over areas to remove them from segmentation.</i>"
            )
        else:
            self.logic.stopEraseMode()
            self.activeToolName = None
            self.toolInstructionsLabel.setText("<i>Select a tool to start segmenting calcium deposits</i>")
            self.updateTerritorySummary()

    def onClickGrowComplete(self):
        """Called when click grow adds a lesion"""
        self.updateTerritorySummary()
        self.checkEnableResultsTab()

    def onClearTerritory(self):
        """Clear segmentation for current territory"""
        segNode = self.segmentationsByTerritory[self.currentTerritory]
        if segNode:
            segmentation = segNode.GetSegmentation()
            segmentation.RemoveAllSegments()
            # Recreate empty segment
            segmentId = segmentation.AddEmptySegment(f"Calcium_{self.currentTerritory}")
            segment = segmentation.GetSegment(segmentId)
            segment.SetColor(*self.TERRITORIES[self.currentTerritory]['color'])

        self.updateTerritorySummary()
        self.segmentationStatusLabel.setText(f"<span style='color:orange'>{self.currentTerritory} cleared</span>")

    def onClearAllTerritories(self):
        """Clear all territory segmentations"""
        for abbrev in self.TERRITORIES:
            segNode = self.segmentationsByTerritory[abbrev]
            if segNode:
                segmentation = segNode.GetSegmentation()
                segmentation.RemoveAllSegments()
                segmentId = segmentation.AddEmptySegment(f"Calcium_{abbrev}")
                segment = segmentation.GetSegment(segmentId)
                segment.SetColor(*self.TERRITORIES[abbrev]['color'])

        self.updateTerritorySummary()
        self.tabWidget.setTabEnabled(3, False)
        self.tabWidget.setTabEnabled(4, False)
        self.segmentationStatusLabel.setText("<span style='color:orange'>All territories cleared</span>")

    def updateTerritorySummary(self):
        """Update the territory summary table"""
        volumeNode = self.volumeSelector.currentNode()

        for i, (abbrev, info) in enumerate(self.TERRITORIES.items()):
            # Territory name with color
            nameItem = qt.QTableWidgetItem(f"{abbrev}")
            nameItem.setBackground(qt.QColor(info['colorHex']))
            if abbrev in ['LAD', 'RCA']:
                nameItem.setForeground(qt.QColor('white'))
            self.territorySummaryTable.setItem(i, 0, nameItem)

            # Check segmentation status
            segNode = self.segmentationsByTerritory.get(abbrev)
            voxelCount = 0
            status = "Empty"

            if segNode and volumeNode:
                try:
                    segmentation = segNode.GetSegmentation()
                    segmentId = f"Calcium_{abbrev}"  # Use ID directly
                    segment = segmentation.GetSegment(segmentId)

                    if segment:
                        # Use arrayFromSegmentBinaryLabelmap for direct access
                        try:
                            segArray = slicer.util.arrayFromSegmentBinaryLabelmap(
                                segNode, segmentId, volumeNode
                            )
                            if segArray is not None:
                                voxelCount = np.sum(segArray > 0)
                                if voxelCount > 0:
                                    status = "Segmented"
                        except Exception:
                            # Fallback: check if segment has any representation
                            if segment.GetRepresentation("Binary labelmap") or segment.GetRepresentation("Closed surface"):
                                status = "Segmented"
                                voxelCount = 1
                except Exception as e:
                    print(f"Error updating summary for {abbrev}: {e}")

            # Display voxel count or dash
            countText = str(voxelCount) if voxelCount > 0 else "-"
            self.territorySummaryTable.setItem(i, 1, qt.QTableWidgetItem(countText))

            # Status with color
            statusItem = qt.QTableWidgetItem(status)
            if status == "Segmented":
                statusItem.setForeground(qt.QColor('green'))
            else:
                statusItem.setForeground(qt.QColor('gray'))
            self.territorySummaryTable.setItem(i, 2, statusItem)

    def checkEnableResultsTab(self):
        """Check if any territory has segmentation and enable results tab"""
        hasSegmentation = False
        for segNode in self.segmentationsByTerritory.values():
            if segNode:
                segmentation = segNode.GetSegmentation()
                if segmentation.GetNumberOfSegments() > 0:
                    hasSegmentation = True
                    break

        self.tabWidget.setTabEnabled(3, hasSegmentation)

    def getPatientInfo(self):
        """Get current patient information"""
        self.patientInfo['name'] = self.patientNameEdit.text
        self.patientInfo['id'] = self.patientIdEdit.text
        self.patientInfo['sex'] = 'M' if self.maleRadio.isChecked() else 'F'
        age = self.patientAgeSpinBox.value
        self.patientInfo['age'] = age if age > 0 else None
        self.patientInfo['ethnicity'] = self.ethnicityComboBox.currentText
        self.patientInfo['date'] = datetime.now().strftime("%Y-%m-%d")
        return self.patientInfo

    def onCalculate(self):
        """Calculate Agatston scores for all territories"""
        volumeNode = self.volumeSelector.currentNode()
        if not volumeNode:
            slicer.util.warningDisplay("Please select a volume first")
            return

        patientInfo = self.getPatientInfo()

        # Calculate scores
        self.currentResults = self.logic.calculateAllTerritoryScores(
            volumeNode,
            self.segmentationsByTerritory,
            patientInfo,
            self.thresholdSpinBox.value,
            self.calibrationFactorSpinBox.value
        )

        # Update results table
        self.updateResultsTable()

        # Update risk classification
        self.updateRiskClassification()

        # Enable report tab
        self.tabWidget.setTabEnabled(4, True)

    def updateResultsTable(self):
        """Update the results table with calculated scores"""
        if not self.currentResults:
            return

        territories = list(self.TERRITORIES.keys()) + ['Total']

        for i, territory in enumerate(territories):
            result = self.currentResults.get(territory, {})

            # Territory name
            if territory == 'Total':
                nameItem = qt.QTableWidgetItem("TOTAL")
                nameItem.setBackground(qt.QColor('#333333'))
                nameItem.setForeground(qt.QColor('white'))
            else:
                nameItem = qt.QTableWidgetItem(territory)
                nameItem.setBackground(qt.QColor(self.TERRITORIES[territory]['colorHex']))
                if territory in ['LAD', 'RCA']:
                    nameItem.setForeground(qt.QColor('white'))

            self.resultsTable.setItem(i, 0, nameItem)

            # Score
            score = result.get('agatston_score', 0)
            self.resultsTable.setItem(i, 1, qt.QTableWidgetItem(f"{score:.1f}"))

            # Volume
            volume = result.get('total_volume_mm3', 0)
            self.resultsTable.setItem(i, 2, qt.QTableWidgetItem(f"{volume:.1f}"))

            # Mass
            mass = result.get('equivalent_mass_mg', 0)
            self.resultsTable.setItem(i, 3, qt.QTableWidgetItem(f"{mass:.1f}"))

            # Lesions
            lesions = result.get('num_lesions', 0)
            self.resultsTable.setItem(i, 4, qt.QTableWidgetItem(str(lesions)))

    def updateRiskClassification(self):
        """Update risk classification display"""
        if not self.currentResults or 'Total' not in self.currentResults:
            return

        totalResult = self.currentResults['Total']

        # MESA/ACC Category
        category = totalResult.get('mesa_category', {})
        categoryName = category.get('category', '--')
        categoryDesc = category.get('description', '')
        risk = category.get('risk', '')

        # Color coding based on risk
        riskColors = {
            'Very Low': '#00AA00',
            'Low': '#88CC00',
            'Low-Moderate': '#FFCC00',
            'Moderate-High': '#FF8800',
            'High': '#FF0000'
        }
        color = riskColors.get(risk, '#000000')

        self.mesaCategoryLabel.setText(
            f"<b>MESA/ACC Category:</b> <span style='color:{color}; font-size:16px;'>{categoryName}</span> ({categoryDesc})"
        )
        self.riskDescriptionLabel.setText(f"<b>Risk Level:</b> {risk}")

        # MESA Percentile
        percentile = totalResult.get('mesa_percentile', {})
        percentileValue = percentile.get('percentile', '--')
        comparison = percentile.get('comparison', '')

        if percentileValue != '--':
            self.mesaPercentileLabel.setText(
                f"<b>MESA Percentile:</b> <span style='font-size:16px;'>{percentileValue}th</span>"
            )
            self.percentileDescriptionLabel.setText(f"{comparison}")
        else:
            self.mesaPercentileLabel.setText("<b>MESA Percentile:</b> Age required for calculation")
            self.percentileDescriptionLabel.setText("")

    def onShow3D(self):
        """Show 3D visualization of calcium by territory"""
        if not self.currentResults:
            slicer.util.warningDisplay("Please calculate scores first")
            return

        volumeNode = self.volumeSelector.currentNode()
        self.logic.create3DVisualization(self.segmentationsByTerritory, self.TERRITORIES, volumeNode)

    def onShowMPR(self):
        """Restore standard MPR (multiplanar) view"""
        volumeNode = self.volumeSelector.currentNode()
        self.logic.showMPRView(volumeNode)

    def onShowCharts(self):
        """Show analysis charts"""
        if not self.currentResults:
            slicer.util.warningDisplay("Please calculate scores first")
            return

        patientInfo = self.getPatientInfo()
        self.logic.createCharts(self.currentResults, self.TERRITORIES, patientInfo)

    def onGenerateReport(self):
        """Generate PDF report"""
        outputDir = self.outputDirEdit.text
        if not outputDir or not os.path.isdir(outputDir):
            slicer.util.warningDisplay("Please select a valid output directory")
            return

        if not self.currentResults:
            slicer.util.warningDisplay("Please calculate scores first")
            return

        self.reportStatusLabel.setText("Generating report...")
        slicer.app.processEvents()

        try:
            patientInfo = self.getPatientInfo()

            options = {
                'includeScreenshots': self.includeScreenshotsCheck.isChecked(),
                'include3D': self.include3DCheck.isChecked(),
                'includeCharts': self.includeChartsCheck.isChecked(),
                'includePercentileChart': self.includePercentileChartCheck.isChecked(),
                'logoPath': self.logoPathEdit.text,
                'companyDesc': self.companyDescEdit.text,
                'threshold': self.thresholdSpinBox.value,
                'calibrationFactor': self.calibrationFactorSpinBox.value
            }

            reportPath = self.logic.generatePDFReport(
                outputDir,
                self.currentResults,
                patientInfo,
                self.TERRITORIES,
                self.segmentationsByTerritory,
                self.volumeSelector.currentNode(),
                options
            )

            self.reportStatusLabel.setText(f"<span style='color:green'>Report saved: {reportPath}</span>")

            # Open report
            qt.QDesktopServices.openUrl(qt.QUrl.fromLocalFile(reportPath))

        except Exception as e:
            self.reportStatusLabel.setText(f"<span style='color:red'>Error: {str(e)}</span>")
            import traceback
            traceback.print_exc()


#
# CoronaryCaScoreLogic
#

class CoronaryCaScoreLogic(ScriptedLoadableModuleLogic):
    """Logic for coronary calcium scoring calculations"""

    # MESA/ACC Risk Categories
    MESA_CATEGORIES = [
        {'min': 0, 'max': 0, 'category': 'Zero', 'description': '0 AU', 'risk': 'Very Low'},
        {'min': 1, 'max': 10, 'category': 'Minimal', 'description': '1-10 AU', 'risk': 'Low'},
        {'min': 11, 'max': 100, 'category': 'Mild', 'description': '11-100 AU', 'risk': 'Low-Moderate'},
        {'min': 101, 'max': 400, 'category': 'Moderate', 'description': '101-400 AU', 'risk': 'Moderate-High'},
        {'min': 401, 'max': float('inf'), 'category': 'Severe', 'description': '>400 AU', 'risk': 'High'}
    ]

    # MESA Percentiles data (simplified - in production use full JSON file)
    MESA_PERCENTILES = {
        'Male': {
            'Caucasian': {
                '45-54': {'25': 0, '50': 4, '75': 95, '90': 297},
                '55-64': {'25': 1, '50': 38, '75': 239, '90': 616},
                '65-74': {'25': 12, '50': 130, '75': 497, '90': 1175},
                '75-84': {'25': 57, '50': 309, '75': 892, '90': 1982}
            },
            'African-American': {
                '45-54': {'25': 0, '50': 0, '75': 16, '90': 110},
                '55-64': {'25': 0, '50': 5, '75': 78, '90': 296},
                '65-74': {'25': 0, '50': 35, '75': 229, '90': 651},
                '75-84': {'25': 8, '50': 115, '75': 509, '90': 1241}
            },
            'Hispanic': {
                '45-54': {'25': 0, '50': 0, '75': 42, '90': 172},
                '55-64': {'25': 0, '50': 14, '75': 133, '90': 403},
                '65-74': {'25': 3, '50': 76, '75': 341, '90': 848},
                '75-84': {'25': 27, '50': 197, '75': 672, '90': 1533}
            },
            'Chinese': {
                '45-54': {'25': 0, '50': 0, '75': 18, '90': 96},
                '55-64': {'25': 0, '50': 3, '75': 66, '90': 247},
                '65-74': {'25': 0, '50': 34, '75': 200, '90': 566},
                '75-84': {'25': 10, '50': 114, '75': 462, '90': 1131}
            }
        },
        'Female': {
            'Caucasian': {
                '45-54': {'25': 0, '50': 0, '75': 3, '90': 40},
                '55-64': {'25': 0, '50': 1, '75': 38, '90': 152},
                '65-74': {'25': 0, '50': 17, '75': 131, '90': 392},
                '75-84': {'25': 2, '50': 75, '75': 321, '90': 829}
            },
            'African-American': {
                '45-54': {'25': 0, '50': 0, '75': 0, '90': 12},
                '55-64': {'25': 0, '50': 0, '75': 10, '90': 65},
                '65-74': {'25': 0, '50': 3, '75': 52, '90': 199},
                '75-84': {'25': 0, '50': 22, '75': 148, '90': 482}
            },
            'Hispanic': {
                '45-54': {'25': 0, '50': 0, '75': 0, '90': 18},
                '55-64': {'25': 0, '50': 0, '75': 16, '90': 79},
                '65-74': {'25': 0, '50': 6, '75': 67, '90': 230},
                '75-84': {'25': 0, '50': 34, '75': 183, '90': 538}
            },
            'Chinese': {
                '45-54': {'25': 0, '50': 0, '75': 0, '90': 6},
                '55-64': {'25': 0, '50': 0, '75': 4, '90': 36},
                '65-74': {'25': 0, '50': 1, '75': 27, '90': 115},
                '75-84': {'25': 0, '50': 13, '75': 93, '90': 313}
            }
        }
    }

    def __init__(self):
        ScriptedLoadableModuleLogic.__init__(self)
        self.roiNode = None
        self.clickGrowObserver = None

    def setupOptimalView(self, volumeNode, enableIntersections=True):
        """Set up optimal view for cardiac calcium scoring"""
        # Set layout to Four-Up
        layoutManager = slicer.app.layoutManager()
        layoutManager.setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutFourUpView)

        # Set volume in all slice views
        sliceCompositeNodes = slicer.util.getNodesByClass("vtkMRMLSliceCompositeNode")
        for sliceCompositeNode in sliceCompositeNodes:
            sliceCompositeNode.SetBackgroundVolumeID(volumeNode.GetID())

        # Apply cardiac window/level (W=1500, L=300)
        displayNode = volumeNode.GetDisplayNode()
        if displayNode:
            displayNode.AutoWindowLevelOff()
            displayNode.SetWindow(1500)
            displayNode.SetLevel(300)

        # Fit slices to volume
        for sliceView in ['Red', 'Yellow', 'Green']:
            sliceLogic = slicer.app.layoutManager().sliceWidget(sliceView).sliceLogic()
            sliceLogic.FitSliceToAll()

        # Enable slice intersections
        if enableIntersections:
            for sliceNode in slicer.util.getNodesByClass("vtkMRMLSliceNode"):
                sliceNode.SetSliceVisible(True)

    def createROI(self, volumeNode):
        """Create ROI box for segmentation"""
        # Remove existing ROI if any
        if self.roiNode:
            slicer.mrmlScene.RemoveNode(self.roiNode)

        # Create new ROI
        self.roiNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLMarkupsROINode", "CalciumROI")

        # Set ROI to volume center with appropriate size
        bounds = [0] * 6
        volumeNode.GetBounds(bounds)

        center = [(bounds[0] + bounds[1]) / 2, (bounds[2] + bounds[3]) / 2, (bounds[4] + bounds[5]) / 2]
        size = [50, 50, 50]  # 50mm cube as default

        self.roiNode.SetCenter(center)
        self.roiNode.SetSize(size)

        return self.roiNode

    def applyThresholdInROI(self, volumeNode, segmentationNode, threshold, territory):
        """Apply threshold segmentation within ROI"""
        if not self.roiNode:
            slicer.util.warningDisplay("Please create an ROI first")
            return

        # Get segment using ID directly
        segmentation = segmentationNode.GetSegmentation()
        segmentId = f"Calcium_{territory}"
        segment = segmentation.GetSegment(segmentId)

        if not segment:
            return

        # Create segment editor widget instance
        segmentEditorWidget = slicer.qMRMLSegmentEditorWidget()
        segmentEditorWidget.setMRMLScene(slicer.mrmlScene)

        segmentEditorNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLSegmentEditorNode")
        segmentEditorWidget.setMRMLSegmentEditorNode(segmentEditorNode)
        segmentEditorWidget.setSegmentationNode(segmentationNode)
        segmentEditorWidget.setSourceVolumeNode(volumeNode)
        segmentEditorWidget.setCurrentSegmentID(segmentId)

        # Apply threshold effect
        segmentEditorWidget.setActiveEffectByName("Threshold")
        effect = segmentEditorWidget.activeEffect()

        if effect:
            effect.setParameter("MinimumThreshold", str(threshold))
            effect.setParameter("MaximumThreshold", "3000")
            effect.self().onApply()

        # Apply masking with ROI
        segmentEditorWidget.setActiveEffectByName("Mask volume")
        # ... ROI masking logic

        # Cleanup
        slicer.mrmlScene.RemoveNode(segmentEditorNode)

    def startClickGrowMode(self, volumeNode, segmentationNode, threshold, territory, callback, allSegmentations=None):
        """Start click and grow mode with automatic 3D region growing on click"""
        self.clickGrowSegmentationNode = segmentationNode
        self.clickGrowVolumeNode = volumeNode
        self.clickGrowThreshold = threshold
        self.clickGrowTerritory = territory
        self.clickGrowCallback = callback
        self.clickGrowAllSegmentations = allSegmentations or {}  # For overwrite functionality
        self.clickGrowObserverTag = None

        # Store volume array for fast access
        self.clickGrowVolumeArray = slicer.util.arrayFromVolume(volumeNode)

        # Set up markup fiducial for clicking
        self.clickGrowFiducialNode = slicer.mrmlScene.AddNewNodeByClass(
            "vtkMRMLMarkupsFiducialNode", f"ClickGrow_{territory}_Seeds"
        )
        self.clickGrowFiducialNode.CreateDefaultDisplayNodes()

        # Make fiducials small and temporary
        displayNode = self.clickGrowFiducialNode.GetDisplayNode()
        displayNode.SetGlyphScale(1.5)
        displayNode.SetSelectedColor(0, 1, 0)  # Green

        # Set interaction mode to place fiducials
        selectionNode = slicer.mrmlScene.GetNodeByID("vtkMRMLSelectionNodeSingleton")
        selectionNode.SetActivePlaceNodeID(self.clickGrowFiducialNode.GetID())
        interactionNode = slicer.mrmlScene.GetNodeByID("vtkMRMLInteractionNodeSingleton")
        interactionNode.SetCurrentInteractionMode(interactionNode.Place)

        # Add observer to detect when fiducial is placed
        self.clickGrowObserverTag = self.clickGrowFiducialNode.AddObserver(
            slicer.vtkMRMLMarkupsNode.PointPositionDefinedEvent,
            self.onFiducialPlaced
        )

        print(f"Click & Grow Mode: Active for {territory} with threshold >= {threshold} HU")
        print("Click on any calcification to segment it automatically")

    def onFiducialPlaced(self, caller, event):
        """Called when user clicks to place a fiducial - perform region growing"""
        try:
            fiducialNode = caller
            numFiducials = fiducialNode.GetNumberOfControlPoints()

            if numFiducials == 0:
                return

            # Get the position of the last placed fiducial
            lastIndex = numFiducials - 1
            seedPos_RAS = [0, 0, 0]
            fiducialNode.GetNthControlPointPosition(lastIndex, seedPos_RAS)

            print(f"\nClick #{lastIndex + 1}: Growing from seed at RAS: {seedPos_RAS}")

            # Convert RAS to IJK
            volumeNode = self.clickGrowVolumeNode
            rasToIjkMatrix = vtk.vtkMatrix4x4()
            volumeNode.GetRASToIJKMatrix(rasToIjkMatrix)

            seedPos_IJK = [0, 0, 0, 1]
            rasPoint = [seedPos_RAS[0], seedPos_RAS[1], seedPos_RAS[2], 1]
            ijkPoint = rasToIjkMatrix.MultiplyPoint(rasPoint)
            seedPos_IJK = [int(round(ijkPoint[0])), int(round(ijkPoint[1])), int(round(ijkPoint[2]))]

            print(f"Seed IJK: {seedPos_IJK}")

            # Perform 3D region growing from this seed
            self.performRegionGrowing3D(seedPos_IJK)

            # Remove the fiducial point after processing
            fiducialNode.RemoveNthControlPoint(lastIndex)

            # Force placement mode to stay active
            qt.QTimer.singleShot(100, self.reactivatePlacementMode)

        except Exception as e:
            print(f"Error in onFiducialPlaced: {str(e)}")
            import traceback
            traceback.print_exc()

    def reactivatePlacementMode(self):
        """Reactivate placement mode for continuous clicking"""
        try:
            if hasattr(self, 'clickGrowFiducialNode') and self.clickGrowFiducialNode:
                selectionNode = slicer.mrmlScene.GetNodeByID("vtkMRMLSelectionNodeSingleton")
                selectionNode.SetActivePlaceNodeID(self.clickGrowFiducialNode.GetID())
                interactionNode = slicer.mrmlScene.GetNodeByID("vtkMRMLInteractionNodeSingleton")
                interactionNode.SetCurrentInteractionMode(interactionNode.Place)
                print(f"-> Ready for next click on {self.clickGrowTerritory}...")
        except Exception as e:
            print(f"Error reactivating placement mode: {str(e)}")

    def performRegionGrowing3D(self, seedIJK):
        """Perform 3D region growing from seed point using connected components"""
        volumeArray = self.clickGrowVolumeArray
        threshold = self.clickGrowThreshold

        # Check if seed is valid
        dims = volumeArray.shape
        i, j, k = seedIJK

        if not (0 <= k < dims[0] and 0 <= j < dims[1] and 0 <= i < dims[2]):
            print(f"Seed outside volume bounds: {seedIJK}")
            return

        seedValue = volumeArray[k, j, i]
        print(f"Seed HU value: {seedValue:.1f}")

        if seedValue < threshold:
            print(f"Seed HU ({seedValue:.1f}) is below threshold ({threshold}). Not segmenting.")
            slicer.util.warningDisplay(
                f"Clicked point has HU = {seedValue:.1f}\n"
                f"This is below the threshold ({threshold} HU).\n\n"
                f"Please click on a calcification (bright white area)."
            )
            return

        # 3D connected component labeling with threshold
        from scipy import ndimage

        print(f"Performing 3D region growing from seed [{i}, {j}, {k}]...")

        # Create binary mask: voxels >= threshold
        binaryMask = (volumeArray >= threshold).astype(np.uint8)

        # Label all connected components in 3D (26-connectivity for full 3D)
        labeled, numFeatures = ndimage.label(binaryMask, structure=ndimage.generate_binary_structure(3, 3))

        # Find which label the seed belongs to
        seedLabel = labeled[k, j, i]

        if seedLabel == 0:
            print("Seed is not part of any connected component above threshold")
            slicer.util.warningDisplay(
                f"No connected calcification found at this point.\n\n"
                f"The clicked voxel may be isolated or below threshold."
            )
            return

        # Extract only the connected component containing the seed
        componentMask = (labeled == seedLabel).astype(np.uint8)

        # Count voxels and calculate volume
        numVoxels = np.sum(componentMask)
        volumeNode = self.clickGrowVolumeNode
        spacing = volumeNode.GetSpacing()
        voxelVolume = spacing[0] * spacing[1] * spacing[2]
        lesionVolume = numVoxels * voxelVolume

        print(f"Segmented lesion: {numVoxels} voxels, {lesionVolume:.1f} mm3")

        # Remove this lesion from other territories (overwrite functionality)
        print(f"Checking overwrite for other territories...")
        if hasattr(self, 'clickGrowAllSegmentations') and self.clickGrowAllSegmentations:
            print(f"  Available territories in dict: {list(self.clickGrowAllSegmentations.keys())}")
            for otherTerritory, otherSegNode in self.clickGrowAllSegmentations.items():
                if otherTerritory != self.clickGrowTerritory:
                    print(f"  Checking {otherTerritory}: segNode={otherSegNode is not None}")
                    if otherSegNode:
                        try:
                            otherSegmentation = otherSegNode.GetSegmentation()
                            # Use segment ID directly (not segment name)
                            otherSegmentId = f"Calcium_{otherTerritory}"
                            print(f"    Using segment ID: '{otherSegmentId}'")
                            if otherSegmentId:
                                otherSegment = otherSegmentation.GetSegment(otherSegmentId)
                                if otherSegment:
                                    # Check if segment has any representation
                                    hasLabelmap = otherSegment.GetRepresentation("Binary labelmap") is not None
                                    print(f"    Has Binary labelmap: {hasLabelmap}")
                                    if hasLabelmap:
                                        try:
                                            otherArray = slicer.util.arrayFromSegmentBinaryLabelmap(
                                                otherSegNode, otherSegmentId, volumeNode
                                            )
                                            print(f"    otherArray shape: {otherArray.shape if otherArray is not None else 'None'}")
                                            print(f"    componentMask shape: {componentMask.shape}")

                                            if otherArray is not None:
                                                # Ensure arrays have same shape
                                                if otherArray.shape != componentMask.shape:
                                                    print(f"    Shape mismatch! Resizing...")
                                                    # Resize otherArray to match componentMask if needed
                                                    if otherArray.shape != componentMask.shape:
                                                        # Create a full-size array
                                                        fullOtherArray = np.zeros_like(componentMask, dtype=np.uint8)
                                                        # Copy whatever fits
                                                        minShape = tuple(min(s1, s2) for s1, s2 in zip(otherArray.shape, componentMask.shape))
                                                        fullOtherArray[:minShape[0], :minShape[1], :minShape[2]] = \
                                                            otherArray[:minShape[0], :minShape[1], :minShape[2]]
                                                        otherArray = fullOtherArray

                                                otherVoxelCount = np.sum(otherArray > 0)
                                                print(f"    Other territory voxel count: {otherVoxelCount}")

                                                overlap = np.logical_and(otherArray > 0, componentMask > 0)
                                                overlapCount = np.sum(overlap)
                                                print(f"    Overlap voxel count: {overlapCount}")

                                                if overlapCount > 0:
                                                    # Remove the overlapping region
                                                    newOtherArray = np.logical_and(otherArray > 0, np.logical_not(componentMask)).astype(np.uint8)
                                                    newVoxelCount = np.sum(newOtherArray)
                                                    print(f"    After removal: {newVoxelCount} voxels")

                                                    # Update the segment binary labelmap
                                                    slicer.util.updateSegmentBinaryLabelmapFromArray(
                                                        newOtherArray, otherSegNode, otherSegmentId, volumeNode
                                                    )
                                                    print(f"  *** Removed {overlapCount} overlapping voxels from {otherTerritory}")

                                                    # Invalidate closed surface representation to force recalculation
                                                    otherSegment.RemoveRepresentation("Closed surface")

                                                    # Force segment and node to refresh display
                                                    otherSegmentation.Modified()
                                                    otherSegNode.Modified()

                                                    # Force display node update
                                                    displayNode = otherSegNode.GetDisplayNode()
                                                    if displayNode:
                                                        displayNode.Modified()

                                                    # Force all views to update
                                                    slicer.util.forceRenderAllViews()
                                                else:
                                                    print(f"    No overlap detected")
                                        except Exception as e:
                                            print(f"    Error reading/updating {otherTerritory}: {e}")
                                            import traceback
                                            traceback.print_exc()
                        except Exception as e:
                            print(f"  Error processing {otherTerritory}: {e}")
        else:
            print(f"  No allSegmentations dict available")

        # Add to territory segmentation
        segmentationNode = self.clickGrowSegmentationNode
        segmentation = segmentationNode.GetSegmentation()
        mainSegmentId = segmentation.GetNthSegmentID(0)

        # Ensure the segmentation has the correct geometry
        segmentationNode.SetReferenceImageGeometryParameterFromVolumeNode(volumeNode)

        # Get current segmentation array
        segment = segmentation.GetSegment(mainSegmentId)
        currentArray = None

        if segment and segment.GetRepresentation("Binary labelmap"):
            try:
                currentArray = slicer.util.arrayFromSegmentBinaryLabelmap(
                    segmentationNode, mainSegmentId, volumeNode
                )
            except Exception as e:
                print(f"Note: Could not export segment (first click?): {e}")

        # If no current array, create empty one
        if currentArray is None:
            imageData = volumeNode.GetImageData()
            dims = imageData.GetDimensions()
            currentArray = np.zeros((dims[2], dims[1], dims[0]), dtype=np.uint8)
            print("Initializing segmentation with first calcification")

        # Replace with new component (not OR - full replacement for clicked lesion)
        currentArray = np.logical_or(currentArray, componentMask).astype(np.uint8)

        # Update segmentation
        slicer.util.updateSegmentBinaryLabelmapFromArray(
            currentArray,
            segmentationNode,
            mainSegmentId,
            volumeNode
        )

        print(f"Successfully added calcification to {self.clickGrowTerritory}!")

        # Show notification
        slicer.util.showStatusMessage(f"Added lesion to {self.clickGrowTerritory}: {lesionVolume:.1f} mm3", 2000)

        # Call callback to update UI
        if self.clickGrowCallback:
            self.clickGrowCallback()

    def stopClickGrowMode(self):
        """Stop click and grow mode"""
        # Remove observer
        if hasattr(self, 'clickGrowFiducialNode') and self.clickGrowFiducialNode and hasattr(self, 'clickGrowObserverTag') and self.clickGrowObserverTag:
            self.clickGrowFiducialNode.RemoveObserver(self.clickGrowObserverTag)
            self.clickGrowObserverTag = None

        # Remove fiducial node
        if hasattr(self, 'clickGrowFiducialNode') and self.clickGrowFiducialNode:
            slicer.mrmlScene.RemoveNode(self.clickGrowFiducialNode)
            self.clickGrowFiducialNode = None

        # Stop placing mode
        interactionNode = slicer.mrmlScene.GetNodeByID("vtkMRMLInteractionNodeSingleton")
        interactionNode.SetCurrentInteractionMode(interactionNode.ViewTransform)

        # Clean up references
        self.clickGrowSegmentationNode = None
        self.clickGrowVolumeNode = None
        self.clickGrowThreshold = None
        self.clickGrowVolumeArray = None
        self.clickGrowCallback = None

        print("Click & Grow Mode: Stopped")

    def startPaintMode(self, segmentationNode, territory, volumeNode=None, threshold=130):
        """Start paint mode for manual segmentation with threshold"""
        try:
            # Create segment editor widget if needed
            if not hasattr(self, 'segmentEditorWidget') or self.segmentEditorWidget is None:
                self.segmentEditorWidget = slicer.qMRMLSegmentEditorWidget()
                self.segmentEditorWidget.setMRMLScene(slicer.mrmlScene)
                self.segmentEditorNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLSegmentEditorNode")
                self.segmentEditorWidget.setMRMLSegmentEditorNode(self.segmentEditorNode)

            self.segmentEditorWidget.setSegmentationNode(segmentationNode)
            if volumeNode:
                self.segmentEditorWidget.setSourceVolumeNode(volumeNode)

            # Get segment ID directly
            segmentation = segmentationNode.GetSegmentation()
            segmentId = f"Calcium_{territory}"
            if segmentation.GetSegment(segmentId):
                self.segmentEditorWidget.setCurrentSegmentID(segmentId)

            # Activate threshold paint effect
            self.segmentEditorWidget.setActiveEffectByName("Paint")
            effect = self.segmentEditorWidget.activeEffect()
            if effect:
                # Set threshold for paint
                effect.setParameter("BrushMinimumAbsoluteDiameter", "1")
                effect.setParameter("BrushMaximumAbsoluteDiameter", "10")

            print(f"Paint mode active for {territory}")
        except Exception as e:
            print(f"Error starting paint mode: {e}")

    def stopPaintMode(self):
        """Stop paint mode"""
        try:
            if hasattr(self, 'segmentEditorWidget') and self.segmentEditorWidget:
                self.segmentEditorWidget.setActiveEffect(None)
        except Exception as e:
            print(f"Error stopping paint mode: {e}")

    def startEraseMode(self, segmentationNode, territory, volumeNode=None):
        """Start erase mode"""
        try:
            # Create segment editor widget if needed
            if not hasattr(self, 'segmentEditorWidget') or self.segmentEditorWidget is None:
                self.segmentEditorWidget = slicer.qMRMLSegmentEditorWidget()
                self.segmentEditorWidget.setMRMLScene(slicer.mrmlScene)
                self.segmentEditorNode = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLSegmentEditorNode")
                self.segmentEditorWidget.setMRMLSegmentEditorNode(self.segmentEditorNode)

            self.segmentEditorWidget.setSegmentationNode(segmentationNode)
            if volumeNode:
                self.segmentEditorWidget.setSourceVolumeNode(volumeNode)

            # Get segment ID directly
            segmentation = segmentationNode.GetSegmentation()
            segmentId = f"Calcium_{territory}"
            if segmentation.GetSegment(segmentId):
                self.segmentEditorWidget.setCurrentSegmentID(segmentId)

            self.segmentEditorWidget.setActiveEffectByName("Erase")
            print(f"Erase mode active for {territory}")
        except Exception as e:
            print(f"Error starting erase mode: {e}")

    def stopEraseMode(self):
        """Stop erase mode"""
        try:
            if hasattr(self, 'segmentEditorWidget') and self.segmentEditorWidget:
                self.segmentEditorWidget.setActiveEffect(None)
        except Exception as e:
            print(f"Error stopping erase mode: {e}")

    def calculateAllTerritoryScores(self, volumeNode, segmentationsByTerritory, patientInfo, threshold, calibrationFactor=0.81):
        """Calculate Agatston scores for all territories"""
        results = {}

        totalAgatston = 0
        totalVolume = 0
        totalMass = 0
        totalLesions = 0

        print(f"\n=== Coronary Calcium Score Calculation ===")
        print(f"Threshold: {threshold} HU, Calibration Factor: {calibrationFactor}")

        # Calculate for each territory
        for territory, segNode in segmentationsByTerritory.items():
            if segNode:
                territoryResult = self.calculateTerritoryScore(volumeNode, segNode, territory, threshold, calibrationFactor)
                results[territory] = territoryResult

                totalAgatston += territoryResult.get('agatston_score', 0)
                totalVolume += territoryResult.get('total_volume_mm3', 0)
                totalMass += territoryResult.get('equivalent_mass_mg', 0)
                totalLesions += territoryResult.get('num_lesions', 0)
            else:
                results[territory] = self.getEmptyResult()

        print(f"\n=== TOTAL: Score={totalAgatston:.1f} AU, Volume={totalVolume:.1f} mm³, Mass={totalMass:.1f} mg ===\n")

        # Calculate total
        results['Total'] = {
            'agatston_score': totalAgatston,
            'total_volume_mm3': totalVolume,
            'equivalent_mass_mg': totalMass,
            'num_lesions': totalLesions,
            'mesa_category': self.getMESACategory(totalAgatston),
            'mesa_percentile': self.getMESAPercentile(totalAgatston, patientInfo)
        }

        return results

    def calculateTerritoryScore(self, volumeNode, segmentationNode, territory, threshold, calibrationFactor=0.81):
        """Calculate Agatston score for a single territory using standardized method

        This implementation follows the Agatston standard:
        - 2D slice-by-slice connected component labeling
        - Minimum 1 mm² area per lesion per slice
        - Density factor determined per slice
        - Handles overlapping slices (skips if spacing < 2.5mm)
        """
        try:
            from scipy import ndimage

            # Get segment array
            segmentId = segmentationNode.GetSegmentation().GetNthSegmentID(0)
            segmentArray = slicer.util.arrayFromSegmentBinaryLabelmap(segmentationNode, segmentId, volumeNode)
            volumeArray = slicer.util.arrayFromVolume(volumeNode)

            if segmentArray is None or not np.any(segmentArray):
                return self.getEmptyResult()

            # Get spacing for area/volume calculation
            spacing = volumeNode.GetSpacing()
            sliceArea = spacing[0] * spacing[1]  # mm² per pixel
            sliceSpacing = spacing[2]

            # Try to get slice thickness from DICOM metadata
            sliceThickness = None
            try:
                import pydicom
                storageNode = volumeNode.GetStorageNode()
                if storageNode:
                    filePath = storageNode.GetFileName()
                    if filePath:
                        dicom = pydicom.dcmread(filePath, stop_before_pixels=True)
                        if hasattr(dicom, 'SliceThickness'):
                            sliceThickness = float(dicom.SliceThickness)
                            print(f"  DICOM Slice Thickness: {sliceThickness:.1f}mm")
            except:
                pass

            if sliceThickness is None:
                sliceThickness = sliceSpacing

            # CRITICAL: Handle overlapping slices (standard Agatston uses 3mm increment)
            STANDARD_SLICE_INCREMENT = 3.0  # mm (Agatston standard)

            if sliceSpacing < 2.5:  # Overlapping slices detected
                sliceStep = max(1, int(round(STANDARD_SLICE_INCREMENT / sliceSpacing)))
                print(f"  {territory}: Overlapping slices (spacing={sliceSpacing:.1f}mm), using step={sliceStep}")
            else:
                sliceStep = 1

            # 2D slice-by-slice connected component labeling (Agatston standard)
            labeled_array = np.zeros_like(segmentArray, dtype=np.int32)
            lesion_counter = 0

            for sliceIdx in range(0, segmentArray.shape[0], sliceStep):
                sliceMask = segmentArray[sliceIdx, :, :]
                if np.any(sliceMask):
                    # Label lesions in this slice only (2D connectivity)
                    labeled_slice, n_lesions_in_slice = ndimage.label(sliceMask)
                    labeled_slice[labeled_slice > 0] += lesion_counter
                    labeled_array[sliceIdx, :, :] = labeled_slice
                    lesion_counter += n_lesions_in_slice

            num_lesions = lesion_counter

            # Agatston standard: minimum 1 mm² area per slice
            MIN_AREA_MM2_STANDARD = 1.0
            minPixelsFor1mm2 = int(np.ceil(MIN_AREA_MM2_STANDARD / sliceArea))
            MIN_PIXELS_PER_SLICE = max(2, minPixelsFor1mm2)
            MIN_AREA_MM2 = MIN_PIXELS_PER_SLICE * sliceArea

            # Calculate metrics
            totalAgatston = 0
            totalVolume = 0
            allDensities = []
            validLesionCount = 0

            for lesion_id in range(1, num_lesions + 1):
                lesionMask = (labeled_array == lesion_id)
                lesionVoxels = volumeArray[lesionMask]

                # Filter by threshold
                lesionVoxels = lesionVoxels[lesionVoxels >= threshold]
                if len(lesionVoxels) == 0:
                    continue

                allDensities.extend(lesionVoxels.tolist())

                # Calculate Agatston score slice-by-slice (STANDARD METHOD)
                lesionScore = 0
                lesionVolume = 0
                lesionHasValidSlice = False

                for sliceIdx in range(lesionMask.shape[0]):
                    sliceMask = lesionMask[sliceIdx, :, :]
                    if np.any(sliceMask):
                        # Get voxels in this slice and filter by threshold
                        sliceVoxels = volumeArray[sliceIdx, :, :][sliceMask]
                        sliceVoxels = sliceVoxels[sliceVoxels >= threshold]

                        if len(sliceVoxels) == 0:
                            continue

                        # Calculate area for this slice
                        pixelCount = len(sliceVoxels)
                        area_mm2 = pixelCount * sliceArea

                        # CRITICAL: Filter by minimum area PER LESION PER SLICE
                        if area_mm2 < MIN_AREA_MM2:
                            continue

                        lesionHasValidSlice = True

                        # Get max density in THIS SLICE (Agatston standard)
                        maxDensityInSlice = np.max(sliceVoxels)

                        # Determine density factor for THIS SLICE
                        if maxDensityInSlice >= 400:
                            densityFactor = 4
                        elif maxDensityInSlice >= 300:
                            densityFactor = 3
                        elif maxDensityInSlice >= 200:
                            densityFactor = 2
                        else:  # 130-199
                            densityFactor = 1

                        # Calculate score for this slice
                        lesionScore += area_mm2 * densityFactor

                        # Volume calculation with slice step adjustment
                        sliceVolumeContribution = pixelCount * sliceArea * (sliceSpacing * sliceStep)
                        lesionVolume += sliceVolumeContribution

                # Only count lesion if it has at least one valid slice
                if lesionHasValidSlice:
                    validLesionCount += 1
                    totalAgatston += lesionScore
                    totalVolume += lesionVolume

            # Calculate statistics
            meanDensity = np.mean(allDensities) if allDensities else 0
            maxDensity = np.max(allDensities) if allDensities else 0

            # Calcium Mass calculation using calibration factor
            # Formula: Mass (mg) = Volume (mm³) × Mean_Density (HU) × Calibration_Factor / 1000
            # Default calibration factor: 0.81 (for 130 HU = 114.5 mg/cm³ CaHA)
            equivalentMass = totalVolume * meanDensity * calibrationFactor / 1000 if totalVolume > 0 else 0

            print(f"  {territory}: Score={totalAgatston:.1f} AU, Lesions={validLesionCount}, Volume={totalVolume:.1f} mm³, Mass={equivalentMass:.1f} mg")

            return {
                'agatston_score': totalAgatston,
                'total_volume_mm3': totalVolume,
                'equivalent_mass_mg': equivalentMass,
                'num_lesions': validLesionCount,
                'mean_density': meanDensity,
                'max_density': maxDensity
            }

        except Exception as e:
            print(f"Error calculating score for {territory}: {str(e)}")
            import traceback
            traceback.print_exc()
            return self.getEmptyResult()

    def getEmptyResult(self):
        """Return empty result dictionary"""
        return {
            'agatston_score': 0,
            'total_volume_mm3': 0,
            'equivalent_mass_mg': 0,
            'num_lesions': 0,
            'mean_density': 0,
            'max_density': 0
        }

    def getMESACategory(self, totalScore):
        """Get MESA/ACC risk category based on total score"""
        for cat in self.MESA_CATEGORIES:
            if cat['min'] <= totalScore <= cat['max']:
                return {
                    'category': cat['category'],
                    'description': cat['description'],
                    'risk': cat['risk']
                }
        return {'category': 'Unknown', 'description': '', 'risk': ''}

    def getMESAPercentile(self, score, patientInfo):
        """Calculate MESA percentile based on age, sex, and ethnicity"""
        age = patientInfo.get('age')
        sex = 'Male' if patientInfo.get('sex') == 'M' else 'Female'
        ethnicity = patientInfo.get('ethnicity', 'Caucasian')

        if not age or age < 45:
            return {'percentile': '--', 'comparison': 'Age 45+ required for percentile calculation'}

        # Determine age group
        if age < 55:
            ageGroup = '45-54'
        elif age < 65:
            ageGroup = '55-64'
        elif age < 75:
            ageGroup = '65-74'
        else:
            ageGroup = '75-84'

        try:
            percentileData = self.MESA_PERCENTILES[sex][ethnicity][ageGroup]

            # Interpolate percentile
            if score == 0:
                percentile = 25 if percentileData['25'] == 0 else 10
            elif score <= percentileData['25']:
                percentile = int(25 * score / max(percentileData['25'], 1))
            elif score <= percentileData['50']:
                percentile = 25 + int(25 * (score - percentileData['25']) / max(percentileData['50'] - percentileData['25'], 1))
            elif score <= percentileData['75']:
                percentile = 50 + int(25 * (score - percentileData['50']) / max(percentileData['75'] - percentileData['50'], 1))
            elif score <= percentileData['90']:
                percentile = 75 + int(15 * (score - percentileData['75']) / max(percentileData['90'] - percentileData['75'], 1))
            else:
                percentile = 90 + int(10 * min((score - percentileData['90']) / percentileData['90'], 1))

            percentile = min(99, max(1, percentile))

            # Interpretation
            if percentile <= 25:
                comparison = f"Below average for {sex.lower()}, {ageGroup} years, {ethnicity}"
            elif percentile <= 50:
                comparison = f"Average for {sex.lower()}, {ageGroup} years, {ethnicity}"
            elif percentile <= 75:
                comparison = f"Above average for {sex.lower()}, {ageGroup} years, {ethnicity}"
            else:
                comparison = f"High for {sex.lower()}, {ageGroup} years, {ethnicity}"

            return {
                'percentile': percentile,
                'ageGroup': ageGroup,
                'comparison': comparison
            }

        except KeyError:
            return {'percentile': '--', 'comparison': 'Percentile data not available'}

    def create3DVisualization(self, segmentationsByTerritory, territories, volumeNode=None):
        """Create 3D visualization with territory-based coloring"""
        # Set up 3D-only layout
        layoutManager = slicer.app.layoutManager()
        layoutManager.setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutOneUp3DView)

        threeDWidget = layoutManager.threeDWidget(0)
        threeDView = threeDWidget.threeDView()

        # Show each territory's segmentation in 3D with full opacity
        for territory, segNode in segmentationsByTerritory.items():
            if segNode:
                segNode.CreateClosedSurfaceRepresentation()
                displayNode = segNode.GetDisplayNode()
                if displayNode:
                    displayNode.SetVisibility3D(True)
                    displayNode.SetOpacity3D(1.0)  # Full opacity for calcium

        # Reset 3D view and set background
        threeDView.resetFocalPoint()
        threeDView.resetCamera()

        # Configure 3D view appearance
        viewNode = threeDView.mrmlViewNode()
        viewNode.SetBackgroundColor(0.1, 0.1, 0.15)
        viewNode.SetBackgroundColor2(0.2, 0.2, 0.25)

        # Hide orientation axes and box
        viewNode.SetOrientationMarkerType(slicer.vtkMRMLViewNode.OrientationMarkerTypeNone)
        viewNode.SetBoxVisible(False)
        viewNode.SetAxisLabelsVisible(False)

    def showMPRView(self, volumeNode=None):
        """Restore standard multiplanar reconstruction (MPR) view"""
        layoutManager = slicer.app.layoutManager()
        layoutManager.setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutFourUpView)

        # Center on volume if provided
        if volumeNode:
            slicer.util.setSliceViewerLayers(background=volumeNode)
            slicer.util.resetSliceViews()

    def createCharts(self, results, territories, patientInfo=None, showInWindow=True):
        """Create analysis charts and optionally display in a window

        Returns the path to the saved chart image for PDF inclusion.
        """
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            import numpy as np

            # Check if we have valid patient info for percentile chart
            hasPercentileData = False
            percentileData = None
            patientScore = results.get('Total', {}).get('agatston_score', 0)

            if patientInfo:
                age = patientInfo.get('age')
                sex = 'Male' if patientInfo.get('sex') == 'M' else 'Female'
                ethnicity = patientInfo.get('ethnicity', 'Caucasian')

                if age and age >= 45:
                    # Determine age group
                    if age < 55:
                        ageGroup = '45-54'
                    elif age < 65:
                        ageGroup = '55-64'
                    elif age < 75:
                        ageGroup = '65-74'
                    else:
                        ageGroup = '75-84'

                    try:
                        percentileData = self.MESA_PERCENTILES[sex][ethnicity][ageGroup]
                        hasPercentileData = True
                    except KeyError:
                        pass

            # Create figure with 2 or 3 subplots depending on percentile data availability
            if hasPercentileData:
                fig, axes = plt.subplots(1, 3, figsize=(14, 4))
            else:
                fig, axes = plt.subplots(1, 2, figsize=(10, 4))

            # Bar chart of scores by territory
            territoryNames = list(territories.keys())
            scores = [results.get(t, {}).get('agatston_score', 0) for t in territoryNames]
            colors = [territories[t]['colorHex'] for t in territoryNames]

            axes[0].bar(territoryNames, scores, color=colors, edgecolor='black', linewidth=0.5)
            axes[0].set_ylabel('Agatston Score (AU)')
            axes[0].set_title('Calcium Score by Territory')
            axes[0].set_xlabel('Coronary Territory')

            # Add value labels on bars
            for i, (name, score) in enumerate(zip(territoryNames, scores)):
                if score > 0:
                    axes[0].text(i, score + max(scores)*0.02, f'{score:.0f}', ha='center', va='bottom', fontsize=9)

            # Pie chart of distribution
            nonzeroScores = [(t, s) for t, s in zip(territoryNames, scores) if s > 0]
            if nonzeroScores:
                labels, values = zip(*nonzeroScores)
                pieColors = [territories[t]['colorHex'] for t in labels]
                wedges, texts, autotexts = axes[1].pie(values, labels=labels, colors=pieColors,
                                                        autopct='%1.1f%%', startangle=90)
                axes[1].set_title('Distribution by Territory')
            else:
                axes[1].text(0.5, 0.5, 'No calcium detected', ha='center', va='center', fontsize=12)
                axes[1].set_title('Distribution by Territory')
                axes[1].axis('off')

            # Percentile chart (if data available)
            if hasPercentileData:
                ax = axes[2]

                # Get percentile values
                p25 = percentileData['25']
                p50 = percentileData['50']
                p75 = percentileData['75']
                p90 = percentileData['90']

                # Create bar chart for percentile reference values
                percentileLabels = ['25th', '50th', '75th', '90th']
                percentileValues = [p25, p50, p75, p90]
                barColors = ['#4CAF50', '#2196F3', '#FF9800', '#F44336']  # Green, Blue, Orange, Red

                bars = ax.bar(percentileLabels, percentileValues, color=barColors, edgecolor='black', linewidth=0.5, alpha=0.7)

                # Add horizontal line for patient's score
                ax.axhline(y=patientScore, color='#9C27B0', linewidth=2.5, linestyle='--', label=f'Your Score: {patientScore:.0f} AU')

                # Calculate and display patient's percentile
                mesaResult = self.getMESAPercentile(patientScore, patientInfo)
                patientPercentile = mesaResult.get('percentile', '--')

                ax.set_ylabel('Agatston Score (AU)')
                ax.set_xlabel('Population Percentile')
                ax.set_title(f'MESA Percentile\n({sex}, {ageGroup} yrs, {ethnicity})')
                ax.legend(loc='upper left', fontsize=9)

                # Add value labels on bars
                for bar, val in zip(bars, percentileValues):
                    if val > 0:
                        ax.text(bar.get_x() + bar.get_width()/2, val + max(percentileValues)*0.02,
                               f'{val:.0f}', ha='center', va='bottom', fontsize=9)

                # Add patient percentile annotation
                if patientPercentile != '--':
                    ax.annotate(f'Your Percentile: {patientPercentile}th',
                               xy=(0.5, 0.95), xycoords='axes fraction',
                               ha='center', fontsize=10, fontweight='bold',
                               bbox=dict(boxstyle='round', facecolor='#E1BEE7', alpha=0.8))

            plt.tight_layout()

            # Save chart
            import tempfile
            chartPath = os.path.join(tempfile.gettempdir(), 'coronary_cac_charts.png')
            plt.savefig(chartPath, dpi=150, bbox_inches='tight', facecolor='white')
            plt.close()

            # Display in a Qt window if requested
            if showInWindow:
                self.showChartWindow(chartPath, results, territories, patientInfo)

            return chartPath

        except ImportError:
            slicer.util.warningDisplay("matplotlib not installed. Please install dependencies.")
            return None

    def showChartWindow(self, chartPath, results, territories, patientInfo=None):
        """Display charts in a popup window"""
        # Create dialog
        dialog = qt.QDialog()
        dialog.setWindowTitle("Coronary Calcium Score - Charts")
        dialog.setMinimumSize(900, 550)

        layout = qt.QVBoxLayout(dialog)

        # Add title
        titleLabel = qt.QLabel("<h2>Calcium Score Distribution & MESA Percentile</h2>")
        titleLabel.setAlignment(qt.Qt.AlignCenter)
        layout.addWidget(titleLabel)

        # Add chart image
        if os.path.exists(chartPath):
            imageLabel = qt.QLabel()
            pixmap = qt.QPixmap(chartPath)
            scaledPixmap = pixmap.scaledToWidth(850, qt.Qt.SmoothTransformation)
            imageLabel.setPixmap(scaledPixmap)
            imageLabel.setAlignment(qt.Qt.AlignCenter)
            layout.addWidget(imageLabel)

        # Add summary
        totalScore = results.get('Total', {}).get('agatston_score', 0)
        summaryText = f"<b>Total Agatston Score: {totalScore:.1f} AU</b>"

        # Add percentile info if available
        if patientInfo:
            mesaResult = self.getMESAPercentile(totalScore, patientInfo)
            percentile = mesaResult.get('percentile', '--')
            if percentile != '--':
                summaryText += f"  |  <b>MESA Percentile: {percentile}th</b>"

        summaryLabel = qt.QLabel(summaryText)
        summaryLabel.setAlignment(qt.Qt.AlignCenter)
        layout.addWidget(summaryLabel)

        # Close button
        closeButton = qt.QPushButton("Close")
        closeButton.clicked.connect(lambda: dialog.close())
        layout.addWidget(closeButton)

        dialog.exec_()

    def generatePDFReport(self, outputDir, results, patientInfo, territories, segmentationsByTerritory, volumeNode, options):
        """Generate comprehensive PDF report"""
        try:
            from reportlab.lib.pagesizes import letter
            from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
            from reportlab.lib.units import inch
            from reportlab.lib import colors
            from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image
            from reportlab.lib.enums import TA_CENTER, TA_LEFT

            # Create filename
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            filename = f"CoronaryCAC_Report_{timestamp}.pdf"
            filepath = os.path.join(outputDir, filename)

            # Create document
            doc = SimpleDocTemplate(filepath, pagesize=letter, rightMargin=72, leftMargin=72, topMargin=72, bottomMargin=72)

            # Styles
            styles = getSampleStyleSheet()
            titleStyle = ParagraphStyle('Title', parent=styles['Heading1'], fontSize=18, alignment=TA_CENTER)

            # Build content
            story = []

            # Logo and header
            if options.get('logoPath') and os.path.exists(options['logoPath']):
                logo = Image(options['logoPath'], width=1.5*inch, height=1.5*inch)
                story.append(logo)

            if options.get('companyDesc'):
                story.append(Paragraph(options['companyDesc'], styles['Normal']))

            story.append(Spacer(1, 0.3*inch))

            # Title
            story.append(Paragraph("CORONARY ARTERY - CALCIUM SCORE", titleStyle))
            story.append(Spacer(1, 0.3*inch))

            # Patient Information
            story.append(Paragraph("<b>PATIENT INFORMATION</b>", styles['Heading2']))

            patientData = [
                ['Name', patientInfo.get('name', 'N/A')],
                ['ID', patientInfo.get('id', 'N/A')],
                ['Sex', 'Male' if patientInfo.get('sex') == 'M' else 'Female'],
                ['Age', str(patientInfo.get('age', 'N/A'))],
                ['Ethnicity', patientInfo.get('ethnicity', 'N/A')],
                ['Analysis Date', patientInfo.get('date', 'N/A')]
            ]

            patientTable = Table(patientData, colWidths=[2*inch, 4*inch])
            patientTable.setStyle(TableStyle([
                ('BACKGROUND', (0, 0), (0, -1), colors.lightgrey),
                ('TEXTCOLOR', (0, 0), (-1, -1), colors.black),
                ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
                ('FONTNAME', (0, 0), (0, -1), 'Helvetica-Bold'),
                ('FONTSIZE', (0, 0), (-1, -1), 10),
                ('BOTTOMPADDING', (0, 0), (-1, -1), 6),
                ('TOPPADDING', (0, 0), (-1, -1), 6),
                ('GRID', (0, 0), (-1, -1), 0.5, colors.grey)
            ]))
            story.append(patientTable)
            story.append(Spacer(1, 0.3*inch))

            # Analysis Parameters
            story.append(Paragraph("<b>ANALYSIS PARAMETERS</b>", styles['Heading2']))
            threshold = options.get('threshold', 130)
            calibFactor = options.get('calibrationFactor', 0.81)
            paramData = [
                ['Threshold', f'{threshold} HU (Standard Agatston)'],
                ['Mass Calibration Factor', f'{calibFactor:.2f}'],
                ['Method', '2D slice-by-slice connected components'],
                ['Minimum Lesion Area', '1.0 mm² per slice']
            ]
            paramTable = Table(paramData, colWidths=[2*inch, 4*inch])
            paramTable.setStyle(TableStyle([
                ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
                ('FONTNAME', (0, 0), (0, -1), 'Helvetica-Bold'),
                ('FONTSIZE', (0, 0), (-1, -1), 9),
                ('BOTTOMPADDING', (0, 0), (-1, -1), 4),
                ('TOPPADDING', (0, 0), (-1, -1), 4),
                ('GRID', (0, 0), (-1, -1), 0.5, colors.grey)
            ]))
            story.append(paramTable)
            story.append(Spacer(1, 0.3*inch))

            # Results by Territory
            story.append(Paragraph("<b>AGATSTON SCORE BY TERRITORY</b>", styles['Heading2']))

            resultsData = [['Territory', 'Score (AU)', 'Volume (mm³)', 'Mass (mg)', 'Lesions']]

            for territory in list(territories.keys()) + ['Total']:
                result = results.get(territory, {})
                row = [
                    territory if territory != 'Total' else 'TOTAL',
                    f"{result.get('agatston_score', 0):.1f}",
                    f"{result.get('total_volume_mm3', 0):.1f}",
                    f"{result.get('equivalent_mass_mg', 0):.1f}",
                    str(result.get('num_lesions', 0))
                ]
                resultsData.append(row)

            resultsTable = Table(resultsData, colWidths=[1.2*inch, 1.2*inch, 1.2*inch, 1.2*inch, 1*inch])
            resultsTable.setStyle(TableStyle([
                ('BACKGROUND', (0, 0), (-1, 0), colors.darkblue),
                ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
                ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
                ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                ('FONTSIZE', (0, 0), (-1, -1), 10),
                ('BOTTOMPADDING', (0, 0), (-1, -1), 8),
                ('TOPPADDING', (0, 0), (-1, -1), 8),
                ('GRID', (0, 0), (-1, -1), 0.5, colors.grey),
                ('BACKGROUND', (0, -1), (-1, -1), colors.lightgrey),
                ('FONTNAME', (0, -1), (-1, -1), 'Helvetica-Bold')
            ]))
            story.append(resultsTable)
            story.append(Spacer(1, 0.3*inch))

            # Risk Classification
            story.append(Paragraph("<b>RISK CLASSIFICATION</b>", styles['Heading2']))

            totalResult = results.get('Total', {})
            category = totalResult.get('mesa_category', {})
            percentile = totalResult.get('mesa_percentile', {})

            riskData = [
                ['MESA/ACC Category', category.get('category', 'N/A')],
                ['Score Range', category.get('description', 'N/A')],
                ['Risk Level', category.get('risk', 'N/A')],
                ['MESA Percentile', f"{percentile.get('percentile', 'N/A')}th"],
                ['Interpretation', percentile.get('comparison', 'N/A')]
            ]

            riskTable = Table(riskData, colWidths=[2*inch, 4*inch])
            riskTable.setStyle(TableStyle([
                ('BACKGROUND', (0, 0), (0, -1), colors.lightgrey),
                ('TEXTCOLOR', (0, 0), (-1, -1), colors.black),
                ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
                ('FONTNAME', (0, 0), (0, -1), 'Helvetica-Bold'),
                ('FONTSIZE', (0, 0), (-1, -1), 10),
                ('BOTTOMPADDING', (0, 0), (-1, -1), 6),
                ('TOPPADDING', (0, 0), (-1, -1), 6),
                ('GRID', (0, 0), (-1, -1), 0.5, colors.grey)
            ]))
            story.append(riskTable)
            story.append(Spacer(1, 0.3*inch))

            # Add charts if option enabled
            if options.get('includeCharts', True):
                # Generate chart image (include percentile chart if option enabled)
                chartPath = self.createCharts(results, territories, patientInfo if options.get('includePercentileChart', True) else None, showInWindow=False)
                if chartPath and os.path.exists(chartPath):
                    story.append(Paragraph("<b>CALCIUM DISTRIBUTION CHARTS</b>", styles['Heading2']))
                    from reportlab.platypus import Image
                    # Wider image to accommodate percentile chart (7 inches for 3 charts, 6 for 2)
                    chartWidth = 7*inch if options.get('includePercentileChart', True) else 6*inch
                    chartImg = Image(chartPath, width=chartWidth, height=2.4*inch)
                    story.append(chartImg)
                    story.append(Spacer(1, 0.3*inch))

            # References
            story.append(Paragraph("<b>REFERENCES</b>", styles['Heading2']))
            references = [
                "1. Agatston AS, et al. J Am Coll Cardiol. 1990;15(4):827-832.",
                "2. McClelland RL, et al. Circulation. 2006;113(1):30-37. (MESA Study)",
                "3. ACC/AHA Guidelines for Cardiovascular Risk Assessment. 2019.",
                "4. Greenland P, et al. J Am Coll Cardiol. 2018;71(22):e267-e346."
            ]
            for ref in references:
                story.append(Paragraph(ref, styles['Normal']))

            story.append(Spacer(1, 0.3*inch))

            # Disclaimer
            disclaimer = """
            <i><b>Clinical Disclaimer:</b> This report is intended for research and educational purposes only.
            The calcium scores and risk assessments should be interpreted by a qualified healthcare professional
            in the context of the patient's complete clinical picture. This software is not a substitute for
            professional medical advice, diagnosis, or treatment.</i>
            """
            story.append(Paragraph(disclaimer, styles['Normal']))

            # Footer
            story.append(Spacer(1, 0.2*inch))
            footer = f"Generated by CoronaryCaScore v1.0 | AITeRTC | {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"
            story.append(Paragraph(footer, ParagraphStyle('Footer', parent=styles['Normal'], fontSize=8, alignment=TA_CENTER)))

            # Build PDF
            doc.build(story)

            return filepath

        except ImportError:
            raise Exception("reportlab not installed. Please install dependencies.")


#
# CoronaryCaScoreTest
#

class CoronaryCaScoreTest(ScriptedLoadableModuleTest):
    """Test suite for Coronary Calcium Score module"""

    def setUp(self):
        slicer.mrmlScene.Clear()

    def runTest(self):
        self.setUp()
        self.test_Basic()

    def test_Basic(self):
        """Basic test - check that module loads"""
        self.delayDisplay("Starting test")

        logic = CoronaryCaScoreLogic()

        # Test MESA category classification
        category = logic.getMESACategory(0)
        self.assertEqual(category['category'], 'Zero')

        category = logic.getMESACategory(50)
        self.assertEqual(category['category'], 'Mild')

        category = logic.getMESACategory(250)
        self.assertEqual(category['category'], 'Moderate')

        category = logic.getMESACategory(500)
        self.assertEqual(category['category'], 'Severe')

        self.delayDisplay('Test passed!')
