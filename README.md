# CoronaryCaScore - Coronary Artery Calcium Scoring Plugin

3D Slicer extension for automated quantification of coronary artery calcification using the Agatston scoring method with per-vessel analysis.

![Version](https://img.shields.io/badge/version-1.0-blue)
![License](https://img.shields.io/badge/license-Restrictive-orange)
![3D Slicer](https://img.shields.io/badge/3D%20Slicer-5.0+-red)

## Video Tutorial

[![Video Tutorial](https://img.youtube.com/vi/Vl_HmS9TFlQ/0.jpg)](https://youtu.be/Vl_HmS9TFlQ)

**[Watch the installation and usage tutorial on YouTube](https://youtu.be/Vl_HmS9TFlQ)**

## Authors

**Vittorio Censullo**
**AITeRTC** - *Associazione Italiana Tecnici di Radiologia Esperti in TC*
2025

## Key Features

### Core Functionality
- **Per-Vessel Analysis**: Separate scoring for LAD, LCx, RCA, and Left Main
- **Guided Tab-Based Workflow**: Intuitive step-by-step process for coronary calcium scoring
- **Advanced Segmentation Methods**:
  - **Point and Click**: One-click calcium island segmentation (Region Growing)
  - **ROI Box**: Rectangular box threshold-based segmentation
  - **Paint Mode**: Manual painting with threshold-guided brush
- **Color-Coded Territories**: Visual distinction between coronary vessels
  - LAD: Red
  - LCx: Blue
  - RCA: Green
  - Left Main: Yellow

### Risk Classification
- **MESA/ACC Categories**:
  - Zero (0 AU): Very Low Risk
  - Minimal (1-10 AU): Low Risk
  - Mild (11-100 AU): Low-Moderate Risk
  - Moderate (101-400 AU): Moderate-High Risk
  - Severe (>400 AU): High Risk
- **MESA Percentiles**: Age, sex, and ethnicity-adjusted percentile calculation

### Settings & Persistence
- **Configurable Threshold**: Default 130 HU for non-contrast CT
- **Persistent Settings**: Logo, company description, and threshold saved across sessions
- **Easy Dependency Management**: One-click installation of required libraries

## Installation

### Method 1: Download ZIP (Easiest)

1. Download the repository as ZIP from [GitHub](https://github.com/vcensullo/CoronaryCaScore/archive/refs/heads/main.zip)

2. Extract the ZIP to a folder of your choice

3. Open 3D Slicer

4. Go to **Edit -> Application Settings -> Modules**

5. In "Additional module paths", click **Add** and select the `CoronaryCaScore-main/CoronaryCaScore` folder

6. Restart 3D Slicer

7. The module will appear under **Modules -> Cardiac -> Coronary Artery Calcium Score**

### Method 2: Git Clone

1. Clone this repository:
   ```bash
   git clone https://github.com/vcensullo/CoronaryCaScore.git
   ```

2. Open 3D Slicer

3. Go to **Edit -> Application Settings -> Modules**

4. In "Additional module paths", click **Add** and select the `CoronaryCaScore/CoronaryCaScore` folder

5. Restart 3D Slicer

6. The module will appear under **Modules -> Cardiac -> Coronary Artery Calcium Score**

## Dependencies

The plugin requires the `reportlab` and `matplotlib` libraries for PDF and chart generation.

1. Open the plugin in 3D Slicer
2. Go to **Tab 0 - Settings**
3. Click **Install/Update Dependencies**

## Usage Workflow

### Patient Information
- Open the accordion menu to insert patient data
- Name, ID, and age are optional
- **Sex and ethnicity are required** for accurate MESA percentile calculation

### Tab 0: Settings
- Install required dependencies (reportlab/matplotlib)
- Configure threshold (default: 130 HU for non-contrast CT)
- Set company logo and description for PDF reports
- View plugin information and version

### Tab 1: Setup View
1. Select CT volume (non-contrast, ECG-gated preferred)
2. Apply optimal layout for cardiac visualization
3. Tab 2 will enable automatically after volume selection

### Tab 2: Calcium Segmentation
1. **Select coronary territory** (LAD, LCx, RCA, or LM)
2. Choose segmentation method:
   - **Click & Grow**: Click on calcium to automatically segment
   - **ROI Box**: Place box around coronary artery -> Apply threshold
   - **Paint Mode**: Manually paint calcium regions
3. Repeat for each territory with calcium
4. Territory summary shows segmentation status
5. Tab 3 enables after segmentation

### Tab 3: Results & Analysis
1. Click **Calculate All Scores**
2. View results table with per-territory breakdown:
   - Agatston Score (AU)
   - Volume (mm³)
   - Equivalent Mass (mg)
   - Number of lesions
3. View **MESA/ACC risk category**
4. View **MESA percentile** (age/sex/ethnicity adjusted)
5. **Show 3D**: Visualize calcium with territory-based coloring
6. **Show Charts**: View distribution charts
7. Tab 4 enables after calculation

### Tab 4: Generate Report
1. Select output directory
2. Configure report options (screenshots, 3D view, charts)
3. Click **Generate PDF Report**
4. Report includes per-territory breakdown, risk classification, and MESA percentile

## Technical Details

### Agatston Score Calculation

```
Agatston Score = Sum (Lesion Area x Density Factor)
```

**Density Factors:** 1 (130-199 HU), 2 (200-299 HU), 3 (300-399 HU), 4 (>=400 HU)

**Minimum Lesion Size:** 1.0 mm² per slice (Agatston standard)

### MESA Percentiles

Based on the Multi-Ethnic Study of Atherosclerosis (MESA), percentiles are calculated considering:
- **Age groups**: 45-54, 55-64, 65-74, 75-84 years
- **Sex**: Male, Female
- **Ethnicity**: Caucasian, African-American, Hispanic, Chinese

### Territory Color Coding

| Territory | Full Name | Color |
|-----------|-----------|-------|
| LAD | Left Anterior Descending | Red |
| LCx | Left Circumflex | Blue |
| RCA | Right Coronary Artery | Green |
| LM | Left Main | Yellow |

## References

1. **Agatston AS, et al.** Quantification of coronary artery calcium using ultrafast computed tomography. *J Am Coll Cardiol.* 1990;15(4):827-832.

2. **McClelland RL, et al.** Distribution of coronary artery calcium by race, gender, and age: results from the Multi-Ethnic Study of Atherosclerosis (MESA). *Circulation.* 2006;113(1):30-37.

3. **Greenland P, et al.** 2018 AHA/ACC/AACVPR/AAPA/ABC/ACPM/ADA/AGS/APhA/ASPC/NLA/PCNA Guideline on the Management of Blood Cholesterol. *J Am Coll Cardiol.* 2019;73(24):e285-e350.

4. **Budoff MJ, et al.** Expert consensus document on coronary artery calcium scoring. *J Cardiovasc Comput Tomogr.* 2023;17(3):151-163.

## Contributing

Contributions welcome! Fork, create feature branch, commit, push, and open a Pull Request.

## License

**Restrictive License** - Free for personal, educational, research, and clinical use. Redistribution and commercial use prohibited without authorization. See [LICENSE](LICENSE) for details.

## Acknowledgments

- **Developed by**: Vittorio Censullo
- **Co-authored by**: **AITeRTC** - *Associazione Italiana Tecnici di Radiologia Esperti in TC*
- Built using **3D Slicer** open-source platform
- Based on **ACC/AHA** clinical guidelines and **MESA Study** data

## Contact

**Vittorio Censullo**
**AITeRTC** - Associazione Italiana Tecnici di Radiologia Esperti in TC
GitHub: [@vcensullo](https://github.com/vcensullo)

Issues: [GitHub Issues](https://github.com/vcensullo/CoronaryCaScore/issues)

## Citation

If you use this plugin in your research, please cite:

```bibtex
@software{censullo2025coronarycascore,
  author = {Censullo, Vittorio and AITeRTC},
  title = {CoronaryCaScore: Coronary Artery Calcium Scoring Plugin for 3D Slicer},
  year = {2025},
  publisher = {GitHub},
  url = {https://github.com/vcensullo/CoronaryCaScore},
  institution = {AITeRTC - Associazione Italiana Tecnici di Radiologia Esperti in TC}
}
```

## Version History

**Version 1.0** (December 2025) - Initial Release
- Per-vessel calcium scoring (LAD, LCx, RCA, LM)
- MESA/ACC risk classification
- MESA percentile calculation
- Color-coded territory visualization
- Comprehensive PDF reporting

---

**Version:** 1.0
**Last Updated:** December 2025
**Maintained by:** Vittorio Censullo & AITeRTC

*Vibe-coded with the aid of Claude*
