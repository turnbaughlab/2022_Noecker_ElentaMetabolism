This document describes the full set of parameters used for processing of stable isotope-resolved metabolomics data in Noecker et al (2022), using Thermo Scientific Compound Discoverer software v3.3.

Input Data: File Name(s) (Hidden)

# Processing node: Select Spectra
1. Spectrum Properties Filter:
- Lower RT Limit: 0
- Upper RT Limit: 0
- First Scan: 0
- Last Scan: 0
- Ignore Specified Scans: (not specified)
- Lowest Charge State: 0
- Highest Charge State: 0
- Min. Precursor Mass: 100 Da
- Max. Precursor Mass: 5000 Da
- Total Intensity Threshold: 0
- Minimum Peak Count: 1
2. Scan Event Filters:
- Mass Analyzer: (not specified)
- MS Order: Any
- Activation Type: (not specified)
- Min. Collision Energy: 0
- Max. Collision Energy: 1000
- Scan Type: Any
- Polarity Mode: (not specified)
- MS1 Mass Range: (not specified)
- FAIMS CV: (not specified)
3. Peak Filters:
- S/N Threshold (FT-only): 1.5
4. Replacements for Unrecognized Properties:
- Unrecognized Charge Replacements: 1
- Unrecognized Mass Analyzer Replacements: ITMS
- Unrecognized MS Order Replacements: MS2
- Unrecognized Activation Type Replacements: CID
- Unrecognized Polarity Replacements: +
- Unrecognized MS Resolution@200 Replacements: 60000
- Unrecognized MSn Resolution@200 Replacements: 30000
5. General Settings:
- Precursor Selection: Use MS(n - 1) Precursor
- Use Isotope Pattern in Precursor Reevaluation: True
- Provide Profile Spectra: Automatic
- Store Chromatograms: False

# Processing node : Align Retention Times (ChromAlign)
1. General Settings:
- Reference File: PETU003_R_2243_TP3_1_Pos_QE2_Hilic_106 or PETU003_R_2243_TP3_1_Neg_QE2_Hilic_106 for positive and negative modes respectively.

# Processing node: Detect Compounds (Legacy)
1. General Settings:
- Intensity Tolerance [%]: 30
- Ions: [M+H]+1 or [M-H]-1 for positive and negative modes respectively.
- Base Ions: [M+H]+1; [M-H]-1
- Min. Element Counts: C H
- Max. Element Counts: C90 H190 Br3 Cl4 K2 N10 Na2 O18 P3 S5
- Mass Tolerance [ppm]: 5 ppm
- Min. Peak Intensity: 10000
- Min. # Scans per Peak: 5
- Use Most Intense Isotope Only: False
2. Trace Detection:
- Max. Number of Gaps to Correct: 2
- Min. Number of Adjacent Non-Zeros: 2
3. Peak Detection:
- Filter Peaks: True
- Max. Peak Width [min]: 0.5
- Remove Singlets: True
4. Isotope Pattern Detection:
- Min. # Isotopes: 2
- Use Peak Quality for Isotope Grouping: False
- Filter out Features with Bad Peaks Only: True
- Zig-Zag Index Threshold: 0.2
- Jaggedness Threshold: 0.4
- Modality Threshold: 0.9
- Min. Spectral Distance Score: 0
- Remove Potentially False Positive Isotopes: True
6. AcquireX Settings:
- Detect Persistent Background Ions: False

# Processing node: Group Compounds
1. General Settings:
- Mass Tolerance: 5 ppm
- RT Tolerance [min]: 0.2
- Align Peaks: False
- Preferred Ions: [M+H]+1; [M-H]-1
- Area Integration: All Ions
2. Peak Rating Contributions:
- Area Contribution: 3
- CV Contribution: 10
- FWHM to Base Contribution: 5
- Jaggedness Contribution: 5
- Modality Contribution: 5
- Zig-Zag Index Contribution: 5
3. Peak Rating Filter:
- Peak Rating Threshold: 4
- Number of Files: 3

# Processing node: Assign Compound Annotations
1. General Settings:
- Mass Tolerance: 5 ppm
2. Data Sources:
- Data Source #1: mzCloud Search
- Data Source #2: MassList Search
- Data Source #3: Predicted Compositions
- Data Source #4: (not specified)
- Data Source #5: (not specified)
- Data Source #6: (not specified)
- Data Source #7: (not specified)
3. Scoring Rules:
- Use mzLogic: True
- Use Spectral Distance: True
- SFit Threshold: 20
- SFit Range: 20
4. Reprocessing:
- Clear Names: False

# Processing node: Analyze Labeled Compounds
1. Label Settings:
- Label Element: [13]C
- Max. Exchange: 25
- Source Efficiency [%]: 100
2. Pattern Analysis:
- Mass Tolerance [ppm]: 5 ppm
- Intensity Tolerance [%]: 30
- Intensity Threshold [%]: 2
- S/N Threshold: 5
3. General Settings:
- Mark Irregular Exchange: True
- Exclude Blanks: True
- Hide Unprocessed: True

# Processing node 40: Search mzCloud
1. General Settings:
- Compound Classes: All
- Precursor Mass Tolerance: 10 ppm- FT Fragment Mass Tolerance: 10 ppm
- IT Fragment Mass Tolerance: 0.4 Da
- Library: Autoprocessed; Reference
- Post Processing: Recalibrated
- Max. # Results: 10
- Annotate Matching Fragments: True
- Search MSn Tree: False
2. DDA Search:
- Identity Search: Cosine
- Match Activation Type: True
- Match Activation Energy: Match with Tolerance
- Activation Energy Tolerance: 20
- Apply Intensity Threshold: True
- Similarity Search: None
- Match Factor Threshold: 85
3. DIA Search:
- Use DIA Scans for Search: False
- Max. Isolation Width [Da]: 500
- Match Activation Type: False
- Match Activation Energy: Any
- Activation Energy Tolerance: 100
- Apply Intensity Threshold: False
- Match Factor Threshold: 20
 
# Processing node: Search Mass Lists
1. Search Settings:
- Mass Lists: Combined Hilic Mass mzRT library.massList
- Mass Tolerance: 5 ppm
- Use Retention Time: True
- RT Tolerance [min]: 0.3

# Processing node: Predict Compositions
1. Prediction Settings:
- Mass Tolerance: 5 ppm
- Min. Element Counts: C H
- Max. Element Counts: C90 H190 Br3 Cl4 N10 O18 P3 S5
- Min. RDBE: 0
- Max. RDBE: 40
- Min. H/C: 0.1
- Max. H/C: 4
- Max. # Candidates: 10
- Max. # Internal Candidates: 200
2. Pattern Matching:
- Intensity Tolerance [%]: 30
- Intensity Threshold [%]: 0.1
- S/N Threshold: 3
- Min. Spectral Fit [%]: 30
- Min. Pattern Cov. [%]: 90
- Use Dynamic Recalibration: True
3. Fragments Matching:
- Use Fragments Matching: True
- Mass Tolerance: 5 ppm
- S/N Threshold: 3

# Processing node: Mark Background Compounds
1. General Settings:
- Max. Sample/Blank: 3
- Max. Blank/Sample: 0
- Hide Background: True

# Processing node: Descriptive Statistics
No parameters

# Processing node: Differential Analysis
1. General Settings:
- Log10 Transform Values: True
2. Peak Rating Contributions:
- Update Peak Rating: True
- Area Contribution: 3
- CV Contribution: 10
- FWHM to Base Contribution: 5
- Jaggedness Contribution: 5
- Modality Contribution: 5
- Zig-Zag Index Contribution: 5
