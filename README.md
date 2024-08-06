# Kendrick Mass Defect Analysis for High-Resolution mass spectrometry data
-----------------------------------------------------------------------------------
Description
-----------
Perfomed Kendrick Mass Defect (KMD) analysis from a feature list previously extracted from a peak picking algorithm (e.g., XCMS)
The workflow combine existing code (step) (see below in citing section) that might be run independently.

The step are the following:
1-PatRoon_XCMS_feature_opti: optimized the parameter for XCMS peak picking
2-PatRoon_featuregroup: 
3-KMD_analysis: 
creat_KMD_susp_list: create a specific KMD suspect list to match the KMD result with some suspect and help the identification of unknow (DO NOT CONSIDR THIS AS A TRUE SUSPECT SCREENING STEP


Dependency
----------------
R(>=4.3.1)
patRoon(>=2.3.0)

Getting started
----------------
1.
2.


Citing
------

The current workflow combine/translate divers code from the following sources

Peak picking with XCMS: 
Helmus, R.; van de Velde, B.; Brunner, A. M.; ter Laak, T. L.; van Wezel, A. P.; Schymanski, E. L. patRoon 2.0: Improved non-target analysis workflows including automated transformation product screening. Journal of Open Source Software 2022, 7 (71), 4029. [DOI: 10.21105/joss.04029](https://doi.org/10.21105/joss.04029); [github patRoon](https://github.com/rickhelmus/patRoon).

KMD analysis:
Zweigle, J.; Bugsel, B.; Fabregat-Palau, J.; Zwiener, C. PFDeltaScreen - an open-source tool for automated PFAS feature prioritization in non-target HRMS data. Anal. Bioanal. Chem. 2023, 16 (2), 349–362. [DOI: 10.1007/s00216-023-05070-2](https://doi.org/10.1007/s00216-023-05070-2); [github PFAScreen](https://github.com/JonZwe/PFAScreen/blob/main/KMD_analysis.py)

Theoretical calculation of the residual based on the formula:
https://github.com/usnistgov/NISTPFAS/tree/main/suspectlist/fn/calculate_residual.R

-------------------------------------------------------
This is a collaborative repository based on the repository by EcoChem-OSU/KMD at https://github.com/EcoChem-OSU/KMD.git

Updated codes for calculation of Kendrick Mass Default
