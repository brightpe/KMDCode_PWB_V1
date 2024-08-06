# Kendrick Mass Defect Analysis for High-Resolution mass spectrometry data
-----------------------------------------------------------------------------------
Description
-----------
Perfomed Kendrick Mass Defect (KMD) analysis from a feature list previously extracted from a peak picking algorithm (e.g., XCMS)
The workflow combine existing code (step) (see below in reference section) that might be run independently.

The step (scripts) are the following:  
1-PatRoon_XCMS_feature_opti: optimized the parameter for XCMS peak picking  
2-PatRoon_featuregroup: Peak picking script prior to KMD analysis  
3-KMD_analysis: KMD analysis script   
creat_KMD_susp_list: create a specific KMD suspect list to match the KMD result with some suspect and help the identification of unknow.  
> [!CAUTION]  
> **DO NOT CONSIDER THIS AS A TRUE SUSPECT SCREENING STEP**


Dependency
----------------
R(>=4.3.1)
patRoon(>=2.3.0)

Getting started
----------------
1.
2.


References
------

The current workflow combine/translate divers code from the following sources

Peak picking with XCMS:  
Helmus, R.; van de Velde, B.; Brunner, A. M.; ter Laak, T. L.; van Wezel, A. P.; Schymanski, E. L. patRoon 2.0: Improved non-target analysis workflows including automated transformation product screening. Journal of Open Source Software 2022, 7 (71), 4029. [DOI: 10.21105/joss.04029](https://doi.org/10.21105/joss.04029); [github patRoon](https://github.com/rickhelmus/patRoon).

KMD analysis:  
Zweigle, J.; Bugsel, B.; Fabregat-Palau, J.; Zwiener, C. PFDeltaScreen - an open-source tool for automated PFAS feature prioritization in non-target HRMS data. Anal. Bioanal. Chem. 2023, 16 (2), 349–362. [DOI: 10.1007/s00216-023-05070-2](https://doi.org/10.1007/s00216-023-05070-2); [github PFAScreen](https://github.com/JonZwe/PFAScreen/blob/main/KMD_analysis.py)

Theoretical calculation of the residual based on the formula:  [gitHub NIST](https://github.com/usnistgov/NISTPFAS/tree/main/suspectlist/fn/calculate_residual.R)

Citing
------
The full work flow is described in:  
coming soon

Other examples using the workflow include:  



License
-------
This work is licensed under a
[Creative Commons Attribution-NonCommercial 4.0 International License][cc-by-nc].

[![CC BY-NC 4.0][cc-by-nc-image]][cc-by-nc]

[cc-by-nc]: https://creativecommons.org/licenses/by-nc/4.0/
[cc-by-nc-image]: https://licensebuttons.net/l/by-nc/4.0/88x31.png
[cc-by-nc-shield]: https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg


-------------------------------------------------------
This is a collaborative repository based on the repository by EcoChem-OSU/KMD at https://github.com/EcoChem-OSU/KMD.git

Updated codes for calculation of Kendrick Mass Default
