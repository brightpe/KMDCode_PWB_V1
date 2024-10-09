# Kendrick Mass Defect Analysis for High-Resolution mass spectrometry data
-----------------------------------------------------------------------------------
Description
-----------
This code performs a Kendrick Mass Defect (KMD) analysis of user input mass spec data by creating a feature list via the XCMS peak picking algorithm
The workflow combines existing code (3 steps, see below in reference section) that might be run independently.

The steps (3 individual scripts) are the following:  
- 1-PatRoon_XCMS_feature_opti: optimizes the parameters used for XCMS peak picking  
- 2-PatRoon_featuregroup: Performs peak picking using XCMS (to fascilitate KMD analysis in step 3)  
- 3-KMD_analysis: Scan featurelist for homologous series (perform a KMD analysis)
- creat_KMD_susp_list: create a specific suspect list to match wtih the KMD result and help the identification of unknown exact masses.  
> [!CAUTION]  
> **DO NOT CONSIDER the creat_KMD_susp_list match as a comprehensive suspect screening method. Results need to be verified**


Dependency and installation
----------------
R(>=4.3.1)    

RTools 

patRoon(>=2.3.0) Check [patron handbook](https://rickhelmus.github.io/patRoon/handbook_bd/index.html)
**Optional**: Text editor such as Rstudio or Visual Studio Code 

MetaboCoreUtils. To install this package, start R and enter:  

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("MetaboCoreUtils")
```

Getting started
----------------
1. Download the codes to your local machine, create a directory (folder), and put all codes in the directory
2. Create two folders in your directory called "input" and "ouptut"
3. In the input folder, put: (1) a sample list (default name is sample_list_KMD), (2) your suspect list (see downloads for: ), (3) all datafiles converted from .raw file format to .mzXML (can convert using MSConvert or another software)
4. Run the scripts in order, taking care to update the parameters in script 2 (2-PatRoon_featuregroup) with the optimized peak picking parameters identified in script 1 (1-PatRoon_XCMS_feature_opti), and to have the featurelist generated from script 2 in the output folder before running script 3 (3-KMD_analysis)


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
