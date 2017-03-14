# SpatialRandomizations
This repository will contain files used to conduct complete and structured statial randomizations from a study conducted by JW Dittel, myself, and SB Vadner Wall.

There are 4 files in total, 2 pertaining to the simulations and 2 for creating figures
1. SpatialRadomizations.R
2. SpatialRadomizations_GRF.R
3. PHatFigures.R
4. PHatFigures_GRF.R

**SpatialRadomizations.R** runs a complete randomization and **PHatFigures.R** plots the complete randomization data. Similarly, **SpatialRadomizations_GRF.R** runs a structured randomization (**GRF** refers to the **G**aussian **R**andom **F**ields used from the package [RandomFields](https://cran.r-project.org/web/packages/RandomFields/index.html)) and **PHatFigures_GRF.R** plots the structured randomization data.  The directory tree, for the time being should look like:

wd
+--SpatialRadomizations.R
+--PHatFigures.R
+--SpatialRadomizations_GRF.R
+--PHatFigures_GRF.R
+--Data
|  +--GRFData
+--GRDs
|  +--AbioticGRDs
|  +--AnimalGRDs
|  +--DifferenceGRDs
|  +--NorthAmericaGRDs
|  +--PlantGRDs
+--Figures
|  +--GRFFigs

Each of the .R files allows users to speficy the directory (e.g., wd, fig.wd, save.wd) as one of the first commands.


More detials will come, but the .R files are here for the time being.
