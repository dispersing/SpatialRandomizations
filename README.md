# SpatialRandomizations
This repository will contain files used to conduct complete and structured statial randomizations from a study conducted by JW Dittel, myself, and SB Vadner Wall.

There are 4 files in total, 2 pertaining to the simulations and 2 for creating figures
1. SpatialRadomizations.R
2. SpatialRadomizations_GRF.R
3. PHatFigures.R
4. PHatFigures_GRF.R

**SpatialRadomizations.R** runs a complete randomization and **PHatFigures.R** plots the complete randomization data. Similarly, **SpatialRadomizations_GRF.R** runs a structured randomization (**GRF** refers to the **G**aussian **R**andom **F**ields used from the package [RandomFields](https://cran.r-project.org/web/packages/RandomFields/index.html)) and **PHatFigures_GRF.R** plots the structured randomization data.

The plots for the complete randomization only include the p-hat value.  I didn't clean them up because they were completely uninformative (yes, spatial points are not completely random---we know that!).  The structured randomizations, however, include the observed test stastic (Spearman's rho, in our case), the p-hat value, and color under the density plot that corresponds to how extreme our observations is, given our data (red is more extreme; blue is less extreme).  Each of the complete and structured randomization correlated the 4 major groups: 1. abiotic v. animal variables, 2. abiotic v. plant variables, 3. abiotic v. difference variables, and 4. animal v. plant variables.  Each of these 4 groups have subgroups, and these compose the panels in the graphs.

The directory tree, for the time being, should look like:

wd  
+--SpatialRadomizations.R  
+--PHatFigures.R  
+--SpatialRadomizations_GRF.R  
+--PHatFigures_GRF.R  
+--Data  
|&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;+--GRFData  
+--GRDs  
|&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;+--AbioticGRDs  
|&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;+--AnimalGRDs  
|&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;+--DifferenceGRDs  
|&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;+--NorthAmericaGRDs  
|&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;+--PlantGRDs  
+--Figures  
|&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;+--GRFFigs

Each of the .R files allows users to speficy the directory (e.g., wd, fig.wd, save.wd) as one of the first commands that correspond to the directory above.

More detials and a bit of reorganization will come, but the .R files are here for the time being.
