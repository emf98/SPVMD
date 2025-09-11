# Welcome!

**SPVDM** is a repository containing code for the best-fit calculation of stratospheric polar vortex ellipse geometries/metrics. This method provides complementary metrics for determining stratospheric polar vortex strength and SSW variability to those established in [*Seviour et al. 2013*](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/grl.50927).

This repository supports a manuscript submitted to the AMS Journal of Applied Meteorology and Climatology by Fernandez, Attard, and Lang 2025 (*under review*). 

Below are the contents of this repository:

* `fitEllipse3_new.py` contains code for mathematical definition statements related to the full geometric calculation of the ellipse metrics.

* `EllipseDef_ERA5.py` contains a definition statement code supporting ellipse metric calculations using the best-fit method for ERA-5 datasets.
    * This code can be used at most pressure levels and can contour other features of interest through a user-defined geopotential height value. The manuscript for this repository supports the use of a climatologically representative contour in defining the vortex.

* `EllipseDef_ERAI.py` contains a definition statement code supporting ellipse metric calculations for ERA-I datasets.

Please check these codes when using them to update the locations of your saved ERA-5 or ERA-I files. 
Additionally, the plotting portion of the code may return issues when attempting to look at other areas globally if you do not modify the cartopy distinctions. 

-> 'EllipseDef_SaveExample.ipynb' provides an example for calculating the ellipse metrics. 

-> 'Fernandez_etal25' folder contains code used for calculating the diagnostics and some of the associated images from the manuscript.

-> 'S2S_Ellipses' folder contains files relevant to calculating these diagnostics with S2S datasets. There is a separate README in this folder. 

Any questions regarding the code files may be directed to Elena Fernandez (emfernandez@albany.edu/elenamf98@gmail.com). 