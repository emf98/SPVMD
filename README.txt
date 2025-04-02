Welcome!

This repo is for code associated with the varied calculation of stratospheric polar vortex ellipse geometries/metrics for application to analyses. 
Code from this folder supports a manuscript submitted to Monthly Weather Review by Fernandez, Attard, and Lang 2025. 

Below are descriptions of the contents of this folder:

-> `fitEllipse3_new.py` contains code for mathematical definition statements related to the full geometric calculation of the ellipse metrics.

-> 'EllipseDef_ERA5.py' contains a definition statement code supporting ellipse metric calculations for ERA-5 datasets.
This code can be used at most pressure levels and can contour other features of interest through gph beyond the stratospheric polar vortex.

-> 'EllipseDef_ERAI.py' contains a definition statement code supporting ellipse metric calculations for ERA-I datasets.

Please check these codes when using them to update the locations of your saved ERA-5 or ERA-I files. 
Additionally, the plotting portion of the code may return issues when attempting to look at other areas globally if you do not modify the cartopy distinctions. 

-> 'Fernandez_etal25' folder contains code used for calculating the diagnostics and associated images from the manuscript.

-> 'S2S_Ellipses' folder contains files relevant to calculating these diagnostics with S2S datasets. There is a separate README in this folder. 

Any questions regarding the code files may be directed to Elena Fernandez (emfernandez@albany.edu/elenamf98@gmail.com). 