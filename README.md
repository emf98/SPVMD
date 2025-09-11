# Welcome!

**SPVDM** is a repository containing code for the best-fit calculation of stratospheric polar vortex ellipse geometries/metrics. This method provides complementary metrics for determining stratospheric polar vortex strength and SSW variability to those established in [*Seviour et al. 2013*](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/grl.50927).

This repository supports a manuscript submitted to the AMS Journal of Applied Meteorology and Climatology by Fernandez, Lang, and Attard 2025 (*under review*). 

Required Python packages include `numpy`, `datetime`, and `xarray`.

## Contents

* `fitEllipse3_new.py` contains code for mathematical definition statements related to the full geometric calculation of the ellipse metrics.

* `EllipseDef_ERA5.py` contains a definition statement code supporting ellipse metric calculations using the best-fit method for ERA-5 datasets.
    * This code can be used at most pressure levels and can contour other features of interest through a user-defined geopotential height value. The manuscript for this repository supports the use of a climatologically representative contour in defining the vortex edge.

* `EllipseDef_ERAI.py` contains a definition statement code supporting ellipse metric calculations for ERA-I datasets.

Please check `EllipseDef_ERA5.py` and `EllipseDef_ERAI.py` to ensure your correctly saved ERA-5 or ERA-I file locations are used. 

Additionally, the code uses an older version (Python 3) of cartopy for plotting circumpolar views of the vortex. If the plotting lines return issues, and you do not wish to modify the cartopy distinction, comment them out or remove them. It will not affect the running of the code. 

### Example Code

* `EllipseDef_SaveExample.ipynb` provides an example for calculating the ellipse metrics for a single month/year. A single instance of this code should not exceed five minutes. 

* `S2S_Ellipses` contains files relevant to calculating these diagnostics with S2S datasets. There is a separate README in this folder. 


## Citation

If you use the Python scripts in this repository, please cite:

Fernandez, E.M., A.L. Lang, and H.E. Attard (2025): “Stratospheric Polar Vortex Ellipse Diagnostics for Realtime and S2S Forecast Analyses.” *Submitted for Peer Review to the Journal of Applied Meteorology and Climatology*. 

Any questions regarding the code files may be directed to Elena Fernandez (emfernandez@albany.edu/elenamf98@gmail.com). 