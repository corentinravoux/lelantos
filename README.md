# lslyatomo

Package to create and manipulate on large-scale Lyman-alpha tomographic maps.

It includes launchers and scripts to modify


## Install

This package can be install with pip:

```
pip install . --user
```

The minimal required python packages are: fitsio, numpy, scipy, matplotlib
They are automaticaly installed if pip is used.

To convert a tomographic box into a format for 3D visualization, you need pyevtk as additional package.

This package use code parts from public packages [picca](https://github.com/igmhub/picca) and [SaclayMocks](https://github.com/igmhub/SaclayMocks) but it is advised to install them.





Finally to run tomographic mapping, you need to put the executable of the tomographic algorithm you want to use in the lslyatomo/exec directory. From now the only algorithm available is [dachshund](https://github.com/caseywstark/dachshund)





## Launchers

/local/home/cravoux/Documents/Python/Packages/lslyatomo/launchers/lslyatomo_convert_delta.py
/local/home/cravoux/Documents/Python/Packages/lslyatomo/launchers/lslyatomo_plot_delta.py

/local/home/cravoux/Documents/Python/Packages/lslyatomo/launchers/lslyatomo_plot_pixel.py

/local/home/cravoux/Documents/Python/Packages/lslyatomo/launchers/lslyatomo_plot_stack_qso.py
/local/home/cravoux/Documents/Python/Packages/lslyatomo/launchers/lslyatomo_plot_stack_void.py
/local/home/cravoux/Documents/Python/Packages/lslyatomo/launchers/lslyatomo_plot_tomography.py


- modify flux contrast files (delta)
- plot flux contrast properties
- convert flux contrast files to pixel files

## lslyatomo


The launcher present in this folder were used for the creation of a Tomographic map obtained from DR16 eBOSS data.

This folder contains :

-	/local/home/cravoux/Documents/Python/Launchers/Tomography/Launch_3D.py

	Program to launch a 3D visualization compatible with Paraview of the box obtained from Tomography


-	/local/home/cravoux/Documents/Python/Launchers/Tomography/Launch_Dachshund.py

	Launcher of Tomographic process on cluster (Task manager)


- 	/local/home/cravoux/Documents/Python/Launchers/Tomography/Launch_Prepare_Daschund.py

	Use picca to prepare all dachshund inputs from deltas


- 	/local/home/cravoux/Documents/Python/Launchers/Tomography/Launch_Ra_Dec_Map.py

	Plot of delta maps obtained at the end of picca in P1D


-	/local/home/cravoux/Documents/Python/Launchers/Tomography/Launch_Stack.py

	Stacking procedure on tomographic maps


-	/local/home/cravoux/Documents/Python/Launchers/Tomography/Launch_Treat_box.py

	Program which analyze Saclay mocks boxes to extract DM density field corresponding to a given tomographic map


-	/local/home/cravoux/Documents/Python/Launchers/Tomography/Launch_TreatClamato.py

	Plot properties from a tomographic map


-	/local/home/cravoux/Documents/Python/Launchers/Tomography/Launch_Treat_stack.py

	Plot properties from stacks of the tomographic map


-	/local/home/cravoux/Documents/Python/Launchers/Tomography/Launch_Treat_Voids.py

	Plot properties from voids or proto-clusters obtained on a tomographic map


-	/local/home/cravoux/Documents/Python/Launchers/Tomography/Launch_VoidFinder.py

	Void and Over-density finder for Tomographic map.



A classical Tomographic map study generally follows :

1) Preparing the datas obtaining from picca (deltas) to prepare the DACHSHUND inputs: Launch_ra_dec_map.py, Launch_Prepare_Daschund.py

2) Launching of DACHSHUND: Launch_Dachshund.py

3) Verifications of the Tomographic map with the initial DM density field (on mocks): Launch_Treat_box.py

4) Supplementary calculations on the map : Launch_VoidFinder.py, Launch_Stack.py

5) Ploting of the tomographic map : Launch_TreatClamato.py, Launch_3D.py

6) Ploting of supplementary studies : Launch_Treat_Voids.py, Launch_Treat_stack.py
