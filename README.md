# LELANTOS - Large-scalE LymAN-alpha TOmography Suite

LELANTOS - Greak titan associated to hunt and air (name found after discussion with [E. Chaussidon](https://github.com/echaussidon))

Package to create and manipulate large-scale Lyman-alpha tomographic maps.

It includes launchers and scripts to modify or process diverse tomographic objects including flux contrast files (delta) lines-of-sight (pixel), tomographic map and void catalog.


## Install

This package can be install with pip:

```
pip install . --user
```

The minimal required python packages are: fitsio, numpy, scipy, matplotlib
They are automaticaly installed if pip is used.

To convert a tomographic box into a format for 3D visualization, you need pyevtk as additional package.

This package use code parts from public packages [picca](https://github.com/igmhub/picca) and [SaclayMocks](https://github.com/igmhub/SaclayMocks) but it is advised to install them.


Finally to run tomographic mapping, you need to put the executable of the tomographic algorithm you want to use in the lelantos/lelantos/exec directory. From now the only algorithm available is [dachshund](https://github.com/caseywstark/dachshund).



## Launchers

The main launchers available in lelantos/launchers folder are:


- **lelantos_subsample_delta.py** Subsampling of flux constrat files (delta), output of picca delta generation routine.

- **lelantos_shuffle_delta.py** Shuffling of flux constrat files (delta).

- **lelantos_convert_delta.py** Conversion of flux constrat files (delta) to pixel file, input of tomographic mapping procedure.

- **lelantos_tomography_manager.py** Running tomography algorithm.

- **lelantos_tomography_process.py** Processing the tomographic map output.

- **lelantos_void_finder.py** Launching of the void finding routine on the tomographic map.

- **lelantos_void_process.py** Post-processing of void catalog.

- **lelantos_plot_delta.py** Plotting routines for flux constrat files (delta).

- **lelantos_plot_pixel.py** Plotting routines for pixel file.

- **lelantos_plot_tomography.py** Plotting routines for tomographic map.

- **lelantos_plot_void.py** Plotting routines for void catalog.




Additional launchers can be used for specific studies:

- **lelantos_stack_void.py** Stacking of void position on a tomographic map.

- **lelantos_stack_qso.py** Stacking of quasar position on a tomographic map.

- **lelantos_plot_stack_void.py** Plotting routines for void stack.

- **lelantos_plot_stack_qso.py** Plotting routines for qso stack.

- **lelantos_tomography_3d.py** Conversion of tomographic objects into a format for 3D visualization.

- **lelantos_boxdm.py** Extraction of the underlying dark matter box of mocks associated to a given box geometry.


## Scripts

The lelantos/scripts folder contains scripts using argparse or config parser.

Most of the previous launchers can be globaly launched by modifying interface_example.ini and running:

```
lelantos_interface.py interface_example.ini
```

The other scripts are used to quickly modify or show information of tomographic objects.
