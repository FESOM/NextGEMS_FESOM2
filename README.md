# NextGEMS FESOM2 examples

This repository holds a collection of notebooks to plot FESOM2 and ICON ocean data (e.g. from the DYAMOND and NextGEMS projects).

- **FESOM2_ICON_grids_easy_plot_and_interpolate.ipynb** - plot and interpolate FESOM2 and ICON data without complicated dependencies (except cartopy, that is optional).
- **FESOM2_ICON_ploting_on_original_grid.ipynb** - plot small regions of FESOM2 and ICON on original unstructured triangular grid.
- **unrotate_UV** - folder that contains examples and utilities to unrotate FESOM2 vector data, ony needed if you want to work with velocity components.
- **arctic1km** - examples of how to work with Arctic 1km resolution FESOM2 data (Wang et al., 2020). Require pyfesom2 installation.

Installation and Dependencies
=============================

Most of the HPC centers will have everything (maybe exept cartopy) installed in the basic python module. For local/personal installation, you should be fine with the installation steps described [here](https://github.com/koldunovn/python_for_geosciences#getting-started-for-linuxmac).

Instructions for [installation of pyfesom2](https://github.com/FESOM/pyfesom2#development-installation), that is needed for some of the examples.

For examples that use containers, you can look on [how to pull and run the containers here](https://github.com/FESOM/FESOM2_Docker/tree/master/pyfesom2#installation-of-pyfesom2-and-fdiag). 

