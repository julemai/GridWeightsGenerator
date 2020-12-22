# Example

This folder contains all data required to run the _GridWeightsGenerator_.

* `input_ERA5`:<br> contains example gridded NetCDF data ([ERA5](https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5)) which has 1D lat/lon coordinates
* `input_MESH`:<br> contains example gridded NetCDF model outputs derived with the [MESH](https://wiki.usask.ca/display/MESH/About+MESH) model
* `input_SWAT`:<br> contains example [SWAT](https://swat.tamu.edu) basin descritization with `nbasins` subbasins as shapefile; model output will then need to be in NetCDF with a 2D variable (`time` x `nbasins`)
* `input_VIC`:<br> contains example gridded NetCDF model outputs derived with the [VIC](https://vic.readthedocs.io/en/master/) model
* `maps`:<br> contains example outputs of the Routing Toolbox specifying the routing network and subbasin discretization
