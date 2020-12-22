# Grid Weights Generator for Raven

Generates weights of grid cells contributing to individual shapes. The grid weights are, for example, helpful to spatially aggregate a gridded model output to a irregular shaped domain or, in general, to map grid cells to a polygon shape. 

In this case here, any given gridded model output (see examples in `example/input_*/*.nc`) shall be mapped to a subbasin discretization for the purpose of routing. The subbasin discretization is derived using the RoutingToolbox described by Han et al. (2020a), Han et al. (2020b) and available through Han et al. (2020c). Parts of the Routing Toolbox output is provided here for convenience for an example watershed draining to the WSC streamflow gauge station 02LE024 (see shapefiles in `example/maps/*.shp`).

<p align="center">
   <img src="https://github.com/julemai/GridWeightsGenerator/wiki/images/logos/RavenBanner.png" width="45%" />
</p>

The generated file containing the grid weights can be used in the hydrologic modeling framework [Raven](http://raven.uwaterloo.ca) to handle gridded NetCDF inputs and map them to the subbasins/ HRUs a distributed model is based on. The grid weights file produced here (based on the Routing Toolbox discretization) is then primarily used to run Raven in routing-mode and route gridded model outputs (from any model).

## Requirements

The script is tested under Python 3.6 and 3.8 on MacOS and Windows. It requires the following Python packages:
* `pip install GDAL` (gets you `osgeo`)
* `pip install argparse`
* `pip install numpy`
* `pip install geopandas`
* `pip install netCDF4`

In case you are working of ComputeCanada's cluster systems, e.g., Graham you can find instructions on how to setup those Python environments [here](https://github.com/julemai/GridWeightsGenerator/blob/main/misc/setup-python-env-graham.md).

## Usage: General

```python
python derive_grid_weights.py -i [your-nc-file] -d [dim-names] -v [var-names] -r [tlbx-shp-file] -s [subbasin-id] -b [gauge-id] -f [shp-attr-name] -o [weights-file]	
```
with
```
[your-nc-file]   ... [in] filename of NetCDF file or shapefile that contains how 
                          you discretized your model
[dim-names]      ... [in] names of NetCDF dimensions of longitude (x) 
                          and latitude (y) in this order, 
                          e.g. "rlon,rlat"
[var-names]      ... [in] names of 2D NetCDF variables containing 
                          longitudes and latitudes (in this order) of 
                          centroids of grid cells, 
                          e.g. "lon,lat"
[tlbx-shp-file]  ... [in] name of shapefile routing toolbox provides; 
                          shapefile contains shapes of all land and lake HRUs,
                          e.g. "HRUs.shp"
[subbasin-id]    ... [in] (either this or [gauge-id] must be set)
                          ID of subbasin  most downstream (likely a subbasin 
                          that contains a streamflow gauge station but can 
                          be any subbasin ID); script will include all subbasins 
                          upstream of the given subbasin automatically; according 
                          attribute in [tlbx-shp-file] is called "SubId"; 
                          e.g. "7202"
[gauge-id]       ... [in] (either this or [subbasin-id] must be set)
                          ID of streamflow gauging station; according attribute in 
                          [tlbx-shp-file] is called "Obs_NM"; 
                          e.g. "02LE024"
[weights-file]	  ... [in] name of file where derived grid weights 
                          should be stored, 
                          e.g. "GridWeights.txt"
[shp-attr-name]  ... [in] in case your model discretization (-i) is a shapefile you 
                          need to provide the attribute name that contains the 
                          numbering from [0...N-1] of your subbasins (this order 
                          needs to match the order of your model outputs in NetCDF file)
```

## Usage: VIC example

Using _coarse_ routing discretization (see file used for `-r`):
```python
python derive_grid_weights.py -i example/input_VIC/VIC_streaminputs.nc -d "lon,lat" -v "lon,lat" -r example/maps/HRUs_coarse.shp -b 02LE024 -o example/input_VIC/GridWeights.txt
```

Using _fine_ routing discretization (see file used for `-r`):
```python
python derive_grid_weights.py -i example/input_VIC/VIC_streaminputs.nc -d "lon,lat" -v "lon,lat" -r example/maps/HRUs_fine.shp -b 02LE024 -o example/input_VIC/GridWeights.txt
```

## Usage: MESH example

Using _coarse_ routing discretization (see file used for `-r`):
```python
python derive_grid_weights.py -i example/input_MESH/RFF_H_GRD.nc -d "rlon,rlat" -v "longitude,latitude" -r example/maps/HRUs_coarse.shp -b 02LE024 -o example/input_MESH/GridWeights_RFF_H_GRD.txt
```

Using _fine_ routing discretization (see file used for `-r`):
```python
python derive_grid_weights.py -i example/input_MESH/RFF_H_GRD.nc -d "rlon,rlat" -v "longitude,latitude" -r example/maps/HRUs_fine.shp -b 02LE024 -o example/input_MESH/GridWeights_RFF_H_GRD.txt
```

## Usage: ERA5 example (1D lat/lon coordinates = regular grid)

Using _coarse_ routing discretization (see file used for `-r`):
```python
python derive_grid_weights.py -i example/input_ERA5/era5-crop.nc -d "longitude,latitude" -v "longitude,latitude" -r example/maps/HRUs_coarse.shp -b 02LE024 -o example/input_ERA5/GridWeights_ERA5.txt
```

Using _fine_ routing discretization (see file used for `-r`):
```python
python derive_grid_weights.py -i example/input_ERA5/era5-crop.nc -d "longitude,latitude" -v "longitude,latitude" -r example/maps/HRUs_fine.shp -b 02LE024 -o example/input_ERA5/GridWeights_ERA5.txt
```

## Usage: SWAT example (model discretization provided as shapefile not NetCDF grid)

Using _coarse_ routing discretization (see file used for `-r`):
```python
python derive_grid_weights.py -i example/input_SWAT/OTT_sub.shp -f "NetCDF_col" -r example/maps/HRUs_coarse.shp -b 02LE024 -o example/input_SWAT/GridWeights_SWAT.txt
```

For several subbasins there is a mismatch between the provided domain of SWAT (`OTT_sub.shp`) and the domain the routing product requires (`HRUs_coarse.shp`). Basins with large errors are reported on the command-line when you run the script. To force the script to resolve these errors by rescaling the weights available such that they add up to one you can use option `-e` and set a threshold to which errors are exceptable for you (default: 5% = 0.05). To resolve all previously reported basins use:
```python
python derive_grid_weights.py -i example/input_SWAT/OTT_sub.shp -f "NetCDF_col" -r example/maps/HRUs_coarse.shp -b 02LE024 -o example/input_SWAT/GridWeights_SWAT.txt -e 0.42
```

Using _fine_ routing discretization (see file used for `-r`):
```python
python derive_grid_weights.py -i example/input_SWAT/OTT_sub.shp -f "NetCDF_col" -r example/maps/HRUs_fine.shp -b 02LE024 -o example/input_SWAT/GridWeights_SWAT.txt
```


## References

Han, M., Mai, J., Tolson, B. A., Craig, J. R., Gaborit, É., Liu, H., and Lee, K. (2020a): <br>
_Subwatershed-based lake and river routing products for hydrologic and land surface models applied over Canada_  <br>
Canadian Water Resources Journal, 0, 1-15. ([publication](https://doi.org/10.1080/07011784.2020.1772116))

Han, M. et al. (2020b): <br>
_An automated GIS toolbox for watershed delineation with lakes_  <br>
In preparation.  

Han, M., Mai, J., Tolson, B. A., Craig, J. R., Gaborit, É., Liu, H., and Lee, K. (2020c): <br>
_A catchment-based lake and river routing product for hydrologic and land surface models in Canada (Dataset)_  <br>
Zenodo. ([dataset](https://doi.org/10.5281/zenodo.3667677))
