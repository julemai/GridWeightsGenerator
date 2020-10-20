# Grid Weights Generator

Generates weights of grid cells contributing to individual shapes. The grid weights are, for example, helpful to spatially aggregate a gridded model output to a irregular shaped domain or, in general, to map grid cells to a polygon shape. 

In this case here, any given gridded model output (see examples in `example/input_*/*.nc`) shall be mapped to a subbasin discretization for the purpose of routing. The subbasin discretization is derived using the RoutingToolbox described by Han et al. (2020a), Han et al. (2020b) and available through Han et al. (2020c). Parts of the Routing Toolbox output is provided here for convenience for an example watershed draining to the WSC streamflow gauge station 02LE024 (see shapefiles in `example/maps/*.shp`).

The generated file containing the grid weights can be used in the hydrologic modeling framework [Raven](http://raven.uwaterloo.ca) to handle gridded NetCDF inputs and map them to the subbasins/ HRUs a distributed model is based on. The grid weights file produced here (based on the Routing Toolbox discretization) is then primarily used to run Raven in routing-mode and route gridded model outputs (from any model).

## Requirements

The script is tested under Pythin 3.8 and it requires the following Python packages:
* `pip install GDAL` (gets you `osgeo`)
* `pip install argparse`
* `pip install numpy`
* `pip install geopandas`
* `pip install netCDF4`

## Usage: General

```python
python derive_grid_weights.py -i [your-nc-file] -d [dim-names] -v [var-names] -r [tlbx-shp-file] -s [subbasin-id] -b [gauge-id] -o [weights-file]	
```
with
```
[your-nc-file]  ... [in] filename of NetCDF file
[dim-names]     ... [in] names of NetCDF dimensions of longitude (x) 
                         and latitude (y) in this order, 
                         e.g. 'rlon,rlat'
[var-names]     ... [in] names of 2D NetCDF variables containing 
                         longitudes and latitudes (in this order) of 
                         centroids of grid cells, 
                         e.g. 'lon,lat'
[tlbx-shp-file] ... [in] name of shapefile routing toolbox provides; 
                         shapefile contains shapes of all land and lake HRUs,
                         e.g. 'gauges_catchment_HRU_withlake.shp'
[subbasin-id]   ... [in] (either this or [gauge-id] must be set)
                         ID of subbasin  most downstream (likely a subbasin 
                         that contains a streamflow gauge station but can 
                         be any subbasin ID); script will include all subbasins 
                         upstream of the given subbasin automatically; according 
                         attribute in [tlbx-shp-file] is called "SubId"; 
                         e.g. '7202'
[gauge-id]      ... [in] (either this or [subbasin-id] must be set)
                         ID of streamflow gauging station; according attribute in 
                         [tlbx-shp-file] is called "Obs_NM"; 
                         e.g. '02LE024'
[weights-file]	... [in] name of file where derived grid weights 
                         should be stored, 
                         e.g. 'GridWeights.txt'
```

## Usage: VIC example

Using _coarse_ routing discretization (see file used for `-r`):
```python
python derive_grid_weights.py -i example/input_VIC/VIC_streaminputs.nc -d 'lon,lat' -v 'lon,lat' -r example/maps/HRUs_coarse.shp -b 02LE024 -o example/input_VIC/GridWeights.txt
```

Using _fine_ routing discretization (see file used for `-r`):
```python
python derive_grid_weights.py -i example/input_VIC/VIC_streaminputs.nc -d 'lon,lat' -v 'lon,lat' -r example/maps/HRUs_fine.shp -b 02LE024 -o example/input_VIC/GridWeights.txt
```

## Usage: MESH example

Using _coarse_ routing discretization (see file used for `-r`):
```python
python derive_grid_weights.py -i example/input_MESH/RFF_H_GRD.nc -d 'rlon,rlat' -v 'longitude,latitude' -r example/maps/HRUs_coarse.shp -b 02LE024 -o example/input_MESH/GridWeights_RFF_H_GRD.txt
```

Using _fine_ routing discretization (see file used for `-r`):
```python
python derive_grid_weights.py -i example/input_MESH/RFF_H_GRD.nc -d 'rlon,rlat' -v 'longitude,latitude' -r example/maps/HRUs_fine.shp -b 02LE024 -o example/input_MESH/GridWeights_RFF_H_GRD.txt
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
