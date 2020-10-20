#!/usr/bin/env python
from __future__ import print_function

# Copyright 2016-2020 Juliane Mai - juliane.mai(at)uwaterloo.ca
#
# License
# This file is part of Juliane Mai's personal code library.
#
# Juliane Mai's personal code library is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Juliane Mai's personal code library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with Juliane Mai's personal code library.  If not, see <http://www.gnu.org/licenses/>.
#

# TESTING
#
# --------------------------------------------------
# GRIP-GL version with streamflow gauge ID given (attribute "Obs_NM" in shapefile) --> -b 02LE024
# --------------------------------------------------

# ------------------------
#        VIC
# ------------------------
#    run derive_grid_weights.py -i example/input_VIC/VIC_streaminputs.nc -d 'lon,lat' -v 'lon,lat' -r example/maps/HRUs.shp -b 02LE024 -o example/input_VIC/GridWeights.txt

# ------------------------
#        MESH
# ------------------------
#    run derive_grid_weights.py -i example/input_MESH/RFF_H_GRD.nc      -d 'rlon,rlat' -v 'longitude,latitude' -r example/maps/HRUs.shp -b 02LE024 -o example/input_MESH/GridWeights_RFF_H_GRD.txt
#    run derive_grid_weights.py -i example/input_MESH/DRAINSOL_H_GRD.nc -d 'rlon,rlat' -v 'longitude,latitude' -r example/maps/HRUs.shp -b 02LE024 -o example/input_MESH/GridWeights_DRAINSOL_H_GRD.txt

# --------------------------------------------------
# GRIP-GL version with subbasin ID given (attribute "SubId" in shapefile)  --> -s 7202
# --------------------------------------------------

# ------------------------
#        VIC
# ------------------------
#    run derive_grid_weights.py -i example/input_VIC/VIC_streaminputs.nc -d 'lon,lat' -v 'lon,lat' -r example/maps/HRUs.shp -s 7202 -o example/input_VIC/GridWeights.txt

# ------------------------
#        MESH
# ------------------------
#    run derive_grid_weights.py -i example/input_MESH/RFF_H_GRD.nc      -d 'rlon,rlat' -v 'longitude,latitude' -r example/maps/HRUs.shp -s 7202 -o example/input_MESH/GridWeights_RFF_H_GRD.txt
#    run derive_grid_weights.py -i example/input_MESH/DRAINSOL_H_GRD.nc -d 'rlon,rlat' -v 'longitude,latitude' -r example/maps/HRUs.shp -s 7202 -o example/input_MESH/GridWeights_DRAINSOL_H_GRD.txt






# command line arguments
import argparse

# to perform numerics
import numpy as np             

# read shapefiles and convert to GeoJSON and WKT
import geopandas as gpd

# get equal-area projection and derive overlay
from   osgeo   import ogr
from   osgeo   import osr

# read netCDF files
import netCDF4 as nc4



input_file  = 'example/input_VIC/VIC_streaminputs.nc'
dimname     = ['rlon','rlat']
varname     = ['lon','lat']
routinginfo = 'example/maps/HRUs.shp'
basin       = None  # e.g. '02LE024'
SubId       = None  # e.g. 7202
output_file = 'GriddedForcings.txt'

parser      = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
              description='''Convert files from ArcGIS raster format into NetDF file usable in CaSPAr.''')
parser.add_argument('-i', '--input_file', action='store',
                    default=input_file, dest='input_file', metavar='input_file', 
                    help='Example NetCDF file containing at least 2D latitudes and 2D longitudes. Grid needs to be representative of model outputs that are then required to be routed.')
parser.add_argument('-d', '--dimname', action='store',
                    default=dimname, dest='dimname', metavar='dimname',
                    help='Dimension names of longitude (x) and latitude (y) (in this order). Example: "rlon,rlat", or "x,y"')
parser.add_argument('-v', '--varname', action='store',
                    default=varname, dest='varname', metavar='varname',
                    help='Variable name of 2D longitude and latitude variables in NetCDF (in this order). Example: "lon,lat"')
parser.add_argument('-r', '--routinginfo', action='store',
                    default=routinginfo, dest='routinginfo', metavar='routinginfo',
                    help='Shapefile that contains all information of the routing toolbox for the catchment of interest (and maybe some more catchments).')
parser.add_argument('-b', '--basin', action='store',
                    default=basin, dest='basin', metavar='basin',
                    help='Basin of interest (corresponds to "Gauge_ID" in shapefile given with -r). Either this or SubId ID (-s) needs to be given.')
parser.add_argument('-s', '--SubId', action='store',
                    default=SubId, dest='SubId', metavar='SubId',
                    help='SubId of most downstream subbasin (containing usually a gauge station) (corresponds to "SubId" in shapefile given with -r). Either this or basin ID (-b) needs to be given.')
parser.add_argument('-o', '--output_file', action='store',
                    default=output_file, dest='output_file', metavar='output_file', 
                    help='File that will contain grid weights for Raven.')

args          = parser.parse_args()
input_file    = args.input_file
dimname       = np.array(args.dimname.split(','))
varname       = np.array(args.varname.split(','))
routinginfo   = args.routinginfo
basin         = args.basin
SubId         = args.SubId
output_file   = args.output_file

if not(SubId is None):
    SubId = np.int(SubId)

if (SubId is None) and (basin is None):
    raise ValueError("Either gauge ID (option -b; e.g., 02AB003) or SubId ID (option -s; e.g., 7173) specified in shapefile needs to be given. You specified none. This basin will be the most downstream gridweights of all upstream subbasins will be added automatically.")

if ( not(SubId is None) ) and ( not(basin is None) ):
    raise ValueError("Either gauge ID (option -b; e.g., 02AB003) or SubId ID (option -s; e.g., 7173) specified in shapefile needs to be specified. You specified both. This basin will be the most downstream gridweights of all upstream subbasins will be added automatically.")

del parser, args


# better dont chnage that ever
crs_lldeg = 4326        # EPSG id of lat/lon (deg) coordinate referenence system (CRS)
crs_caea  = 3573        # EPSG id of equal-area    coordinate referenence system (CRS)


def create_gridcells_from_centers(lat, lon):
    
    # create array of edges where (x,y) are always center cells
    nlon = np.shape(lon)[1]
    nlat = np.shape(lat)[0]
    lonh = np.empty((nlat+1,nlon+1), dtype=np.float)
    lath = np.empty((nlat+1,nlon+1), dtype=np.float)
    tmp1 = [ [ (lat[ii+1,jj+1]-lat[ii,jj])/2 for jj in range(nlon-1) ] + [ (lat[ii+1,nlon-1]-lat[ii,nlon-2])/2 ] for ii in range(nlat-1) ]
    tmp2 = [ [ (lon[ii+1,jj+1]-lon[ii,jj])/2 for jj in range(nlon-1) ] + [ (lon[ii+1,nlon-1]-lon[ii,nlon-2])/2 ] for ii in range(nlat-1) ]
    dlat = np.array(tmp1 + [ tmp1[-1] ])
    dlon = np.array(tmp2 + [ tmp2[-1] ])
    lonh[0:nlat,0:nlon] = lon - dlon
    lath[0:nlat,0:nlon] = lat - dlat

    # make lat and lon one column and row wider such that all
    lonh[nlat,0:nlon] = lonh[nlat-1,0:nlon] + (lonh[nlat-1,0:nlon] - lonh[nlat-2,0:nlon]) 
    lath[nlat,0:nlon] = lath[nlat-1,0:nlon] + (lath[nlat-1,0:nlon] - lath[nlat-2,0:nlon]) 
    lonh[0:nlat,nlon] = lonh[0:nlat,nlon-1] + (lonh[0:nlat,nlon-1] - lonh[0:nlat,nlon-2]) 
    lath[0:nlat,nlon] = lath[0:nlat,nlon-1] + (lath[0:nlat,nlon-1] - lath[0:nlat,nlon-2]) 
    lonh[nlat,nlon]   = lonh[nlat-1,nlon-1] + (lonh[nlat-1,nlon-1] - lonh[nlat-2,nlon-2]) 
    lath[nlat,nlon]   = lath[nlat-1,nlon-1] + (lath[nlat-1,nlon-1] - lath[nlat-2,nlon-2])

    return [lath,lonh]

def shape_to_geometry(shape_from_jsonfile,epsg=None):
    
    # converts shape read from shapefile to geometry
    # epsg :: integer EPSG code
    
    ring_shape = ogr.Geometry(ogr.wkbLinearRing)
    
    for ii in shape_from_jsonfile:
        ring_shape.AddPoint_2D(ii[0],ii[1])
    # close ring
    ring_shape.AddPoint_2D(shape_from_jsonfile[0][0],shape_from_jsonfile[0][1])

    poly_shape = ogr.Geometry(ogr.wkbPolygon)
    poly_shape.AddGeometry(ring_shape)

    if not( epsg is None):
        source = osr.SpatialReference()
        source.ImportFromEPSG(crs_lldeg)       # usual lat/lon projection

        target = osr.SpatialReference()
        target.ImportFromEPSG(epsg)       # any projection to convert to
    
        transform = osr.CoordinateTransformation(source, target)
        poly_shape.Transform(transform)
        
    return poly_shape

def check_proximity_of_envelops(gridcell_envelop, shape_envelop):

    # checks if two envelops are in proximity (intersect)
    
    # minX  --> env[0]
    # maxX  --> env[1]
    # minY  --> env[2]
    # maxY  --> env[3]

    if  ((gridcell_envelop[0] <= shape_envelop[1]) and (gridcell_envelop[1] >= shape_envelop[0]) and
         (gridcell_envelop[2] <= shape_envelop[3]) and (gridcell_envelop[3] >= shape_envelop[2])):

        grid_is_close = True

    else:

        grid_is_close = False

    return grid_is_close

def check_gridcell_in_proximity_of_shape(gridcell_edges, shape_from_jsonfile):

    # checks if a grid cell falls into the bounding box of the shape
    # does not mean it intersects but it is a quick and cheap way to
    # determine cells that might intersect
    
    # gridcell_edges = [(lon1,lat1),(lon2,lat2),(lon3,lat3),(lon4,lat4)]
    # shape_from_jsonfile

    min_lat_cell  = np.min([ii[1] for ii in gridcell_edges])
    max_lat_cell  = np.max([ii[1] for ii in gridcell_edges])
    min_lon_cell  = np.min([ii[0] for ii in gridcell_edges])
    max_lon_cell  = np.max([ii[0] for ii in gridcell_edges])

    lat_shape = np.array([ icoord[1] for icoord in shape_from_jsonfile ])     # is it lat???
    lon_shape = np.array([ icoord[0] for icoord in shape_from_jsonfile ])     # is it lon???

    min_lat_shape = np.min(lat_shape)
    max_lat_shape = np.max(lat_shape)
    min_lon_shape = np.min(lon_shape)
    max_lon_shape = np.max(lon_shape)

    if  ((min_lat_cell <= max_lat_shape) and (max_lat_cell >= min_lat_shape) and
         (min_lon_cell <= max_lon_shape) and (max_lon_cell >= min_lon_shape)):

        grid_is_close = True

    else:

        grid_is_close = False

    return grid_is_close

# -------------------------------
# Read NetCDF
# -------------------------------
print('   (1) Reading NetCDF data ...')

nc_in = nc4.Dataset(input_file, "r")
lon      = nc_in.variables[varname[0]][:]
lon_dims = nc_in.variables[varname[0]].dimensions
lat      = nc_in.variables[varname[1]][:]
lat_dims = nc_in.variables[varname[1]].dimensions
nc_in.close()

# Raven numbering is:
#
#      [      1      2      3   ...     1*nlon
#        nlon+1 nlon+2 nlon+3   ...     2*nlon
#           ...    ...    ...   ...     ...
#           ...    ...    ...   ...  nlat*nlon ]
#
# --> Making sure shape of lat/lon fields is like that
#
if np.all(np.array(lon_dims) == dimname[::1]):
    lon = np.transpose(lon)
    print('   >>> switched order of dimensions for variable "{0}"'.format(varname[0]))
elif np.all(np.array(lon_dims) == dimname[::-1]):
    print('   >>> order of dimensions correct for variable "{0}"'.format(varname[0]))
else:
    print('   >>> Dimensions found {0} does not match the dimension names specified with (-d): {1}'.format(lon_dims,dimname))
    raise ValueError('STOP')

if np.all(np.array(lat_dims) == dimname[::1]):
    lat = np.transpose(lat)
    print('   >>> switched order of dimensions for variable "{0}"'.format(varname[1]))
elif np.all(np.array(lat_dims) == dimname[::-1]):
    print('   >>> order of dimensions correct for variable "{0}"'.format(varname[1]))
else:
    print('   >>> Dimensions found {0} does not match the dimension names specified with (-d): {1}'.format(lat_dims,dimname))
    raise ValueError('STOP')

lath, lonh    = create_gridcells_from_centers(lat, lon)

nlon       = np.shape(lon)[1]
nlat       = np.shape(lat)[0]


# -------------------------------
# Read Basin shapes and all subbasin-shapes
# -------------------------------
print('   (2) Reading routing toolbox data ...')

shape     = gpd.read_file(routinginfo)
# shape     = shape.to_crs(epsg=crs_lldeg)        # this is lat/lon in degree
shape     = shape.to_crs(epsg=crs_caea)           # WGS 84 / North Pole LAEA Canada


# select only relevant basins/sub-basins

if not(basin is None):    # if gauge ID is given
    
    idx_basin = list(np.where(shape['Obs_NM']==basin)[0]) # list(np.where(shape['Gauge_ID']==basin)[0])
    # print('   >>> subbasins found = ',list(np.sort(np.array(shape.loc[idx_basin].HRU_ID,dtype=np.int))),'  (total: ',len(idx_basin),')')

    # find corresponding SubId
    SubId = np.int(shape.loc[idx_basin].SubId)
    print("   >>> found gauge at SubId = ",SubId)

if not(SubId is None): # if routing toolbox basin ID is given

    # maybe replace "SubId" with "HRU_ID" -->> needs to be tested

    old_SubId     = []
    new_SubId     = [ SubId ]

    while len(new_SubId) > 0:
        
        old_SubId.append(new_SubId)
        new_SubId = [ list(shape.loc[(np.where(shape['DowSubId']==ii))[0]].SubId) for ii in new_SubId ]  # find all upstream catchments of these new basins
        new_SubId = list(np.unique([item for sublist in new_SubId for item in sublist])) # flatten list and make entries unique

        # print("new_SubId = ",new_SubId)

    old_SubId = np.array([item for sublist in old_SubId for item in sublist],dtype=np.int)  # flatten list
    
    idx_basin = [ list(np.where(shape['SubId']==oo)[0]) for oo in old_SubId ]
    idx_basin = [ item for sublist in idx_basin for item in sublist ]  # flatten list
    idx_basin = list(np.unique(np.sort(idx_basin)))                    # getting only unique list indexes

    # get HRU_ID associated to those SubIds (SubIds can appear multiple times- one for "lake" and one for "land")
    
    print('   >>> HRU_IDs found = ',list(shape.loc[idx_basin].HRU_ID),'  (total: ',len(idx_basin),')')


shape     = shape.loc[idx_basin]

# indexes of all lines in df
keys = shape.index
nsubbasins = len(keys) 

# initialize
coord_catch_wkt = {}

# loop over all subbasins and transform coordinates into equal-area projection
for kk in keys:

    ibasin = shape.loc[kk]

    poly      = ibasin.geometry
    coord_catch_wkt[kk] = ogr.CreateGeometryFromWkt(poly.to_wkt())

# -------------------------------
# construct all grid cell polygons
# -------------------------------
print('   (3) Generate shapes for NetCDF grid cells ...')

grid_cell_geom_gpd_wkt = [ [ [] for ilon in range(nlon) ] for ilat in range(nlat) ]
for ilat in range(nlat):
    if ilat%10 == 0:
        print('   >>> Latitudes done: {0} of {1}'.format(ilat,nlat))
            
    for ilon in range(nlon):

        # -------------------------
        # EPSG:3035   needs a swap before and after transform ...
        # -------------------------
        # gridcell_edges = [ [lath[ilat,ilon]    , lonh[ilat,  ilon]    ],            # for some reason need to switch lat/lon that transform works
        #                    [lath[ilat+1,ilon]  , lonh[ilat+1,ilon]    ],
        #                    [lath[ilat+1,ilon+1], lonh[ilat+1,ilon+1]  ],
        #                    [lath[ilat,ilon+1]  , lonh[ilat,  ilon+1]  ]]

        # tmp = shape_to_geometry(gridcell_edges, epsg=crs_caea)
        # tmp.SwapXY()              # switch lat/lon back
        # grid_cell_geom_gpd_wkt[ilat][ilon] = tmp

        # -------------------------
        # EPSG:3573   does not need a swap after transform ... and is much faster than transform with EPSG:3035
        # -------------------------
        gridcell_edges = [ [lath[ilat,ilon]    , lonh[ilat,  ilon]    ],            # for some reason need to switch lat/lon that transform works
                           [lath[ilat+1,ilon]  , lonh[ilat+1,ilon]    ],
                           [lath[ilat+1,ilon+1], lonh[ilat+1,ilon+1]  ],
                           [lath[ilat,ilon+1]  , lonh[ilat,  ilon+1]  ]]

        tmp = shape_to_geometry(gridcell_edges, epsg=crs_caea)
        grid_cell_geom_gpd_wkt[ilat][ilon] = tmp

        
# -------------------------------
# Derive overlay and calculate weights
# -------------------------------
print('   (4) Deriving weights ...')

filename = output_file
ff       = open(filename,'w')
ff.write(':GridWeights                     \n')
ff.write('   #                                \n')
ff.write('   # [# HRUs]                       \n')
ff.write('   :NumberHRUs       {0}            \n'.format(nsubbasins))
ff.write('   :NumberGridCells  {0}            \n'.format(nlon*nlat))
ff.write('   #                                \n')
ff.write('   # [HRU ID] [Cell #] [w_kl]       \n')
            

for ikk,kk in enumerate(keys):

    ibasin = shape.loc[kk]

    area_basin = coord_catch_wkt[kk].Area()
    enve_basin = coord_catch_wkt[kk].GetEnvelope()   # bounding box around basin (for easy check of proximity)

    area_all = 0.0
    ncells   = 0

    for ilat in range(nlat):
        for ilon in range(nlon):
            
            enve_gridcell  = grid_cell_geom_gpd_wkt[ilat][ilon].GetEnvelope()   # bounding box around grid-cell (for easy check of proximity)
            grid_is_close  = check_proximity_of_envelops(enve_gridcell, enve_basin)

            if grid_is_close: # this check decreases runtime DRASTICALLY (from ~6h to ~1min)

                grid_cell_area = grid_cell_geom_gpd_wkt[ilat][ilon].Area()
                
                inter = grid_cell_geom_gpd_wkt[ilat][ilon].Intersection(coord_catch_wkt[kk].Buffer(0.0)) # "fake" buffer to avoid invalid polygons and weirdos dumped by ArcGIS
                area_intersect = inter.Area()
                
                area_all += area_intersect
                if area_intersect > 0:
                    ncells += 1

                    print("   >>> {0},{1},{2},{3},{4}".format(int(ibasin.HRU_ID),ilat,ilon,ilat*nlon+ilon,area_intersect/area_basin))
                    ff.write("   {0}   {1}   {2}\n".format(int(ibasin.HRU_ID),ilat*nlon+ilon,area_intersect/area_basin))

    print('   >>> (Sub-)Basin: {0} ({1} of {2})'.format(int(ibasin.HRU_ID),ikk+1,nsubbasins))
    print('   >>> Derived area of {0}  cells: {1}'.format(ncells,area_all))
    print('   >>> Read area from shapefile:   {0}'.format(area_basin))  
    print('   >>> error:                      {0}%'.format((area_basin - area_all)/area_basin*100.))
    print('   ')


ff.write(':EndGridWeights \n')
ff.close()

print('')
print('Wrote: ',filename)
print('')









