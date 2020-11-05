# PYTHON 3.6 

```bash
cd ~
mkdir env-3.6

module load mpi4py
module load proj
module load python/3.6

virtualenv --no-download ~/env-3.6
source ~/env-3.6/bin/activate

pip install --no-index --upgrade pip
pip install numpy --no-index
pip install GDAL --no-index
pip install argparse --no-index
pip install geopandas --no-index
pip install netCDF4 --no-index

# change to the folder containing the script "derive_grid_weights.py"
python derive_grid_weights.py -i example/input_MESH/RFF_H_GRD.nc -d "rlon,rlat" -v "longitude,latitude" -r example/maps/HRUs_coarse.shp -b 02LE024 -o example/input_MESH/GridWeights_RFF_H_GRD.txt
```


# PYTHON 3.8

```bash
cd ~
mkdir env-3.8

module load StdEnv/2020
module load gcc/9.3.0 gdal/3.0.4
module load mpi4py
module load proj
module load python/3.8

virtualenv --no-download ~/env-3.8
source ~/env-3.8/bin/activate

pip install --no-index --upgrade pip
pip install numpy --no-index
pip install GDAL --no-index
pip install argparse --no-index
pip install geopandas --no-index
pip install netCDF4 --no-index

# change to the folder containing the script "derive_grid_weights.py"
python derive_grid_weights.py -i example/input_MESH/RFF_H_GRD.nc -d "rlon,rlat" -v "longitude,latitude" -r example/maps/HRUs_coarse.shp -b 02LE024 -o example/input_MESH/GridWeights_RFF_H_GRD.txt
```
