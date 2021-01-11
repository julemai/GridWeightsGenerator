#!/bin/bash

# submit with:
#       sbatch submit-to-graham.sh    

#SBATCH --account=<your-account>                   # your group rpp-kshook
#SBATCH --mem-per-cpu=30GB                         # memory; default unit is megabytes
#SBATCH --mail-user=<your-email>                   # email address for notifications
#SBATCH --mail-type=FAIL                           # email send only in case of failure
#SBATCH --time=0-01:00:00                          # time (DD-HH:MM:SS); 
#SBATCH --job-name=test-grid-weights               # name of job in queque


module load netcdf
module load StdEnv/2020
module load gcc/9.3.0 gdal/3.0.4
module load mpi4py
module load proj
module load python/3.8

# change PATH to your Python environment
source /home/julemai/env-3.8/bin/activate     

# chnage path to your cloned GitHub repo
cd /project/6008034/julemai/GridWeightsGenerator/

# chnage if you want to run other example
python derive_grid_weights.py -i example/input_NOAH-MP/Roff_small.nc -d "x,y" -v "LONGITUDE,LATITUDE" -r example/maps/HRUs_coarse.shp -b 02LE024 -o example/input_NOAH-MP/GridWeights_Roff_small.txt

