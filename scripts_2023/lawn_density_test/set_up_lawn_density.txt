# Set up to run lawn_density.py

# for this github example, gunzip the scanner images. the uncompressed files are each about 45Mb so the files get big quickly

# change this to your path -->
cd path/to/your/lawn_density_test/example_images

for file in *.tif.gz; do gunzip "$file"; done

# go back to your main directory
cd

##############
# If you're using a cluster load conda, if not, disregard this step
module load conda

# create a new conda environment 
conda create --name lawn_density_env python=3.9.15

# load the environment you just made
source activate lawn_density_env

# install the dependencies
conda install -c conda-forge opencv

conda install numpy

conda install matplotlib

##############

# write a slurm script to run the python script. If you are running this on a local machine just activate your 
# conda environment and run the python script (assuming the python script is in the enclosing folder).

# lawn_density.py has 3 arguments that are in the following order:
#
# directory --> directory with the scanner experiment images
# endpoint_dir --> directory where the outputs will go (see description of outputs below)
# min_radius --> minimum radius of the petri dish in pixels. for this experiment 662 worked
# well for 100 mm plates. You may need to adjust this for your experiment. 
# max_radius --> maximum radius of the petri dish in pixels. For this experiment 676 worked for 100mm plates. 


#!/bin/bash -l
#SBATCH --time=1:00:00
#SBATCH --ntasks=8
#SBATCH --mem=16g
#SBATCH --tmp=16g

module load conda 

source activate lawn_density_env

# change this to your directory --> 
cd path/to/your/lawn_density_test

python lawn_density.py example_images example_results 662 676


# assuming you are using slurm scheduler on a cluster, you'd save the script above as something like run_lawn_density.slurm in a relevant directory
# and then run the slurm script like this (without a hash, of course):
# sbatch path/to/run_lawn_density.slurm

##

# NOTE: There are 3 outputs:
# 
# 1) image_with_rectangles.tsv --> this is a photo with labeled circles where the 
# the program identified petri dishes locations. The labels are not in any particular order, 
# so you must know the locations of your plates on the scanner and then analyze your samples with
# with the relevant plate locations:sample ids.
#
# 2) density_data.csv --> mean grey values across the plate at a particular timepoint
#
# 3) variation_data.csv --> sd of grey values across the plate at a particular timepoint -- this 
# could be useful if the plate has substantial variation in growth across the plate (spatial 
# heterogeneity).

# Depending on your experiment, you may need to analyze the results differently from the way 
# I did, but a good starting place may be to look at:
# https://github.com/JonMartinson/ecology_DFE/blob/main/scripts_2023/09_agar_scanner_growth.R


# to avoid any confusion this is an example file structure
# needed to run this program

lawn_density_test/
|-- lawn_density.py
|-- example_images
|-- example_results