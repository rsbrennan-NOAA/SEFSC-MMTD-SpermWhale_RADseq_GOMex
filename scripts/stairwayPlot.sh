#!/bin/bash
#SBATCH --job-name=StairwayPlot
#SBATCH --mail-user=reid.brennan@noaa.gov
#SBATCH --mail-type=END
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH -D /home/rbrennan/Tursiops-RAD-popgen/logout
#SBATCH -n 1
#SBATCH --mem=32G
#SBATCH --partition=standard
#SBATCH --time=48:00:00

cd ~/Tursiops-RAD-popgen/analysis/moments/stairway_plot_v2.1.1

BPPrefix=fourpop_CGulf

echo "Population is $BPPrefix"

java -cp stairway_plot_es Stairbuilder ${BPPrefix}.blueprint

bash ${BPPrefix}.blueprint.sh 

