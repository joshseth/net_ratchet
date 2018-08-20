#!/bin/bash
#SBATCH -p long
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --mem-per-cpu=4G

module load R

set -u

echo "Git hash: $(git rev-parse HEAD)"
pwd

echo "Running:"
echo "   ./evolve_sexual.R $BASEDIR $POPSIZE $MAXGEN $PMUT $SMUT $NCORES"

/usr/bin/time --format='elapsed: %E / kernel: %S / user: %U / mem: %M' \
    ./evolve_sexual.R $BASEDIR $POPSIZE $MAXGEN $PMUT $SMUT $NCORES

echo "Done!"


