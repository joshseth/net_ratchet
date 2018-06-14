#!/bin/bash

# Make a random system and run our scripts on it,
# to test that the scripts run

set -eu

SYSDIR=$(./generate_system.R 4 1 1)

if [ ! -d $SYSDIR -o ! -f $SYSDIR/params.R ]
then
    echo "Failed to generate a system in $SYSDIR."
    exit 1
fi

./make_plots.sh $SYSDIR 4 0.1 3

EVOLDIR=$(./evolve_population.R $SYSDIR 20 10 0.1 0.01 0.05 0.05)
echo $EVOLDIR

echo "Done!"
