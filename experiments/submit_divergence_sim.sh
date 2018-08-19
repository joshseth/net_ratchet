#!/bin/bash

USAGE="
Usage:
    $0 (directory with system) (population_size) (max_generation) (p_mut) (sigma_mut) (ncores)
"

if [ $# -lt 6 ]
then
    echo "$USAGE"
    exit 0
fi

BASEDIR="$1"
export BASEDIR
shift
POPSIZE="$1"
export POPSIZE
shift
MAXGEN="$1"
export MAXGEN
shift
PMUT="$1"
export PMUT
shift
SMUT="$1"
export SMUT
shift
NCORES="$1"
export NCORES
shift


if [ ! -e $BASEDIR ]
then
    echo "$BASEDIR does not exist."
    exit 1
fi

TAG=$(printf "%06d" $RANDOM); 
RUNDIR=$(dirname "${BASH_SOURCE[0]}")

sbatch -o $BASEDIR/slurm_${TAG}.out -e $BASEDIR/slurm_${TAG}.out \
    --ntasks-per-core=$NCORES \
    $RUNDIR/run_divergence_sim.sbatch 

