#!/bin/bash

USAGE="
Runs descriptive scripts on a system:
    ./make_plots.sh (directory with system) (max number of extra dimensions) (sigma) (number of replicates)
"

if [ $# -ne 4 ]
then
    echo "$USAGE"
    exit 1
fi

SYSDIR=$1
shift
MAXDIM=$1
shift
SIGMA=$1
shift
NREPS=$1

set -eu

if [ ! -d $SYSDIR -o ! -f $SYSDIR/params.R ]
then
    echo "System does not exist in $SYSDIR."
    exit 1
fi

for K in $(seq $MAXDIM)
do
    ./do_plots.R $SYSDIR $K $SIGMA $NREPS
done

pdfjoin --outfile $SYSDIR/kryptotypes.pdf $SYSDIR/kryptotypes/*.pdf &>/dev/null \
    && rm -rf $SYSDIR/kryptotypes
pdfjoin --outfile $SYSDIR/simultaneous_mutations.pdf $SYSDIR/simultaneous_mutations/*.pdf &>/dev/null \
    && rm -rf $SYSDIR/simultaneous_mutations
pdfjoin --outfile $SYSDIR/gene_deletions.pdf $SYSDIR/gene_deletions/*.pdf &>/dev/null \
    && rm -rf $SYSDIR/gene_deletions
pdfjoin --outfile $SYSDIR/eigenvalues.pdf $SYSDIR/eigenvalues/*.pdf &>/dev/null \
    && rm -rf $SYSDIR/eigenvalues

echo $SYSDIR

