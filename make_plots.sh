#!/bin/bash

USAGE="
Runs descriptive scripts on a system:
    ./make_plots.sh (directory with system) (max number of extra dimensions) (sigma)
"

if [ $# -ne 3 ]
then
    echo "$USAGE"
    exit 1
fi

SYSDIR=$1
shift
MAXDIM=$1
shift
SIGMA=$1

set -eu

if [ ! -d $SYSDIR -o ! -f $SYSDIR/params.R ]
then
    echo "System does not exist in $SYSDIR."
    exit 1
fi

for K in $(seq $MAXDIM)
do
    ./do_kryptotype_plots.R $SYSDIR $K $SIGMA
    ./do_simultaneous_mutations.R $SYSDIR $K $SIGMA
    ./do_gene_deletions.R $SYSDIR $K $SIGMA
done

pdfjoin --outfile $SYSDIR/kryptotype_plots.pdf $SYSDIR/kryptotype_plots/*.pdf &>/dev/null \
    && rm -rf $SYSDIR/kryptotype_plots
pdfjoin --outfile $SYSDIR/simultaneous_mutations_phenotypes.pdf $SYSDIR/simultaneous_mutations_phenotypes/*.pdf &>/dev/null \
    && rm -rf $SYSDIR/simultaneous_mutations_phenotypes
pdfjoin --outfile $SYSDIR/deletion_phenotypes.pdf $SYSDIR/deletion_phenotypes/*.pdf &>/dev/null \
    && rm -rf $SYSDIR/deletion_phenotypes

echo $SYSDIR

