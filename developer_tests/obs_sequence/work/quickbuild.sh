#!/usr/bin/env bash -v

# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

main() {


export DART=$(git rev-parse --show-toplevel)
source "$DART"/build_templates/buildfunctions.sh

MODEL="lorenz_96"
EXTRA="$DART"/assimilation_code/programs/obs_sequence_tool

echo $EXTRA
ls -l $EXTRA

LOCATION="oned"
dev_test=1
TEST="obs_sequence"


serial_programs=(
obs_sequence_tool
)

echo $serial_programs

# quickbuild arguments
arguments "$@"

# clean the directory
\rm -f -- *.o *.mod Makefile .cppdefs

# build and run preprocess before making any other DART executables
buildpreprocess

# build 
buildit

# clean up
\rm -f -- *.o *.mod

}

main "$@"
