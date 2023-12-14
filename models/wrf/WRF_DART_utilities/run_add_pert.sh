#!/bin/sh

# test the program that addds perturbations to the U, V, T and QV fields where
# reflectivity is high.  if all state variables are 0, the assimilation will not
# be able to change the values even if the observations show it should.
#
# this test program requires a wrfinput_d01 that will be modified.
# the rest of the inputs can be used unchanged for testing.
#
#
# the inputs needed to run (from a comment in the file):
#
#
# input parameters from command line:
# (1) refl_ob_file  -- name of text file containing WRF grid indices where observed reflectivity is high
# (2) wrf_file      -- path name of WRF netcdf file
# (3) lh            -- horizontal length scale (m) for perturbations
# (4) lv            -- vertical length scale (m) for perturbations
# (5) u_sd          -- std. dev. of u noise (m/s), before smoothing
# (6) v_sd          -- std. dev. of v noise (m/s), before smoothing
# (7) w_sd          -- std. dev. of w noise (m/s), before smoothing
# (8) t_sd          -- std. dev. of potential temperature noise (K), before smoothing
# (9) td_sd         -- std. dev. of dewpoint noise (K), before smoothing
# (10) qv_sd        -- std. dev. of water vapor mixing ratio noise, before smoothing
#                         (input value is in g/kg, value after conversion is in kg/kg)
# (11) ens_num      -- ensemble number for consistently seeding the random number generator
# (12) gdays        -- gregorian day number of wrf analysis time
# (13) gsecs        -- gregorian seconds of day of wrf analysis time
#
#
# 
# need to adjust last 2 parms to match the date of the given wrfinput_d01 file
# for now, use jan 1, 2022 midnight:
#
# 153767     0 == 2022/01/01 00:00:00
#

cp -f wrfinput_d01 wrfinput_d01.orig

ln -sf ../work/input.nml .
ln -sf ../work/add_pert_where_high_refl .

./add_pert_where_high_refl add_pert_test_input.txt wrfinput_d01 100 100 0.01 0.01 0.01 0.1 0.01 0.01  1  153767 0

ncdiff -O wrfinput_d01 wrfinput_d01.orig pert_diffs.nc
echo look at U,V,T and QV in pert_diffs.nc to see differences added by this program

exit 0

