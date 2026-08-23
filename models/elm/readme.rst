ELM
===

|ELM gridcell breakdown|
ELM retains the gridcell, landunit, column, and plant-functional-type hierarchy
shown in this figure from the CLM5.0 Technical Note. The image links to the
original description.

Overview
--------

This is the DART interface to the E3SM Land Model (ELM). The maintained target is
E3SM ``maint-3.0``; the interface, SourceMods, and workflow in this directory were
checked against E3SM revision ``1eda4d4dcacd9ff603764264d5bf8ca861b97ac1``.

E3SM/CIME advances the ensemble and writes restart files at each assimilation
time. DART reads those files, performs the land analysis, and writes updated ELM
restart states for the next coupled forecast. Become familiar with creating and
running an E3SM multi-instance case before configuring data assimilation, then
follow `Configuring an E3SMv3 Experiment`_ below.

The bundled ``tutorial`` directory is retained as legacy CLM5/CESM material and
has not been ported to this ELM/E3SMv3 workflow. Do not use its scripts to create
an E3SMv3 experiment.

Land observations and forward operators are defined by DART modules such as
``obs_def_land_mod.f90``, ``obs_def_tower_mod.f90``, and
``obs_def_COSMOS_mod.f90``. ELM-DART applications have included snow, soil
moisture, leaf area index, biomass, carbon flux, and related observations. The
historical CLM-DART publications in `References`_ remain useful scientific
background for this descendant interface.


Important Features
------------------

Land DA is extremely diverse. The support for Land DA as pertains to ELM-DART
has some features that need to be described in some detail.

SourceMods
~~~~~~~~~~

The version-matched ELM overrides are in
``DART_SourceMods/e3sm_maint_3.0/src.elm``. They were verified against E3SM
revision ``1eda4d4dcacd9ff603764264d5bf8ca861b97ac1``. CIME requires the
Fortran files directly under ``SourceMods/src.elm``; do not reproduce the ELM
source-tree subdirectories.

``EcosystemBalanceCheckMod.F90`` suppresses expected C/N/P balance failures on
the first restart step after a DART update. ``SurfaceRadiationMod.F90`` exposes
``PARVEG`` for the absorbed-PAR forward operator. The E3SM setup workflow
installs these files before ``case.setup`` and compilation. See the versioned
SourceMods README for limitations, including the SIF files that are not yet
ported.

ELM indeterminate values
~~~~~~~~~~~~~~~~~~~~~~~~

ELM variables are encoded as rectangular arrays in the netCDF files.
However, this means that some variables have space for layers that are unused.
Anything with snow layers, for example. ELM has the *SNLSNO* variable to indicate
which snow layers are active. The *unused* layers may not have the *_FillValue*
value, but can have 'indeterminate' values. The :doc:`elm_to_dart`
must be run to convert these indeterminate values to *_FillValue* to be
interpreted correctly by DART.  After the assimilation is complete, the
:doc:`dart_to_elm` must be called to replace the *_FillValue* with whatever
is originally in that slot. This approach preserves the 'indeterminate' value
for *unused* snow layers and prevents DART from adjusting the value during
the *filter* step. If the surface snow layer has a *trace* of snow this is
considered an active snow layer, and we allow DART to adjust this value.
See the *Discussion of Indeterminate Values*
section of :doc:`elm_to_dart` for more details.


Model Interpolate - The Forward Operator
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Since the subgridscale components of ELM have no explicit location associated
with them, the location of every component in the gridcell is the same as the
gridcell itself. The DART forward operators fundamentally rely on
interpolating the model state to some arbitrary location. At present, the best we
can do is to create an area-weighted average of all components in the gridcell.
This is sub-optimal because it introduces representation mismatch between the
grid cell and observation spatial resolution. A nice project would be to use a lookup
table for the observation location to determine the dominant PFT (or relevant metadata) at
that location and only average the PFTs specifically associated with
the observation within the gridcell. This will allow the forward operator to be
more accurate and might have a discernable impact on the regression relationship
(i.e. ensemble covariance) between the variables in the DART state vector.

The *model_interpolate* function in DART achieves efficiency by interpolating
all the ensemble members at the same time. This gives rise to some challenging
problems when interpolating values for variables with changing numbers of active layers.
For example, some ensemble members may only have 2 active snow layers, some may have 3.
This is an untenable situation when asked for the snow temperature or water
content in layer 3, for example. Consequently - *model_interpolate* will fail
and return an error code - the forward operator will fail - and the observation
is rejected and the DART QC is marked as such. Be aware.

Localization
~~~~~~~~~~~~

Localization is the term used to restrict the portion of the state to the portion
believed to be related to the observation. Most often, this is a spatial argument
but it does not need to be restricted to that. In some way, even the selection
of the ELM variables to include in the DART state is a de-facto localization.
Since ELM has such a rich description of land unit types: urban columns, glaciers,
lakes, etc. it is also possible (and probably desirable) to explicitly declare
some columns and/or PFTs to be unaffected by the assimilation - i.e., we
declare that soil moisture observations should not impact urban columns
or deep lakes or ... The **get_close_state()** function employs a routine to
explicitly declare what subgridscale components are allowed to be modified by
the assimilation. This routine can easily be customized to suit your purpose.
The code segment below should make this clear.

.. code-block:: fortran

  ! Determine if state_index is a variable from a column (or whatever is of interest).
  ! Determine what dimension is of interest, need to know to index into
  ! cols1d_ityplun(ncolumn) array (for example).

  RELATEDLOOP: do jdim = 1, get_num_dims(dom_id, var_id)

     dimension_name = get_dim_name(dom_id, var_id, jdim)
     select case ( trim(dimension_name) )
            case ("gridcell","lon","lat")
               related = .true.
            case ("lndgrid")
               related = .true.
            case ("landunit")
               if ( land1d_ityplun(indices(jdim)) == ilun_vegetated_or_bare_soil ) related = .true.
               if ( land1d_ityplun(indices(jdim)) == ilun_crop                   ) related = .true.
            case ("column")
               if ( cols1d_ityplun(indices(jdim)) == icol_vegetated_or_bare_soil ) related = .true.
               if ( cols1d_ityplun(indices(jdim)) == icol_crop                   ) related = .true.
            case ("pft")
               related = .true.
            case default
     end select

     ! Since variables can use only one of these dimensions,
     ! there is no need to check the other dimensions.
     if (related) exit RELATEDLOOP

  enddo RELATEDLOOP


Snow Data Assimilation
~~~~~~~~~~~~~~~~~~~~~~

The *prognostic* variables for snow (i.e. the ones that impact the forecast)
are the ones that have layers. The snow observations are typically without
explicit depths and are essentially column-integrated quantities like snow
water equivalent (SWE - ELM variable *H2OSNO*) or snow depth
(ELM variable *SNOW_DEPTH*).  These ELM *diagnostic* variables
simplify the forward operator and have been part of the DART state.
However, updating the prognostic variables
(*T_SOISNO, H2OSOI_LIQ, H2OSOI_ICE, DZSNO, ZSNO, ZISNO*)
through their ensemble covariance with the update to *H2OSNO*, for example
will generally not result in a posterior SWE (calculated from the prognostic
variables) that matches the posterior SWE in *H2OSNO*.

**In order to address this challenge,** a snow repartition function has been
created in ``dart_to_elm`` that redistributes the posterior SWE into
the appropriate prognostic variables. This **guarantees** that the posterior
SWE of the prognostic snow variables matches the posterior SWE in H2OSNO.
**When snow related variables are being updated within an assimilation it is
recommended to invoke this repartitioning function by setting the namelist
option ``repartition_swe = 1 or 2`` within ``&dart_to_elm_nml``**. See the
:doc:`dart_to_elm` for more details describing the repartitioning function including
guidance on how to set up a case that repartitions snow. Note that we have not
attempted to include any of the snow property variables most important to controlling
albedo (eg. grain radius, carbon, dust) within the DART state.  To what extent adjusting
mass and dimensional properties of snow layers indirectly influences the
albedo properties is an active scientific question. See the :doc:`dart_to_elm`
for more details on how to implement ``repartition_swe`` if conserving albedo
is important for your application.

The snow formulation in ELM is complex. Reducing the amount of snow through
assimilation is well-defined. Creating snow when there is none is
**a limited capability** in ELM-DART. If snow exists for a subset of ensemble
members at a given location, then it is possible to adjust ensemble members
with a value of zero to a non-zero value.  On the other hand,
**if all ensemble members do not have snow, or at least one member has a FillValue**,
the statistical assumptions for ensemble data assimilation are
not valid and the snow variables remain at zero. The best method would be to alter the
amount of snow *from the forcing file* and let ELM manage the snow. This is
beyond the scope of ELM-DART. We have thought that if one member does not have
snow - maybe we should just use the values from some other member - but when
does that stop being acceptable? 10 ensemble members? 20? The distributions
become multimodal, and the logical end result is that you could wind up using
1 ensemble member to declare the snow for all the remaining members. That seems
like a bad idea.

Similar logic applies to the variables related to plant growth. If the LAI
observations indicates there should be something growing and nothing has
sprouted yet, DART does nothing to the variables.


Configuring an E3SMv3 Experiment
---------------------------------

The maintained case and cycling workflow is in
``shell_scripts/e3sm_maint_3.0``. It was derived from the operational
``v3_dart_cda/3_ne30pg2_dart_cpl_en40`` experiment and supports coupled E3SMv3
forecasts with independently configurable EAM and ELM DART analyses.

Read ``shell_scripts/e3sm_maint_3.0/README.md`` before use. The main entry
points are:

+------------------------------------+-----------------------------------------------------------+
| Script                             | Purpose                                                   |
+====================================+===========================================================+
| ``create_and_setup_case.sh``       | Shared E3SM, ensemble, restart, DART, and data settings.  |
+------------------------------------+-----------------------------------------------------------+
| ``1_run_e3sm_ensemble_setup.sh``   | Build the base E3SM case and create ensemble cases.       |
+------------------------------------+-----------------------------------------------------------+
| ``2_run_dart_e3sm_icbc.sh``        | Prepare and validate coupled initial conditions.          |
+------------------------------------+-----------------------------------------------------------+
| ``3_run_dart_eam_perturb.sh``      | Create the initial EAM ensemble perturbations.            |
+------------------------------------+-----------------------------------------------------------+
| ``4_run_dart_e3sm_cycleda.sh``     | Run restart-safe coupled forecast-assimilation cycles.    |
+------------------------------------+-----------------------------------------------------------+
| ``5`` through ``8``                | Optional compression, diagnostics, and post-processing.   |
+------------------------------------+-----------------------------------------------------------+

Step 1 installs the version-matched ELM files from
``DART_SourceMods/e3sm_maint_3.0/src.elm`` before ``case.setup`` and
compilation. Review all paths in ``create_and_setup_case.sh`` and all Slurm
``#SBATCH`` directives before submitting the numbered scripts.

Declaring the Variables in the DART State
-----------------------------------------

The DART state vector is constructed in a very flexible manner.
A namelist is used to relate the netCDF variable name, the netCDF file
type [restart, (XY) history, or vector history] with a DART QUANTITY.
Including variables from an 'XY' ELM history file allows the
inclusion of diagnostic variables that can speed up the forward
observation operators if gridcell averages are appropriate.

It is also possible to read some variables from the restart file,
and some from a 'vector-based' history file that has the same
structure (gridcell/landunit/column/pft) as the restart file - but may be
temporal averages instead of instantaneous quantities.
Care must be taken to assign the proper DART QUANTITY to the variables.
Any variable in the DART state is updated, but the forward operator
looks for specific QUANTITIES. If you want to use the vector-based history
file for the forward operator - make sure you declare it to be of the
QUANTITY used by the forward operator code.

.. "Simple" observations like snowcover fraction come directly from
   the DART state. It is possible to configure the ELM history files
   to contain the ELM estimates of some quantities (mostly flux tower
   observations e.g, net ecosystem production, sensible heat flux,
   latent heat flux) that are very complicated combinations of portions
   of the ELM state.  The forward observation operators for these flux tower
   observations read these quantities from the ELM ``.h1.`` history file.
   The smaller the ELM gridcell, the more likely it seems that these
   values will agree with point observations. Be advised that the
   **obs_def_tower_mod.f90** is **not supported in this version**.

The namelist specification of what goes into the DART state vector
includes the ability to specify if the quantity should have a lower
bound, upper bound, or both, what file the variable should be read
from, and if the variable should be modified by the assimilation or not.
Make sure you read the `Inflation`_ section to fully understand what
happens when you designate a variable 'NO_COPY_BACK'.

.. attention::

   It is important to know that the variables in the DART diagnostic files
   ``preassim``, ``postassim``, ``analysis``, and ``output`` will contain
   the unbounded versions of ALL the variables specified in ``elm_variables``.
   Only the files specified in the ``filter_nml:output_state_file_list``
   will have the 'clamped' values.

The example ``input.nml`` ``model_nml`` demonstrates how to construct the
DART state vector. The following table explains in detail each entry
for ``elm_variables``:

.. container::

   ======== ==============================================================
    Column  Description
   ======== ==============================================================
    **1**   The ELM variable name as it appears in the ELM netCDF file.
    **2**   The corresponding DART QUANTITY.
    **3**   | Minimum value of the posterior.
            | If set to 'NA' there is no minimum value.
            | The DART diagnostic files will not reflect this value, but
            | the file used to restart ELM will.
    **4**   | Maximum value of the posterior.
            | If set to 'NA' there is no maximum value.
            | The DART diagnostic files will not reflect this value, but
            | the file used to restart ELM will.
    **5**   | Specifies which file should be used to obtain the variable.
            | ``'restart'`` => elm_restart_filename
            | ``'history'`` => elm_history_filename
            | ``'vector'``  => elm_vector_history_filename
    **6**   | Should ``filter`` update the variable in the specified file.
            | ``'UPDATE'`` => the variable is updated.
            | ``'NO_COPY_BACK'`` => the variable remains unchanged.
   ======== ==============================================================

The following are only meant to be examples - they are not scientifically validated.
Some of these that are UPDATED are probably diagnostic quantities, Some of these that
should be updated may be marked NO_COPY_BACK.  This list is by no means complete.

::

   elm_variables  = 'leafc',       'QTY_LEAF_CARBON',            '0.0', 'NA', 'restart' , 'UPDATE',
                    'frac_sno',    'QTY_SNOWCOVER_FRAC',         '0.0', '1.', 'restart' , 'UPDATE',
                    'SNOW_DEPTH',  'QTY_SNOW_THICKNESS',         '0.0', 'NA', 'restart' , 'NO_COPY_BACK',
                    'H2OSOI_LIQ',  'QTY_SOIL_LIQUID_WATER',      '0.0', 'NA', 'restart' , 'UPDATE',
                    'H2OSOI_ICE',  'QTY_SOIL_ICE',               '0.0', 'NA', 'restart' , 'UPDATE',
                    'T_SOISNO',    'QTY_TEMPERATURE',            '0.0', 'NA', 'restart' , 'UPDATE',
                    'livestemc',   'QTY_LIVE_STEM_CARBON',       '0.0', 'NA', 'restart' , 'UPDATE',
                    'deadstemc',   'QTY_DEAD_STEM_CARBON',       '0.0', 'NA', 'restart' , 'UPDATE',
                    'NEP',         'QTY_NET_CARBON_PRODUCTION',  'NA' , 'NA', 'history' , 'NO_COPY_BACK',
                    'H2OSOI',      'QTY_SOIL_MOISTURE',          '0.0', 'NA', 'history' , 'NO_COPY_BACK',
                    'SMINN_vr',    'QTY_SOIL_MINERAL_NITROGEN',  '0.0', 'NA', 'history' , 'NO_COPY_BACK',
                    'LITR1N_vr',   'QTY_NITROGEN',               '0.0', 'NA', 'history' , 'NO_COPY_BACK',
                    'TSOI',        'QTY_SOIL_TEMPERATURE',       'NA' , 'NA', 'history' , 'NO_COPY_BACK',
                    'FSDSVDLN',    'QTY_PAR_DIRECT',             '0.0', 'NA', 'history' , 'NO_COPY_BACK',
                    'FSDSVILN',    'QTY_PAR_DIFFUSE',            '0.0', 'NA', 'history' , 'NO_COPY_BACK',
                    'PARVEGLN',    'QTY_ABSORBED_PAR',           '0.0', 'NA', 'history' , 'NO_COPY_BACK',
                    'NEE',         'QTY_NET_CARBON_FLUX',        'NA' , 'NA', 'vector'  , 'NO_COPY_BACK',
                    'H2OSNO',      'QTY_SNOW_WATER',             '0.0', 'NA', 'vector'  , 'NO_COPY_BACK',
                    'TLAI',        'QTY_LEAF_AREA_INDEX',        '0.0', 'NA', 'vector'  , 'NO_COPY_BACK',
                    'TWS',         'QTY_TOTAL_WATER_STORAGE',    'NA' , 'NA', 'vector'  , 'NO_COPY_BACK',
                    'SOILC_vr',    'QTY_SOIL_CARBON',            '0.0', 'NA', 'vector'  , 'NO_COPY_BACK',
                    'SOIL1N_vr',   'QTY_SOIL_NITROGEN',          '0.0', 'NA', 'vector'  , 'NO_COPY_BACK',
                    'SMP',         'QTY_SOIL_MATRIC_POTENTIAL',  '0.0', 'NA', 'vector'  , 'NO_COPY_BACK'
      /


**Only the first variable for a DART QUANTITY in the elm_variables list will
be used for the forward observation operator.**
The following is perfectly legal:

::

   elm_variables = 'LAIP_VALUE', 'QTY_LEAF_AREA_INDEX', 'NA', 'NA', 'restart' , 'UPDATE',
                   'tlai',       'QTY_LEAF_AREA_INDEX', 'NA', 'NA', 'restart' , 'UPDATE',
                   'elai',       'QTY_LEAF_AREA_INDEX', 'NA', 'NA', 'restart' , 'UPDATE',
                   'ELAI',       'QTY_LEAF_AREA_INDEX', 'NA', 'NA', 'history' , 'NO_COPY_BACK',
                   'LAISHA',     'QTY_LEAF_AREA_INDEX', 'NA', 'NA', 'history' , 'NO_COPY_BACK',
                   'LAISUN',     'QTY_LEAF_AREA_INDEX', 'NA', 'NA', 'history' , 'NO_COPY_BACK',
                   'TLAI',       'QTY_LEAF_AREA_INDEX', 'NA', 'NA', 'history' , 'NO_COPY_BACK',
                   'TLAI',       'QTY_LEAF_AREA_INDEX', 'NA', 'NA', 'vector'  , 'NO_COPY_BACK'
      /

however, only **LAIP_VALUE** will be used to calculate the LAI when an
observation of LAI is encountered. **All** (the other LAI) variables in
the DART state will be modified by the assimilation based on the
relationship of LAIP_VALUE and the observation. It is possible that
several ELM variables could serve as the input for the forward operator,
however, in practice, the user should choose the variable that best
matches the observation (temporal/spatial resolution, units etc), to help
limit the complexity of the forward operator.

Inflation
---------

Inflation has been shown to be quite useful in our experience of
DA with ELM and DART. The model is strongly influenced by the
atmospheric forcing and will cause the ELM ensemble to
relax to a state consistent with the forcing when the assimilation
stops. Depending on the forecast length between assimilations, and
sometimes just to restore the variance lost during an assimilation,
inflation should be used.

The 'NO_COPY_BACK' designation has some side effects when it
comes to state-space inflation (inf_flavor 2,4 or 5 -
'VARYING_SS_INFLATION','RELAXATION_TO_PRIOR_SPREAD',
or 'ENHANCED_SS_INFLATION' - respectively).  State-space inflation
requires an inflation value for everything in the DART state.
If the variable has been designated as 'NO_COPY_BACK'
the DART write routine (when called from ``filter``) simply
skips the variable and nothing is written.
This is a problem for inflation files that need to adapt.

The solution is to run
:doc:`../../assimilation_code/programs/fill_inflation_restart/fill_inflation_restart`
to create an initial inflation file with inflation values of 1.0 (i.e.
no inflation). ``fill_inflation_restart`` has been specially designed
to output inflation values for every variable in the DART state.
The idea is to copy the *input* inflation file to the *output* inflation
file name *before each assimilation cycle*. No new values will be written
for the variables designated 'NO_COPY_BACK', the original values will persist.

It remains a scientific question as to whether or not this is the **right** thing
to do! The 'NO_COPY_BACK' mechanism was initially intended to simply avoid
writing variables that did not impact the next model forecast. Since inflation
is a powerful mechanism to overcome observation-model bias, it might be
perfectly warranted to 'UPDATE' these diagnostic variables. Be warned, if
you do 'UPDATE' the diagnostic variables, you may want to create copies
of the prior so you explore exactly what happens during an assimilation.

When land inflation is enabled, the maintained E3SM workflow initializes the
first-cycle files with ``fill_inflation_restart`` in
``workflow_lib/cycle/elm_dart_assimilation.sh`` and carries their names forward
between cycles.


.. attention::

   It is recommended to apply no inflation during the first assimilation step. In other
   words within ``input.nml`` and namelist ``&fill_inflation_restart_nml``
   set ``prior_inf_mean = 1.00`` and ``post_inf_mean = 1.00``.  Otherwise, a spatially
   uniform inflation will be applied to the entire spatial domain of the assimilation
   which can make ELM unstable. In general, inflation is intended to account for biases
   between the observation and model-estimated observation, as well as to restore ensemble
   spread after an observation has been assimilated.


Namelist
--------

Namelists start with an ampersand (``&``) and terminate with a slash (``/``).
Character strings containing a slash must be quoted. The ``model_mod.f90``
defaults are:

::

   &model_nml
      elm_restart_filename         = 'elm_restart.nc'
      elm_history_filename         = 'elm_history.nc'
      elm_vector_history_filename  = 'elm_vector_history.nc'
      assimilation_period_days     = 0
      assimilation_period_seconds  = 60
      calendar                     = 'Gregorian'
      debug                        = 0
      elm_variables                = ''
   /

The maintained workflow supplies an experiment configuration in
``shell_scripts/e3sm_maint_3.0/workflow_lib/namelists/elm/filter.nml``. Treat
that file as a starting point, not as a scientifically universal state vector.

``elm_restart_filename``
   ELM restart file used for sparse-grid metadata and restart variables. Only
   variables originating in this file may be copied back after assimilation.

``elm_history_filename``
   ELM history file used for the full grid, vertical-coordinate metadata, and
   any diagnostic variables whose ``elm_variables`` origin is ``history``.

``elm_vector_history_filename``
   Vector-form ELM history file used when an ``elm_variables`` origin is
   ``vector``.

``assimilation_period_days`` and ``assimilation_period_seconds``
   Together define the DART model time step and expected assimilation cadence.
   The workflow value must match the interval between ELM restart states.

``calendar``
   Calendar used to interpret ELM times. Use ``Gregorian`` for real-date
   observations unless the E3SM case intentionally uses another calendar.

``debug``
   Runtime diagnostic level. Zero produces minimal output; larger values enable
   progressively more interface diagnostics.

``elm_variables``
   A flat list of six strings per state variable:

   #. ELM netCDF variable name.
   #. DART quantity name generated in ``obs_kind_mod.f90`` by ``preprocess``.
   #. Posterior minimum, or ``NA`` for no lower bound.
   #. Posterior maximum, or ``NA`` for no upper bound.
   #. Origin file: ``restart``, ``history``, or ``vector``.
   #. Copy-back policy: ``UPDATE`` or ``NO_COPY_BACK``.

   All listed variables participate in DART internally. ``UPDATE`` writes a
   posterior variable back only when its origin is ``restart``;
   ``NO_COPY_BACK`` leaves the source file unchanged. Bounds are applied to
   files written by ``filter``, not to the unconstrained diagnostic copies.


Modules used
-----------------------------

::

   default_model_mod
   distributed_state_mod
   ensemble_manager_mod
   mpi_utilities_mod
   netcdf_utilities_mod
   obs_def_utilities_mod
   obs_kind_mod
   options_mod
   state_structure_mod
   threed_sphere/location_mod
   time_manager_mod
   types_mod
   utilities_mod


Files
-----

====================== ===========================================================================
filename               purpose
====================== ===========================================================================
input.nml              to read the model_mod namelist
elm_restart.nc         both read and modified by the ELM model_mod
elm_history.nc         read by the ELM model_mod for metadata and possible diagnostic variables.
elm_vector_history.nc  read by the ELM model_mod for possible diagnostic variables.
dart_log.out           the run-time diagnostic output
dart_log.nml           the record of all the namelists actually USED - contains the default values
====================== ===========================================================================


Error codes and conditions
--------------------------

+---------------------+---------------------------------------------+---------------------------------------------------+
|       Routine       |                   Message                   |                      Comment                      |
+=====================+=============================================+===================================================+
| nc_write_model_atts | Various netCDF-f90 interface error messages | From one of the netCDF calls in the named routine |
+---------------------+---------------------------------------------+---------------------------------------------------+


Future plans:
-------------

1. Implement a lookup table that relates the observation location to a dominant PFT or COLUMN
   so the *model_interpolate* code can average quantities from similar PFTs or COLUMNs instead
   of everything in the entire grid cell.
2. Implement a fast way to get the quantities needed for the calculation of
   radiative transfer models - needs a whole column of ELM variables, redundant if
   multiple frequencies are used.
3. Figure out what to do when one or more of the ensemble members does not have
   snow/leaves/etc. when the observation indicates there should be. Ditto for removing
   snow/leaves/etc. when the observation indicates otherwise.
4. Right now, the soil moisture observation operator is used by the COSMOS code to
   calculate the expected neutron intensity counts. This is the right idea, however,
   the COSMOS forward operator uses m3/m3 and the ELM units are kg/m2. I have not
   checked to see if they are, in fact, identical. This brings up a bigger issue in
   that the soil moisture observation operator would also be used to calculate whatever
   a TDT probe or ??? would measure. What units are they in? Can one operator support both?



References
----------

The
`CTSM Documentation <https://escomp.github.io/ctsm-docs/versions/master/html/index.html>`__
documents the CLM/CTSM heritage shared by ELM. The following historical CLM-DART publications provide scientific background:

       Zhang, Y.-F., T. J. Hoar, Z.-L. Yang, J. L. Anderson, A. M. Toure and M. Rodell, 2014:
       Assimilation of MODIS snow cover through the Data Assimilation Research Testbed
       and the Community Land Model version 4.
       *Journal of Geophysical Research: Atmospheres*, **142** 1489-1508,
       `doi:10.1002/2013JD021329 <https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2013JD021329>`__

       Lin, P., J. Wei, Z. -L. Yang, Y. Zhang, K. Zhang, 2016:
       Snow data assimilation‐constrained land initialization improves seasonal
       temperature prediction.
       *Geophysical Research Letters* **43** (21), 11,423-11,432
       `doi:10.1002/2016GL070966 <https://doi.org/10.1002/2016GL070966>`__

       Zhao, L., Z. -L. Yang and T. J. Hoar, 2016:
       Global soil moisture estimation by assimilating AMSR-E brightness temperatures
       in a coupled CLM4-RTM-DART system.
       *Journal of Hydrometeorology*, **17**, 2431-2454,
       `doi:10.1175/JHM-D-15-0218.1 <https://doi.org/10.1175/JHM-D-15-0218.1>`__

       Kwon, Y., Z. -L. Yang, T. J. Hoar and A. M. Toure, 2017:
       Improving the radiance assimilation performance in estimating snow water storage across
       snow and land-cover types in North America.
       *Journal of Hydrometeorology*, **18**, 651-668,
       `doi:10.1175/JHM-D-16-0102.1 <https://doi.org/10.1175/JHM-D-16-0102.1>`__

       Fox, A. M., Hoar, T. J., Anderson, J. L., Arellano, A. F., Smith, W. K., Litvak, M. E., et al., 2018:
       Evaluation of a data assimilation system for land surface models using CLM4.5.
       *Journal of Advances in Modeling Earth Systems*, **10**, 2471–2494,
       `doi.org/10.1029/2018MS001362 <https://doi.org/10.1029/2018MS001362>`__

       Ling, X. L., Fu, C. B., Yang, Z. L., & Guo, W. D., 2019:
       Comparison of different sequential assimilation algorithms for satellite-derived leaf area
       index using the Data Assimilation Research Testbed (version Lanai).
       *Geoscientific Model Development*, 12(7), 3119-3133.
       `doi.org/10.5194/gmd-12-3119-2019 <https://doi.org/10.5194/gmd-12-3119-2019>`__

       Bian, Q., Xu, Z., Zhao, L., Zhang, Y. F., Zheng, H., Shi, C., … & Yang, Z. L., 2019:
       Evaluation and intercomparison of multiple snow water equivalent products over the Tibetan Plateau.
       *Journal of Hydrometeorology*, 20(10), 2043-2055.
       `doi.org/10.1175/JHM-D-19-0011.1 <https://doi.org/10.1175/JHM-D-19-0011.1>`__

       Raczka, B., Hoar T.J., Duarte H.F., Fox A.M., Anderson J.L., Bowling D.R., & Lin J.C., 2021
       Improving CLM5.0 Biomass and Carbon Exchange across the Western US Using a Data Assimilation System.
       *Journal of Advances in Modeling Earth Systems*, `doi.org/10.1029/2020MS002421 <https://doi.org/10.1029/2020MS002421>`__


.. |ELM gridcell breakdown| image:: ../../guide/images/clm_landcover.png
   :height: 600px
   :target: https://escomp.github.io/ctsm-docs/versions/release-clm5.0/html/tech_note/Ecosystem/CLM50_Tech_Note_Ecosystem.html#surface-characterization
