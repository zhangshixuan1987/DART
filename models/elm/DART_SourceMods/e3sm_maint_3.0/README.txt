DART SourceMods for ELM-DART on E3SM maint-3.0
============================

These source modifications override ELM source files during an E3SM build.
Place this directory structure under your E3SM case SourceMods directory:

  <caseroot>/SourceMods/src.elm/

These SourceMods were verified directly against the local E3SM checkout at
/qfs/people/zhan391/e3sm_dart_work/code/E3SMv3. Relative to that source,
they contain only the DART changes documented below.

Verified E3SM branch: v3_dart_cda
Verified E3SM revision: 1eda4d4dcacd9ff603764264d5bf8ca861b97ac1
Nearest version: v3.0.2-96-g1eda4d4dca
Reference ELM path: components/elm/src/


Files Modified
--------------

src.elm/EcosystemBalanceCheckMod.F90
  PURPOSE: Skip C/N/P balance checks on the first restart step after a DART
  assimilation. DART may create or destroy mass during the update step, so
  the balance error check on startup is expected and should not abort the run.
  MODIFICATION: Added is_first_restart_step() guard to all err_found checks
  in ColCBalanceCheck, ColNBalanceCheck, ColPBalanceCheck, and GridCBalanceCheck.

src.elm/SurfaceRadiationMod.F90
  PURPOSE: Expose the PARVEG (absorbed PAR by vegetation) history variable.
  DART's QTY_ABSORBED_PAR forward operator (obs_def_land_mod.f90) requires
  the 'PARVEG' history variable. ELM computes this value internally but only
  stores 'PARVEGLN' (local noon value). This modification adds 'PARVEG' as
  an every-timestep history output variable.
  MODIFICATION: Added parveg_patch member, allocation, history registration,
  and assignment from the existing parveg local variable.


Files Required but Not Yet Ported (TODO)
-----------------------------------------

src.elm/biogeophys/CanopyFluxesMod.F90
src.elm/biogeophys/PhotosynthesisMod.F90
  PURPOSE: Required to support the Solar-Induced Fluorescence (SIF) DART
  forward operator (QTY_SOLAR_INDUCED_FLUORESCENCE). In the legacy CLM SourceMods (recoverable from Git commit d4009d7c4),
  these files were modified to expose the
  per-patch sun-lit and shaded net assimilation rates (anetsun_patch,
  anetsha_patch) used by the SIF forward operator.
  ELM does not currently have these variables. Porting these modifications
  requires adding the SIF-related state variables to ELM's photosynthesis
  data structures.
  Skip these if you do not plan to assimilate SIF observations.


Installation
------------

To use these SourceMods in an E3SM assimilation case, copy the files into
the case SourceMods directory before building:

  cp -p <dartroot>/models/elm/DART_SourceMods/e3sm_maint_3.0/src.elm/*.F90 \
        <caseroot>/SourceMods/src.elm/
