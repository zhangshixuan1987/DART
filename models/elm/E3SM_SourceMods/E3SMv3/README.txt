E3SM_SourceMods for ELM-DART
============================

These source modifications override ELM source files during an E3SM build.
Place this directory structure under your E3SM case SourceMods directory:

  <caseroot>/SourceMods/src.elm/

The modifications are based on E3SMv3. Verify compatibility if using a
different E3SM tag.

Reference E3SM tag: E3SMv3
Reference ELM path: components/elm/src/


Files Modified
--------------

src.elm/biogeochem/EcosystemBalanceCheckMod.F90
  PURPOSE: Skip C/N/P balance checks on the first restart step after a DART
  assimilation. DART may create or destroy mass during the update step, so
  the balance error check on startup is expected and should not abort the run.
  MODIFICATION: Added is_first_restart_step() guard to all err_found checks
  in ColCBalanceCheck, ColNBalanceCheck, ColPBalanceCheck, and GridCBalanceCheck.

src.elm/biogeophys/SurfaceRadiationMod.F90
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
  forward operator (QTY_SOLAR_INDUCED_FLUORESCENCE). In the CLM SourceMods
  (DART_SourceMods/cesm2_2_0), these files were modified to expose the
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

  cp -r <dartroot>/models/elm/E3SM_SourceMods/E3SMv3/src.elm/* \
        <caseroot>/SourceMods/src.elm/
