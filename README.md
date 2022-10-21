# schism_Nea

SCHISM, Original Code cloned from https://github.com/schism-dev/schism
with Modifications started on 2022 October 16 by J.Lefevre, IRD :

Updates : 21-Oct-2022
---------------------
- EXPERIMENTAL Branch: Add WindDragPowell as one option to compute Cdrag by storm sectors in PaHM. Induce some changes in Hydro/schism_step.F90
-                      Add USE_POWELL in include_modules
-                      See rationale from Powell et al 2006 : “Final Report to the National Oceanic and Atmospheric Administration (NOAA) Joint 
                           Hurricane Testbed (JHT) Program.” 26 pp
-                      Add one png file in Core/Pahm/inputs with Niran using PAHM+POWELL (using "Prate" as output Bucket for Cdrag) 
- This is not a recommended commit, only to keep a trace of WindDragPowell() somewhere. 
- Keep care on this : 
"The application of Powell wind drag in modelling storm surge in “small islands” can lead to weird variation of water levels, specially unrealistic negative wind setup on the leeward side. This is because Powell always assume opposing waves in the rear sector, and consequently large drag values, which is not the case on the leeward of small islands. source A.E.Hermosa 

Updates : 17-Oct-2022
---------------------
- Add GAHM (Generalized Asymmetric Holland Model) in PaHM/parwind.F90 as a second alternative to Holland Model 
 (see details https://wiki.adcirc.org/Generalized_Asymmetric_Holland_Model and https://noaa-ocs-modeling.github.io/PaHM/pahm_manual.pdf),
- Add support for Hurricanes in both South and North Hemisphere (in HM and GAHM models),
- Test using Niran, SW Pacific: (See inside src/Core/PaHM/inputs/ : a track for Cyclone Niran "niran2021-bdeck.dat" and my comparison with HM and GAHM versus a weather model output)
- CAUTION // CAUTION : 
  To switch from HM (1) or GAHM (10), the user still need to change the value of "modelType" in Pahm_Utilities.F90 (line 3210) and recompile schism 
- CAUTION / CAUTION : Unlike in noaa-ocs-modeling.github.io/PaHM, there isnot Control File support in SCHISM/PaHM yet
---------------------
Updates End

# SCHISM

The **S**emi-implicit **C**ross-scale **H**ydroscience **I**ntegrated **S**ystem **M**odel (SCHISM) is an open-source community-supported modeling system based on unstructured grids and designed for the seamless simulation of 3D baroclinic circulation across creek-lake-river-estuary-shelf-ocean scales.

# Building and documentation

The manual may be found on the SCHISM wiki at http://ccrm.vims.edu/schismweb/. Build instructions are described in Chapter 1.

The online documentation can be accessed at https://schism-dev.github.io/schism.

# Developing and contributing

When using the development version, note changes in flags and features described in `src/Readme.beta_notes` and `sample_inputs/param.nml`, `sample_inputs/bctides.in`, etc.

Please refer to `CONTRIBUTING.md` for more information on contributing to SCHISM.
