#
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 17O_CANOPS
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
Stochastic inversion of 17O for atmospheric pO2 using CANOPS model output
#
These scripts will perform a filtered inversion for estimating atmospheric pO2 based on the oxygen isotope composition of lacustrine or marine sulfate minerals.  The model parameters and equations are described in N. Planavsky, C. Reinhard, T. Isson, K. Ozaki, and P. Crockford [2020] "Large mass-independent oxygen isotope fractionations in mid-Proterozoic sediments: Strong evidence for a low-oxygen atmosphere?, Astrobiology, doi:10.1089/ast.2019.2060.  The scripts can be used to produce Figures 1, 4, and 5 in Planavsky et al. [2020], along with an additional figure comparing filtered and unfiltered results and some summary statistics.
#
The main model script [canops_17O.m] first evaluates an inclusive range of atmospheric pO2 and gross primary productivity [GPP] values that is consistent with a given sulfate oxygen isotope value [using paleoGPP.m].  The model then uses results from the CANOPS global biogeochemical cycle model to filter the mass balance results baed on combinations of atmospheric pO2 and GPP that result in a redox-balanced Earth system [see Ozaki et al, 2019, doi:10.1111/gbi.12317].
#
As a courtesy, we request that others who use this code [or implement the underlying approach] please cite Planavsky et al. [2020].  We also request that those who use/modify the code please send publications and/or modified code to NJP [noah.planavsky@yale.edu] or CTR [chris.reinhard@eas.gatech.edu].
#
# REQUIREMENTS: Written/tested on Matlab R2018a/R2018b/R2019a [but should be backwards-compatible]
#
TO RUN THE CODE
#
User can call canops_17O from the working directory, which will load sulfate 17O data and the required CANOPS outputs and run an analysis using default values for all parameters [these are specified near the top of canops_17O.m].
#
Alternatively, user-specified parameters can be supplied according to: canops_17O(params)
#
User must choose:
- num_run              --> number of samples for stochastic routine 
- f_O_sulfate_min/max  --> minimum/maximum O atom incorporation into sulfate
- logPALO2_min/max     --> minimum/maximum atmospheric pO2 range for mass balance calculations [PAL]
- CO2_min/max          --> minimum/maximum atmospheric pCO2 range for mass balance calculations [PAL]
- GPP_mod              --> modern global net primary productivity [GtC y-1]
- f_ex                 --> export ratio reltaive to marine GPP [dimensionless]
- p_fit                --> prediction interval for CANOPS filter
#
Finally, user can instead call make_clean, which will remove all existing output data and figures, then call make_inversion, which will run the analysis from scratch
#
------------------------------------------------------------------------
last run on: 2020-03-18 by CTR
------------------------------------------------------------------------
median_pO2 = 3.5E-03
q90_pO2    = [1.4E-03 8.9E-03]
------------------------------------------------------------------------
#
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++