# Single-Elastic-Constant-SEC-Determination
This repostory contains an algorithm of a method to determine single elastic constants 
(SEC) also known as (directional moduli) of each phase within a dual-phase material and 
also the global material SEC of a mono-phase or dual-phase material using data from 
high energy Xray Diffraction (HEXRD). 


It was developed by Muritala Oluwaseun Arowolo during his PhD work at the 
University of Toulouse, Toulouse France within two CNRS research laboratories:

a) CIRIMAT - situated in the site of Institute National Polytechnic (Toulouse INP) &
b) ICA - situated in the site of Ecole des Mines Albi (IMT Albi)


###............................NOTE.....................................
This algorithm originally handles two phase materials with HCP and BCC crystals. 
It can be easily adapted to other crystals like FCC, tetragonal, etc (to be updated later). 
All the codes are written in python programming language.

............................PROCEDURE...................................

(1) Obtain the directional d-spacings along the loading direction (LD) and 
transverse direction (TD) of XRD patterns. 

A jupyer notebook script is included in the file named "2D-Integration" which maybe 
used to perform 2D-Integration/partial integration/peak-to-peak analyses of the XRD 
patterns. This script can return the following from the Pseudovoigt fittings of the 
diffractograms: 
[a] 2θ - positions of the spectral peaks, 
[b] peaks intensities, 
[c] FWHM, 
[d] d-spacings and 
[d] the images of the fitted patterns for visual observation.
 
NB: only d-spacings are needed from the patterns for the SEC determination.

(2) Calculate the hkl lattice micro-strain from the the d-spacings obtained from (1)
along 0°, 90°, 180° and 270°. Take the average over 0° and 180°, likewise 90° and 270°.

A script named "Compute_mean_microstrain_direction.py" may be used to compute the 
average. A file named "Average_microstrain_data.xlsx" is returned which contains 
the results.

(3) These average lattice micro-strain should be provided in the sheet named 
'Microstrain_FIT_ELASTIC' inside an excel file named "SEC_d-spacing_stress_strain_data.xlsx" 
in the columns with the appropriate hkl headers as shown in the file. 
Also, the macroscopic stress vs strain data should also be provided in the sheet named 
'Stress-exp-elastic'. The macroscopic strain should be included in the a column header 
'Macro-True_Strain' in the 'Microstrain_FIT_ELASTIC' sheet.

(4) The SEC determination algorithm is in a file named "Genetic_alg_python_elastic_constants.py".
The "SEC_d-spacing_stress_strain_data.xlsx" should be provided with the same folder as the 
SEC determination algorithm. Further, the following should be provied in the algorithm:
[a] The SEC bounds (upper and lower) of each phase in the material 
[b] The volume fraction of each phase.
[b] For HCP crytal, the c/a ratio.
