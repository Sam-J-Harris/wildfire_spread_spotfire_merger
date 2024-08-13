# Wildfire spread and spotfire merger
Supplementary material to [1]. 
This is preliminary code and test data that has been provided, please contact the author (sam.harris.16@ucl.ac.uk) if you have any issues.

---- DATA ----

The .mat and .fig files for all figures in [1] can be found in the ``00 Data'' file.
The code used to generate each figure is given in the corresponding .m files.

---- CODE ----

Two main codes are created: ROSS and ROSA.
ROSS is the code to model single wildfire spread. It uses conformal mapping based methods.
ROSA is the code to model multiple spotfires spread and merger. It uses the AAA-least squares algorithm.

General code is given in the A1_run_ROSS.m and A2_run_ROSA.m files. 

References 
[1] Harris S.J., McDonald N.R. "Modelling wildfire spread and spotfire merger using conformal mapping and AAA-least squares methods."
