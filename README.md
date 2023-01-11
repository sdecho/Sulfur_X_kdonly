# Sulfur_X_kdonly
This package can be used to calculate the sulfur partition coefficients with given P, T composition, H2O, and H2O mole fraction in the vapor.
Input file
sample_name: samples' names
pressure: P in MPa
temperature: T in C
fO2: fO2 relative to FMQ
SiO2... H2O: major elements and H2O in wt.%
XH2Of: mole fraction of h2o in the vapor

Output: partition_coefficient.csv
All Kds are calculated as kd = XS_vapor/XS_melt, XS is mole fraction
