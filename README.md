# Sulfur_X_kdonly

This is a sub-package of Sulfur_X that can be used to calculate the sulfur partition coefficients with given P, T composition, H2O, and H2O mole fraction in the vapor.

Please cite as follows: Ding, S., Plank, T., Wallace, P., Rasmussen, D. J., in review. Sulfur_X: A model of sulfur degassing during magma ascent. Geochemistry, Geophysics, Geosystems. https://doi.org/10.31223/X56H0F. The related manuscritp is available on EarthArxiv: https://eartharxiv.org/repository/view/3559/

Input file
sample_name: samples' names
pressure: P in MPa
temperature: T in C
fO2: fO2 relative to FMQ
SiO2... H2O: major elements and H2O in wt.%
XH2Of: mole fraction of h2o in the vapor

Output: partition_coefficient.csv
All Kds are calculated as kd = XS_vapor/XS_melt, XS is mole fraction
