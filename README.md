# SiZer_manifold

Repo for codes used for "Inference on the shape of densities on Riemannian manifolds via SiZer" (Journal of Korean Statistical Society, DOI 10.1007/s42952-024-00301-3)  

Currently four R codes are available:

 - sphere_code.R: draw a SiZer map if spherical data is inputted. Before run this file, make sure to run all of sphere_function.R.

 - sphere_function.R: contains functions for performing SiZer inference on spherical data, along with auxiliary functions. After run the code in this file, visualization can be performed using sphere_code.R.

 - torus_code.R: draw a SiZer map if toroidal data is inputted. Before run this file, make sure to run all of torus_function.R.

 - torus_function.R: contains functions for performing SiZer inference on toroidal data, along with auxiliary functions. After run the code in this file, visualization can be performed using torus_code.R.

query.csv contains the seismic activity data retrieved from U.S. Geological Survey (USGS, https://www.usgs.gov/).
