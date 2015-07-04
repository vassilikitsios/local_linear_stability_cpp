----------------------------------------------------------------------
1) Code overview:
----------------------------------------------------------------------

Author: Vassili Kitsios

This C++ code performs local stability calculations for a series of analytical profiles and profiles read in from file.

This README file contains installation instructions and an overview of the code. If you are to use this code in your own projects please cite the following documents:

Kitsios, V., Cordier, L., Bonnet, J.-P., Ooi, A. & Soria, J., 2010, Development of a nonlinear eddy viscosity closure for the stability analysis of a turbulent channel flow, Journal of Fluid Mechanics, Vol. 664, pp 74-107.
http://journals.cambridge.org/action/displayAbstract?fromPage=online&aid=7928161

Kitsios, V., 2010, Recovery of fluid mechanical modes in unsteady separated flows, PhD Thesis, The University of Melbourne
https://minerva-access.unimelb.edu.au/handle/11343/35705


----------------------------------------------------------------------
2) Installation instructions:
----------------------------------------------------------------------

Prior to compliling this code the following libraries are required on the system:

	bliz++		http://www.oonumerics.org/blitz/

	clapack		http://www.netlib.org/clapack/

	gnuplot 	http://www.gnuplot.info/		(required only for plotting the output)


Once these components are installed see "Makefile" for details on how to compile the present code.


----------------------------------------------------------------------
3) List of files and directories with brief explanations: 
----------------------------------------------------------------------

drivers/					- directory containing executable directory and Makefile

	/Makefile				- Makefile to jump make executable directory

	/local_linear_stability/ 		- directory containing executable code


include/					- directory containing include files from libraries made from src directory


lib/						- directory containing libraries made from src directory


Makefile					- main makefile

Makefile.in					- input makefile with system specific settings. See "Makefile" to see how to create your own "Makefile.in"

makefiles/					- directory containing system specific makefiles


src/						- directory containing the library source directories

	/libEigenSolver/			- directory containing the source for the non-Hermitian complex eigensolver library

	/libLocalStability/			- directory containing the source for the local stability eigenvalue problem (EVP) library

	/libUtils/				- directory containing the source for a general Utilities library


tests/						- directory containing Blasius boundary layer and turbulent channel test cases


tests/blasius/					- stability analysis of a laminar Blasius boundary layer, see PhD thesis appendix B

	/1mapping/				- results directory containing mapping of EVP from the complex wavenumber to frequency plane

	/2most_unstable_mode/			- results directory containing the most spatial unstable wave calculated at a higher resolution


tests/turbulent_channel/			- triple decomposition stability analysis of a turbulent channel, see PhD thesis chapter 6

	/1raw_profiles/ 			- directory containing raw base flow

	/2interpolated_profiles/		- raw base flow interpolated onto Chebyshev profiles of various resolutions

	/3stability_laminar/ 			- results directory containing stability analysis using the laminar closure

	/4stability_linear_eddy/		- results directory containing stability analysis using the linear eddy viscosity closure 

	/5stability_nonlinear_eddy/		- results directory containing stability analysis using the non-linear eddy viscosity closure 


----------------------------------------------------------------------
3) List of files and directories within the results directories: 
----------------------------------------------------------------------

local_linear_stability				- executable

local_linear_stability.in			- input deck containing parameters specifying the EVP

local_linear_stability.run			- batch file to execute the code


results/					- directory containing the eigvalue and eigenvector files


images/						- directory containing the images created from the results


----------------------------------------------------------------------
