#ed_ising

Exact diagonalization procedures for the Ising model (can be easily extended to
other 1d - Hamiltonians).

Sebastian Dick; Dec 2016

##How-to

- add desired hamiltonian terms to `ed_hamiltonian.cpp` and adjust the global variables accordingly
- adjust the input file `ed_config_readable.dat` to include your desired parameters and symmetries (so far only non-negative integers possible)
- switch to */Source* and run `make`, `make build_dict` and `make ed_update`, this will create the required executables
- the loop in the script `run.sh` controls the number of computed momenta. If there is no translational symmetry, make sure to set the upper limit to 1
- execute `run.sh` 
- eigenvalues are stored in `out.dat`

##Dependencies

This program uses the linear algebra package **Armadillo**, which is available as open source, using the Mozilla Public License (MPL) 2.0. and can be downloaded at http://arma.sourceforge.net. For more information please see:

Conrad Sanderson and Ryan Curtin. 
*Armadillo: a template-based C++ library for linear algebra.*
Journal of Open Source Software, Vol. 1, pp. 26, 2016.



