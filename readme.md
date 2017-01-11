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


##Files

`ed_config_readable.dat`

parameter input file; readable

(file is altered by routines and then saved in ed_config.dat which is not readable but used by other routines)

**./Source**	source code


- `complile.sh` Compiles the ed_main.cpp with all its dependencies TODO: Solve with makefile

- /cpp_routines
 - `ed_main.cpp` Builds Hamiltonian and diagonalizes it.

 - `/ed_hamiltonian.cpp` Saves routines that correspond to terms in Hamiltonian (e.g. interaction, NN hopping etc.)

 - `ed_build_dict_gen.cpp`
Builds a state 'dictionary' by applying the symmetries saved in the config file.
Run after modifying symmetry.			

 - `ed_shared.cpp`
Routines that are shared by ed_main.cpp and ed_build_dict_gen.cpp

 - `ed_update_config.cpp`
Updates the configuration file. This is necessary if several representations 
of a symmetry group want to be computed in a row (e.g. different momenta).


**./Utils**		different utilities

- `sym.py`
Create character string that describes symmetry operations 
in "config_readable.dat"	

```
For example mirror symmetry for 6 sites is described by
1 2 3 4 5 6
6 5 4 3 2 1
```

- **/R**		R scripts that I use to plot the spectra	

