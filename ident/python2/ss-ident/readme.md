# Identifiability Analysis and Experimental Design:
###### Generating *in silico* data for identifiability analysis:
- run `create_experiment_data` with desired options active to generate initial set of experimental data
`mpirun -np 4 python create_experiment_data.py`
- all generated data files are stored in *exp* subdirectory

- run `ident_exp_data` with desired options to generate data needed for identifiability analysis
`python ident_exp_data.py`
- all generated data files are stored in *exp* subdirectory

###### Perform identifiability analysis:
- run `parallel_ident` with desired functions to perform identifiability analysis on the small network
- run the script with `mpi` to enable multithreaded operation
- e.g., `mpirun -np 4 python parallel_ident.py`
- all generated data files are stored in *ident* subdirectory
- all generated figures are stored in the *results* subdirectory

###### Validate estimated parameters:
- all estimated parameters from the identifiability analysis in the previous step can be validated
- run `parallel_validate` with desired function to perform validation of identified and estimated parameters
- run the script with `mpi` to enable multithreaded operation
- e.g., `mpirun -np 4 python parallel_validate.py`

# Prerequisites:
- Python 3.4
- mpi4py `conda install mpi4py`
- [mpi4py master-slave library](https://github.com/luca-s/mpi-master-slave/blob/master/examples/example4.py) 