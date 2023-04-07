.. title:: files_and_modules

.. _files_and_modules:

============
File and Modules
============

vanDANA is a python package with main executable (vanDANA.py) and three submodules (common, utilities, user_inputs).

vanDANA.py
==========

This is the main executable file that pulls information from other modules (namely: user inputs, utilities, common) and runs the main function vanDANA_solver (args) which sets the entire workflow. The user need not require making any changes to vanDANA.py. Before the time loop, the following preliminary operations are conducted in order:

#. Set the current directory and MPI controls.
#. Read mesh files
#. Load physical problem classes.
#. Setup initial and boundary conditions.
#. Preassemble the flow and temperature problems
#. Create output files for postprocessing
#. Read restart files (only if restarting the simulation)
