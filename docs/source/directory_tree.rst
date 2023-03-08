.. title:: Directory Tree

.. _directory_tree:

============
Directory Tree
============

"""
    vanDANA-IB-FSI
     ┃
     ┣ __init__.py
     ┣ vanDANA.py
     ┣ LICENSE.md
     ┃
     ┣ common
     ┃ ┣ __init__.py
     ┃ ┣ constitutive_eq.py
     ┃ ┣ delta_interpolation.py
     ┃ ┣ fem_stabilizations.py
     ┃ ┣ flow_temperature_variational_problem.py
     ┃ ┣ flow_variational_problem.py
     ┃ ┣ functions.py
     ┃ ┣ lagrange_variational_problem.py
     ┃ ┣ solid_variational_problem.py
     ┃ ┗ solver_options.py
     ┃
     ┣ user_inputs
     ┃ ┣ __init__.py
     ┃ ┣ boundary_initial_conditions.py
     ┃ ┣ problem_specific.py
     ┃ ┣ user_parameters.py
     ┃ ┗ ...
     ┃
     ┣ utilities
     ┃ ┣ __init__.py
     ┃.┣ utils.py
     ┃ ┣ write.py
     ┃ ┗ read.py   
     ┃ 
     ┣ src
     ┃ ┣ flag.geo
     ┃ ┣ mesh.py
     ┃ ┣ fenics.sb
     ┃ ┣ fenics_2019_dev 
     ┃ ┗ ...
     ┃
     ┗ results
       ┣ HDF5_files
       ┣ mesh_files
       ┣ restart_variables
       ┣ XDMF_files 
       ┗ text_files
         ┣ flow_data.txt
         ┣ flow_temp_data.txt
         ┣ lagrange_data.txt
         ┣ log_info.txt
         ┣ restart.txt
         ┣ runtime_stats.txt
         ┣ solid_data.txt
         ┣ solid_mesh_quality.txt
         ┗ solid_temp_data.txt
"""
