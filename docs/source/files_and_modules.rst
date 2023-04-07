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

To obtain more control at the command prompt, we use the `argparse library <https://docs.python.org/3/library/argparse.html>`__ which allows us to add keywords to the executable, for eg: ``mpirun.mpich -n 64 python3 vanDANA.py -restart=True -T=20 -velocity_degree=1``. Note that we also provide a separate batch script (`src/fenics.sb <https://github.com/patelte8/vanDANA/blob/IB-FSI/src/fenics.sb>`__) to submit jobs on the HPC cluster.

Our Immersed Boundary FSI formulation is based on the Distributed Lagrange Multiplier (DLM) based Fictitious Domain (FD) method and the physical problem is decomposed into the fluid, solid and Lagrange multiplier sub-problems. The entire code is non-dimensional, and we use the Finite Element Method (FEM) for spatial discretization of the governing equations. For a detailed understanding of the numerical scheme, we advise the reader to refer to Yu :sup:`[1]` and Yu et al :sup:`[2]`.

**Flow problem: IPCS scheme**

**Step 1 - Velocity prediction**: solve for  solve for :math:`u^{*}`

.. math::

  \begin{aligned}
  \left\langle\frac{\boldsymbol{u}^{*}-\boldsymbol{u}^{\boldsymbol{\theta}}}{\Delta \tau}, \boldsymbol{v}\right\rangle+ & \left\langle\frac{3 \boldsymbol{u}^{\boldsymbol{\theta}}-\boldsymbol{u}^{\boldsymbol{\theta}-1}}{2} \cdot \nabla \boldsymbol{u}_{\boldsymbol{c k}}^{*}, \boldsymbol{v}\right\rangle \\
  & +\left\langle\frac{\nabla \boldsymbol{u}_{\boldsymbol{c k}}^{*}+\nabla \boldsymbol{u}_{\boldsymbol{c k}}^{* T}}{R e}-p^{\theta} \boldsymbol{I}, \nabla \boldsymbol{v}\right\rangle+\left\langle p^{\theta} \cdot \boldsymbol{n}, \boldsymbol{v}\right\rangle_{\partial \boldsymbol{\Omega}}-\left\langle\frac{\widehat{\boldsymbol{g}}}{F r^{2}}, \boldsymbol{v}\right\rangle \\
  & +\left\langle\gamma_{S U P G}\left(\boldsymbol{u}^{\boldsymbol{\theta}}\right) P\left(\boldsymbol{u}^{\boldsymbol{\theta}}, \boldsymbol{v}\right), \boldsymbol{R}^{\boldsymbol{\theta}}\right\rangle+\left\langle\gamma_{C W}\left(\boldsymbol{u}^{\boldsymbol{\theta}}\right) \Lambda\left(\boldsymbol{u}^{\boldsymbol{\theta}}, \boldsymbol{u}^{*}\right), \nabla \boldsymbol{v}\right\rangle=\left\langle\lambda^{\boldsymbol{\theta}}, \boldsymbol{v}\right\rangle_{\boldsymbol{P}}


  \end{aligned}

(1)