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

**Step 1 - Velocity prediction**: solve for :math:`u^{*}`

.. math::
  
  \begin{aligned}
  \left\langle\frac{\boldsymbol{u}^{*}-\boldsymbol{u}^{\boldsymbol{\theta}}}{\Delta \tau}, \boldsymbol{v}\right\rangle+ & \left\langle\frac{3 \boldsymbol{u}^{\boldsymbol{\theta}}-\boldsymbol{u}^{\boldsymbol{\theta}-1}}{2} \cdot \nabla \boldsymbol{u}_{\boldsymbol{c k}}^{*}, \boldsymbol{v}\right\rangle \\
  & +\left\langle\frac{\nabla \boldsymbol{u}_{\boldsymbol{c k}}^{*}+\nabla \boldsymbol{u}_{\boldsymbol{c k}}^{* T}}{R e}-p^{\theta} \boldsymbol{I}, \nabla \boldsymbol{v}\right\rangle+\left\langle p^{\theta} \cdot \boldsymbol{n}, \boldsymbol{v}\right\rangle_{\partial \boldsymbol{\Omega}}-\left\langle\frac{\widehat{\boldsymbol{g}}}{F r^{2}}, \boldsymbol{v}\right\rangle \\
  & +\left\langle\gamma_{S U P G}\left(\boldsymbol{u}^{\boldsymbol{\theta}}\right) P\left(\boldsymbol{u}^{\boldsymbol{\theta}}, \boldsymbol{v}\right), \boldsymbol{R}^{\boldsymbol{\theta}}\right\rangle+\left\langle\gamma_{C W}\left(\boldsymbol{u}^{\boldsymbol{\theta}}\right) \Lambda\left(\boldsymbol{u}^{\boldsymbol{\theta}}, \boldsymbol{u}^{*}\right), \nabla \boldsymbol{v}\right\rangle=\left\langle\lambda^{\boldsymbol{\theta}}, \boldsymbol{v}\right\rangle_{\boldsymbol{P}} 
  \end{aligned} 

  
where, :math:`\boldsymbol{u}^{*}` is the unknown intermediate flow velocity and :math:`\boldsymbol{u}_{\boldsymbol{c k}}^{*}=\frac{\boldsymbol{u}^{*}+\boldsymbol{u}^{\boldsymbol{\theta}}}{2}` is the Crank-Nicolson velocity. In the above equation, the vector Lagrange multiplier :math:`\boldsymbol{\lambda}` serves as a pseudo body force only in the overlapping solid region :math:`\boldsymbol{P}` and needs to be interpolated from the Lagrangian mesh via the smeared delta functions.

**Step 2 - Pressure projection :math:`\&` velocity correction**:
solve for :math:`p^{\theta+1}` and correct :math:`u^{\Phi}`.

.. math::

  \begin{gathered}
  \left\langle\nabla\left(p^{\theta+1}-p^{\theta}\right), \nabla q\right\rangle+\frac{1}{\Delta \tau}\left\langle\nabla \cdot \boldsymbol{u}^{*}, q\right\rangle+\left\langle\gamma_{P S P G}\left(\boldsymbol{u}^{\theta}\right) \nabla q, \boldsymbol{R}^{\boldsymbol{\theta}}\right\rangle=0 \end{gathered} \label{b}   \tag{2}


.. math::

  \begin{gathered}
  \left\langle\frac{\boldsymbol{u}^{\boldsymbol{\phi}}-\boldsymbol{u}^{*}}{\Delta \tau}, \boldsymbol{v}\right\rangle+\left\langle\nabla\left(p^{\theta+1}-p^{\theta}\right), \boldsymbol{v}\right\rangle=0
  \end{gathered} \label{c}   \tag{3}


**Solid problem (compressible)**: solve for :math:`\Delta x^{\theta+1}` on the reference configuration.

.. math::

  \begin{aligned}
  &\left\langle\rho_{r} \frac{\Delta \boldsymbol{x}^{\boldsymbol{\theta}+\mathbf{1}}}{\Delta \tau^{2}}, \boldsymbol{\omega}\right\rangle_{\boldsymbol{P}_{\mathbf{0}}} \\
  &=\left\langle\frac{\boldsymbol{u}^{\boldsymbol{\phi}}}{\Delta \tau}, \boldsymbol{\omega}\right\rangle_{\boldsymbol{P}_{\mathbf{0}}}+\left\langle\rho_{r}-1 \frac{\Delta \boldsymbol{x}^{\boldsymbol{\theta}}}{\Delta \tau^{2}}, \boldsymbol{\omega}\right\rangle_{\boldsymbol{P}_{\mathbf{0}}}-\left\langle\nabla \boldsymbol{\omega}^{T}, \lambda_{0} \ln J \boldsymbol{F}^{-\mathbf{1}}+G\left(\boldsymbol{F}^{\boldsymbol{T}}-\boldsymbol{F}^{-\mathbf{1}}\right)\right\rangle_{\mathbf{P}_{\mathbf{0}}}^{\boldsymbol{\theta}+\mathbf{1}} \\
  &+\left\langle\rho_{r}-1 \frac{\widehat{\boldsymbol{g}}}{F r^{2}}, \boldsymbol{\omega}\right\rangle_{\boldsymbol{P}_{\mathbf{0}}}-\left\langle\lambda^{\boldsymbol{\theta}}, \boldsymbol{\omega}\right\rangle_{\boldsymbol{P}}
  \end{aligned}

**Lagrange multiplier problem**: solve for :math:`\lambda^{\theta+1}` on the current configuration.

.. math::


  \left\langle\frac{\Delta \boldsymbol{x}^{\boldsymbol{\theta}+1}}{\Delta \tau^{2}}-\frac{\boldsymbol{u}^{\boldsymbol{\phi}}}{\Delta \tau}, \zeta\right\rangle_{P}=\left\langle\lambda^{\boldsymbol{\theta}+1}-\lambda^{\boldsymbol{\theta}}, \boldsymbol{\zeta}\right\rangle_{P} \label{e}   \tag{5}

In the above equations :math:`\langle`,\ :math:`\rangle` represents the Euclidean inner product in the respective domain, :math:`\boldsymbol{u}` and :math:`p` are the fluid velocity and pressure, :math:`\boldsymbol{x}` is the solid displacement, :math:`\boldsymbol{F}` is the deformation gradient tensor and :math:`\rho_{r}` is the ratio of solid to fluid densities. Here, :math:`\Omega` is the overall fluid domain including the fictitious fluid and :math:`P_{0}` and :math:`P` are the reference and current solid configurations. The above DLM/FD algorithm is implemented in the time loop in **vanDANA_solver (args)**, see Figure 2 below.

.. code:: python

      # Time loop
      try:

          while T > tsp and t < T:
              
              timer_dt.start()
              update_counter(counters)

              # Update current time
              t += tsp   

              # Update boundary conditions : only if time-dependent
              # parabolic_profile.t = t; tim.t = t; num_cycle.cycle = int(t / t_period)     
              # for ui, value in inflow.items():     
                  #  inflow[ui][0].v = evaluate_boundary_val(param_LSPV); inflow[ui][1].v = evaluate_boundary_val(param_LIPV)
                  #  inflow[ui][2].v = evaluate_boundary_val(param_RSPV); inflow[ui][3].v = evaluate_boundary_val(param_RIPV)

              if problem_physics['solve_FSI'] == True:
                  timer_si.start()
                  Lm_f.assign(interpolate_nonmatching_mesh_delta(fsi_interpolation, Lm_[1], FS['fluid'][2], interpolation_fx, "F"))
                  si += timer_si.stop()
                  
              timer_s1.start()
              # print(BLUE % "1: Predict tentative velocity step", flush = True)
              A1, b1 = flow.assemble_tentative_velocity(u_, p_, Lm_f, dt)
              flow.solve_tentative_velocity(A1, u_[0], b1, bcs['velocity'])
              s1 += timer_s1.stop()

              timer_s2.start()
              # print(BLUE % "2: Pressure correction step", flush = True)
              b2 = flow.assemble_pressure_correction(u_, p_, Lm_f, dt)
              flow.solve_pressure_correction(p_[0], b2, bcs['pressure'])
              s2 += timer_s2.stop()

              timer_s3.start()
              # print(BLUE % "3: Velocity correction step", flush = True)
              b3 = flow.assemble_velocity_correction(u_, p_, dt)
              flow.solve_velocity_correction(u_[0], b3, bcs['velocity'])
              s3 += timer_s3.stop()

              assigner_uv.assign(uv, [u_[0][ui] for ui in range(u_components)])

              # --------------------------------------------------------------------------------- 

              if problem_physics['solve_FSI'] and problem_physics['solve_temperature'] == True:
                  timer_si.start()
                  LmTf_.assign(interpolate_nonmatching_mesh_delta(fsi_interpolation, LmTs_[1], FS['fluid_temp'][0], interpolation_fx, "F"))
                  si += timer_si.stop()

              timer_s4.start()
              # print(BLUE % "4: Energy conservation step", flush = True)
              if problem_physics['solve_temperature'] == True:
                  A4, b4 = flow_temp.assemble_temperature(T_, uv, LmTf_, dt)
                  flow_temp.solve_temperature(A4, T_[0], b4, bcs['temperature'])
              s4 += timer_s4.stop()	    

              # --------------------------------------------------------------------------------- 

              if problem_physics['solve_FSI'] == True:
                  timer_si.start()
                  uf_.assign(interpolate_nonmatching_mesh_delta(fsi_interpolation, uv, FS['lagrange'][0], interpolation_fx, "S"))
                  si += timer_si.stop()

              timer_s5.start()    
              # print(BLUE % "5: Solid momentum eq. step", flush = True)    
              if problem_physics['solve_FSI'] == True:    
                  a5 = solid.assemble_solid_problem(problem_physics['compressible_solid'], Dp_, mix, uf_, Lm_[1], dt)
                  try:
                      solid.solve_solid_displacement(solid_mesh_R.mesh, problem_physics['compressible_solid'], a5, Dp_[1], mix, ps_, p_[0], bcs['solid'])
                  except:
                      solid.change_initial_guess(Dp_[1], mix)	        		        	
                      solid.solve_solid_displacement(solid_mesh_R.mesh, problem_physics['compressible_solid'], a5, Dp_[1], mix, ps_, p_[0], bcs['solid'])

                  Dp_[0].vector().axpy(1.0, Dp_[1].vector())
                  # solid.compute_jacobian(J_, Dp_[0])

                  us_.vector().zero()
                  us_.vector().axpy(1/float(dt), Dp_[1].vector())
              s5 += timer_s5.stop()
              
              # --------------------------------------------------------------------------------- 

              timer_s6.start()
              # print(BLUE % "6: Lagrange multiplier (fictitious force) step", flush = True)
              if problem_physics['solve_FSI'] == True:
                  a6, b6 = lagrange.assemble_lagrange_multiplier(Lm_, us_, uf_, dt)
                  lagrange.solve_lagrange_multiplier(a6, Lm_[0], b6)
              s6 += timer_s6.stop()    

              # --------------------------------------------------------------------------------- 

              if problem_physics['solve_FSI'] and problem_physics['solve_temperature'] == True:
                  timer_si.start()
                  Ts_[0].assign(interpolate_nonmatching_mesh_delta(fsi_interpolation, T_[0], FS['solid_temp'][1], interpolation_fx, "S"))
                  si += timer_si.stop()

              timer_s7.start()
              # print(BLUE % "7: Solid temperature based lagrange multiplier step", flush = True)
              if problem_physics['solve_FSI'] and problem_physics['solve_temperature'] == True:
                  a7, b7 = solid_temp.assemble_solid_temperature_lagrange_multiplier(Ts_, uf_, dt)
                  solid_temp.solve_solid_temperature_lagrange_multiplier(a7, LmTs_[0], b7)
              s7 += timer_s7.stop()

**Figure 2: Time loop in vanDANA.py which runs the IB-FSI solver algorithm. All steps are timed using separate timers and the timings are listed in log_info.txt.**

Our flow solver uses the Incremental Pressure Correction Scheme (IPCS) and is solved in step 1,2 and 3. The solid momentum equation is solved in step 5 and the Lagrange multiplier problem is solved in step 6. In the solid equation, one needs to note that we solve for :math:`\Delta \boldsymbol{x}` which is the incremental displacement instead of the current solid position :math:`\boldsymbol{x}` (see Figure 3).


