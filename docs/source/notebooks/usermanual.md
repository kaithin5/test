---
jupyter:
  deepnote: {}
  deepnote_app_layout: article
  deepnote_notebook_id: aef04ddaa00640aa833f91ec7203a184
  kernelspec:
    display_name: Python 3.8.10 64-bit
    language: python
    name: python3
  language_info:
    codemirror_mode:
      name: ipython
      version: 3
    file_extension: .py
    mimetype: text/x-python
    name: python
    nbconvert_exporter: python
    pygments_lexer: ipython3
    version: 3.8.10
  nbformat: 4
  nbformat_minor: 0
  orig_nbformat: 4
  vscode:
    interpreter:
      hash: 916dbcbb3f70747c44a77c7bcd40155683ae19c65e1c03b4aa3499c5328201f1
---

::: {.cell .markdown cell_id="0593d8eaf22c46c6b20e6a9da78dc6e7" deepnote_app_coordinates="{\"h\":5,\"w\":12,\"x\":0,\"y\":0}" deepnote_cell_type="markdown" tags="[]"}
**Kai Thin, Tejas Patel - Department of Mechanical Engineering, Michigan
State University, East Lansing, MI, 2022**
:::

::: {.cell .markdown cell_id="2d1891bbdfbb42a089fe07ad6bb389de" deepnote_app_coordinates="{\"h\":11,\"w\":12,\"x\":0,\"y\":97}" deepnote_cell_type="markdown"}
vanDANA is a highly efficient FEM Immersed Boundary (IB) based
Flow-thermal FSI solver utilizing the
[FEniCS](https://fenicsproject.org/) library (version 2019.2.0). The
solver is based on the [Distibuted Langrange Multiplier based Fictiious
Domain
method](https://www.sciencedirect.com/science/article/pii/S0021999105000148)
and is extended to deal with [heat
transfer](https://www.sciencedirect.com/science/article/pii/S0021999106000167).
The interpolation of variables is conducted using the smeared
[delta-functions](https://www.sciencedirect.com/science/article/pii/S0021999109004136).
Additionally, the flow solver is incompressible and has the option of
choosing from various stabilization schemes : SUPG, PSPG and Crosswind;
and the structure can be set as either incompressible/compressible.

This manual intends to create a basic understanding of the numerical
algorithm and also provides a detailed explanation of the workflow for
the reader to set up a user-defined problem. Here we expect the reader
to have some basic knowledge of coding partial differential equations in
[FEniCS](https://fenicsproject.org).

`<b>`{=html}Note:`</b>`{=html} vanDANA is licensed under the GNU GPL,
version 3 or any later version and is Copyright (2022) by the authors.
:::

::: {.cell .markdown cell_id="24029acd542840ddb239b4d11ce03cff" deepnote_app_coordinates="{\"h\":5,\"w\":12,\"x\":0,\"y\":7}" deepnote_cell_type="markdown"}
# Directory Tree

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
:::

::: {.cell .markdown cell_id="d4acef93d0ef49369856e8f9969f7735" deepnote_app_coordinates="{\"h\":5,\"w\":12,\"x\":0,\"y\":0}" deepnote_cell_type="markdown" tags="[]"}
**Figure 1: Directory tree for vanDANA IB-FSI code. Note: /src folder
contains sample Gmsh (geo) files for reference.**
:::

::: {.cell .markdown cell_id="a26f6452dd564797b0804839bf25cdee" deepnote_app_coordinates="{\"h\":5,\"w\":12,\"x\":0,\"y\":0}" deepnote_cell_type="text-cell-h1" formattedRanges="[]" is_collapsed="false" tags="[]"}
# Files and modules
:::

::: {.cell .markdown cell_id="64cda579c24e496e86cd7c583b652e43" deepnote_app_coordinates="{\"h\":2,\"w\":8,\"x\":0,\"y\":0}" deepnote_cell_type="text-cell-p" formattedRanges="[{\"fromCodePoint\":9,\"marks\":{\"bold\":true},\"toCodePoint\":22,\"type\":\"marks\"}]" is_collapsed="false" tags="[]"}
vanDANA (IB-FSI branch is most recent) is a python package with main
executable (vanDANA.py) and three submodules (common, utilities,
user_inputs).
:::

::: {.cell .markdown cell_id="151d6d35a83d481488fa6ae32c29b92f" deepnote_app_coordinates="{\"h\":5,\"w\":12,\"x\":0,\"y\":13}" deepnote_cell_type="markdown"}
## vanDANA.py {#vandanapy}

This is the main executable file that pulls information from other
modules (namely: user inputs, utilities, common) and runs the main
function `<b>`{=html}vanDANA_solver (args)`</b>`{=html} which sets the
entire workflow. The user need not require making any changes to
vanDANA.py. Before the time loop, the following preliminary operations
are conducted in order:

1.  Set the current directory and MPI controls.
2.  Read mesh files
3.  Load physical problem classes.
4.  Setup initial and boundary conditions.
5.  Preassemble the flow and temperature problems
6.  Create output files for postprocessing
7.  Read restart files (only if restarting the simulation)

To obtain more control at the command prompt, we use the [argparse
library](https://docs.python.org/3/library/argparse.html) which allows
us to add keywords to the executable, for eg:
`mpirun.mpich -n 64 python3 vanDANA.py -restart=True -T=20 -velocity_degree=1`.
Note that we also provide a separate batch script
([src/fenics.sb](https://github.com/patelte8/vanDANA/blob/IB-FSI/src/fenics.sb))
to submit jobs on the HPC cluster.

Our Immersed Boundary FSI formulation is based on the Distributed
Lagrange Multiplier (DLM) based Fictitious Domain (FD) method and the
physical problem is decomposed into the fluid, solid and Lagrange
multiplier sub-problems. The entire code is non-dimensional, and we use
the Finite Element Method (FEM) for spatial discretization of the
governing equations. For a detailed understanding of the numerical
scheme, we advise the reader to refer to Yu ${ }^{[1]}$ and Yu et al
${ }^{[2]}$.
:::

::: {.cell .markdown cell_id="0e26da99175b4785b86fc8edd28d7066" deepnote_app_coordinates="{\"h\":41,\"w\":12,\"x\":0,\"y\":109}" deepnote_cell_type="markdown" tags="[]"}
**Flow problem: IPCS scheme**

**Step 1 - Velocity prediction**: solve for $u^{*}$

$$
\begin{aligned}
\left\langle\frac{\boldsymbol{u}^{*}-\boldsymbol{u}^{\boldsymbol{\theta}}}{\Delta \tau}, \boldsymbol{v}\right\rangle+ & \left\langle\frac{3 \boldsymbol{u}^{\boldsymbol{\theta}}-\boldsymbol{u}^{\boldsymbol{\theta}-1}}{2} \cdot \nabla \boldsymbol{u}_{\boldsymbol{c k}}^{*}, \boldsymbol{v}\right\rangle \\
& +\left\langle\frac{\nabla \boldsymbol{u}_{\boldsymbol{c k}}^{*}+\nabla \boldsymbol{u}_{\boldsymbol{c k}}^{* T}}{R e}-p^{\theta} \boldsymbol{I}, \nabla \boldsymbol{v}\right\rangle+\left\langle p^{\theta} \cdot \boldsymbol{n}, \boldsymbol{v}\right\rangle_{\partial \boldsymbol{\Omega}}-\left\langle\frac{\widehat{\boldsymbol{g}}}{F r^{2}}, \boldsymbol{v}\right\rangle \\
& +\left\langle\gamma_{S U P G}\left(\boldsymbol{u}^{\boldsymbol{\theta}}\right) P\left(\boldsymbol{u}^{\boldsymbol{\theta}}, \boldsymbol{v}\right), \boldsymbol{R}^{\boldsymbol{\theta}}\right\rangle+\left\langle\gamma_{C W}\left(\boldsymbol{u}^{\boldsymbol{\theta}}\right) \Lambda\left(\boldsymbol{u}^{\boldsymbol{\theta}}, \boldsymbol{u}^{*}\right), \nabla \boldsymbol{v}\right\rangle=\left\langle\lambda^{\boldsymbol{\theta}}, \boldsymbol{v}\right\rangle_{\boldsymbol{P}}
\end{aligned}
$$ `<div style="text-align: right">`{=html} (1) `</div>`{=html}

where, $\boldsymbol{u}^{*}$ is the unknown intermediate flow velocity
and
$\boldsymbol{u}_{\boldsymbol{c k}}^{*}=\frac{\boldsymbol{u}^{*}+\boldsymbol{u}^{\boldsymbol{\theta}}}{2}$
is the Crank-Nicolson velocity. In the above equation, the vector
Lagrange multiplier $\boldsymbol{\lambda}$ serves as a pseudo body force
only in the overlapping solid region $\boldsymbol{P}$ and needs to be
interpolated from the Lagrangian mesh via the smeared delta functions.

**Step 2 - Pressure projection $\&$ velocity correction**: solve for
$p^{\theta+1}$ and correct $u^{\Phi}$.

$$
\begin{gathered}
\left\langle\nabla\left(p^{\theta+1}-p^{\theta}\right), \nabla q\right\rangle+\frac{1}{\Delta \tau}\left\langle\nabla \cdot \boldsymbol{u}^{*}, q\right\rangle+\left\langle\gamma_{P S P G}\left(\boldsymbol{u}^{\theta}\right) \nabla q, \boldsymbol{R}^{\boldsymbol{\theta}}\right\rangle=0 \end{gathered}
$$ `<div style="text-align: right">`{=html} (2) `</div>`{=html}

$$
\begin{gathered}
\left\langle\frac{\boldsymbol{u}^{\boldsymbol{\phi}}-\boldsymbol{u}^{*}}{\Delta \tau}, \boldsymbol{v}\right\rangle+\left\langle\nabla\left(p^{\theta+1}-p^{\theta}\right), \boldsymbol{v}\right\rangle=0
\end{gathered}
$$ `<div style="text-align: right">`{=html} (3) `</div>`{=html}

**Solid problem (compressible)**: solve for $\Delta x^{\theta+1}$ on the
reference configuration.

$$
\begin{aligned}
&\left\langle\rho_{r} \frac{\Delta \boldsymbol{x}^{\boldsymbol{\theta}+\mathbf{1}}}{\Delta \tau^{2}}, \boldsymbol{\omega}\right\rangle_{\boldsymbol{P}_{\mathbf{0}}} \\
&=\left\langle\frac{\boldsymbol{u}^{\boldsymbol{\phi}}}{\Delta \tau}, \boldsymbol{\omega}\right\rangle_{\boldsymbol{P}_{\mathbf{0}}}+\left\langle\rho_{r}-1 \frac{\Delta \boldsymbol{x}^{\boldsymbol{\theta}}}{\Delta \tau^{2}}, \boldsymbol{\omega}\right\rangle_{\boldsymbol{P}_{\mathbf{0}}}-\left\langle\nabla \boldsymbol{\omega}^{T}, \lambda_{0} \ln J \boldsymbol{F}^{-\mathbf{1}}+G\left(\boldsymbol{F}^{\boldsymbol{T}}-\boldsymbol{F}^{-\mathbf{1}}\right)\right\rangle_{\mathbf{P}_{\mathbf{0}}}^{\boldsymbol{\theta}+\mathbf{1}} \\
&+\left\langle\rho_{r}-1 \frac{\widehat{\boldsymbol{g}}}{F r^{2}}, \boldsymbol{\omega}\right\rangle_{\boldsymbol{P}_{\mathbf{0}}}-\left\langle\lambda^{\boldsymbol{\theta}}, \boldsymbol{\omega}\right\rangle_{\boldsymbol{P}}
\end{aligned}
$$ `<div style="text-align: right">`{=html} (4) `</div>`{=html}

**Lagrange multiplier problem**: solve for $\lambda^{\theta+1}$ on the
current configuration.

$$
\left\langle\frac{\Delta \boldsymbol{x}^{\boldsymbol{\theta}+1}}{\Delta \tau^{2}}-\frac{\boldsymbol{u}^{\boldsymbol{\phi}}}{\Delta \tau}, \zeta\right\rangle_{P}=\left\langle\lambda^{\boldsymbol{\theta}+1}-\lambda^{\boldsymbol{\theta}}, \boldsymbol{\zeta}\right\rangle_{P}
$$ `<div style="text-align: right">`{=html} (5) `</div>`{=html}

In the above equations $\langle$,$\rangle$ represents the Euclidean
inner product in the respective domain, $\boldsymbol{u}$ and $p$ are the
fluid velocity and pressure, $\boldsymbol{x}$ is the solid displacement,
$\boldsymbol{F}$ is the deformation gradient tensor and $\rho_{r}$ is
the ratio of solid to fluid densities. Here, $\Omega$ is the overall
fluid domain including the fictitious fluid and $P_{0}$ and $P$ are the
reference and current solid configurations. The above DLM/FD algorithm
is implemented in the time loop in **vanDANA_solver (args)**, see Figure
2 below.
:::

::: {.cell .code cell_id="76e3dda7e8524ed4beff5432f342f66e" deepnote_app_coordinates="{\"h\":73,\"w\":12,\"x\":0,\"y\":151}" deepnote_cell_type="code" output_cleared="true" tags="[]"}
``` python
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
```
:::

::: {.cell .markdown cell_id="37a361d91de94bc59c0a84e46db0ecf2" deepnote_app_coordinates="{\"h\":3,\"w\":12,\"x\":0,\"y\":225}" deepnote_cell_type="markdown" tags="[]"}
`<b>`{=html}Figure 2:$\quad$ Time loop in vanDANA.py which runs the
IB-FSI solver algorithm. All steps are timed using separate timers and
the timings are listed in log_info.txt.`</b>`{=html}
:::

::: {.cell .markdown cell_id="b233c99f7d224552bf3cc7c12e1dc962" deepnote_app_coordinates="{\"h\":5,\"w\":12,\"x\":0,\"y\":229}" deepnote_cell_type="markdown" tags="[]"}
Our flow solver uses the Incremental Pressure Correction Scheme (IPCS)
and is solved in step 1,2 and 3. The solid momentum equation is solved
in step 5 and the Lagrange multiplier problem is solved in step 6. In
the solid equation, one needs to note that we solve for
$\Delta \boldsymbol{x}$ which is the incremental displacement instead of
the current solid position $\boldsymbol{x}$ (see Figure 3).
:::

::: {.cell .markdown cell_id="e0b6bb7ecc5f4a37b7df41857c6044ab" deepnote_app_coordinates="{\"h\":13,\"w\":12,\"x\":0,\"y\":235}" deepnote_cell_type="markdown" tags="[]"}
```{=html}
<center><img src="image-20230123-193321.png" width="500" height="500"></center>
```
:::

::: {.cell .markdown cell_id="3b8746ac51584ced8b121ba05adaf8a6" deepnote_app_coordinates="{\"h\":3,\"w\":12,\"x\":0,\"y\":249}" deepnote_cell_type="markdown" tags="[]"}
`<b>`{=html}Figure 3:$\quad$ Deformation of the solid continuum in
space. Here, $\boldsymbol{X}$ is the reference configuration,
$\boldsymbol{x}$ is the current configuration and
$\Delta \boldsymbol{x}$ is the incremental displacement at any
particular time step.`</b>`{=html}
:::

::: {.cell .markdown cell_id="245874898e5d4ac0b16890b6336bd212" deepnote_app_coordinates="{\"h\":37,\"w\":12,\"x\":0,\"y\":253}" deepnote_cell_type="markdown" tags="[]"}
**Interpolation**: The velocity $\boldsymbol{u}^{\boldsymbol{\phi}}$
calculated at step 3 is interpolated to the solid current configuration
using the smeared delta functions. The same delta function is used to
interpolate the Lagrange multiplier onto the fluid Eulerian mesh.

$$
u^{\phi}(X)=\int_{\Omega} u^{\phi}(x) \delta(x-X) d x  \quad\quad \quad \lambda(x)=\int_{P} \lambda(X) \delta(x-X) d X
$$ `<div style="text-align: right">`{=html} (6) `</div>`{=html}

where
$\delta(x-X)=\frac{1}{h^{3}} \Phi_{4}\left(\frac{x-X}{h}\right) \Phi_{4}\left(\frac{y-Y}{h}\right) \Phi_{4}\left(\frac{z-Z}{h}\right)$
and the 4-point piecewise function is:

$$
\Phi_{4}(r)=\left\{\begin{array}{cc}
\frac{1}{8}\left(3-2|r|-\sqrt{1+4|r|-4 r^{2}}\right) & |r| \leq 1 \\
\frac{1}{8}\left(5-2|r|-\sqrt{-7+12|r|-4 r^{2}}\right) & 1 \leq|r| \leq 2 \\
0 & 2 \leq|r|
\end{array}\right\}
$$ `<div style="text-align: right">`{=html} (7) `</div>`{=html}

Furthermore, our **vanDANA_solver** is also extended to deal with heat
transfer and the mathematical formulation is given as:

**Flow temperature problem**: solve for $\mathrm{T}^{\theta+1}$.

$$
\begin{aligned}
\left\langle\frac{\mathrm{T}^{\theta+1}-\mathrm{T}^{\theta}}{\Delta \tau}, \gamma\right\rangle & +\left\langle\boldsymbol{u}^{\boldsymbol{\theta}+\mathbf{1}} \cdot \nabla \mathrm{T}_{c k}, \gamma\right\rangle+\frac{1}{P e}\left\langle\nabla \mathrm{T}_{c k}, \nabla \gamma\right\rangle-\left\langle\frac{2E c}{R e}\langle\boldsymbol{E}(\boldsymbol{u}),\boldsymbol{\nabla u}\rangle^{\boldsymbol{\theta}+\mathbf{1}}, \gamma\right\rangle \\
& +\left\langle\gamma_{S U P G}\left(\boldsymbol{u}^{\boldsymbol{\theta}+\mathbf{1}}\right) P\left(\boldsymbol{u}^{\boldsymbol{\theta}+\mathbf{1}}, \gamma\right), \boldsymbol{R}^{\boldsymbol{\theta}}\right\rangle+\left\langle\gamma_{C W}\left(\mathrm{~T}^{\theta}\right) \Lambda\left(\boldsymbol{u}^{\boldsymbol{\theta}+\mathbf{1}}, \mathrm{T}^{\theta+1}\right), \nabla \gamma\right\rangle \\
& =\left\langle\lambda_{T}^{\theta}, \gamma\right\rangle_{P}
\end{aligned}
$$ `<div style="text-align: right">`{=html} (8) `</div>`{=html}

set $\mathrm{T}_{s}^{\theta+1}=\mathrm{T}^{\theta+1}$ using
$\mathrm{T}_{s}(\boldsymbol{X})=\int_{\Omega} \mathrm{T}(\boldsymbol{x}) \boldsymbol{\delta}(\boldsymbol{x}-\boldsymbol{X}) \boldsymbol{d} \boldsymbol{x}$.

**Solid temperature Lagrange multiplier problem**: solve for
$\lambda_{T}{ }^{\theta+1}$.

$$
-\left\langle\lambda_{T}^{\theta+1}, \gamma_{s}\right\rangle_{P}=\left(\rho_{r} C_{p r}-1\right)\left\langle\frac{\mathrm{T}_{s}^{\theta+1}-\mathrm{T}_{s}^{\theta}}{\Delta \tau}, \gamma_{s}\right\rangle_{P}+\left\langle\frac{2E c}{R e}\langle\boldsymbol{E}(\boldsymbol{u}),\boldsymbol{\nabla u}\rangle^{\boldsymbol{\theta}+\mathbf{1}}, \gamma_{s}\right\rangle_{P}+\frac{k_{r}-1}{P e}\left\langle\nabla \mathrm{T}_{s, c k}, \nabla \gamma_{s}\right\rangle_{P}
$$ `<div style="text-align: right">`{=html} (9) `</div>`{=html}

where, $\mathrm{T}$ is the fluid temperature and
$\mathrm{T}_{c k}=\frac{\mathrm{T}^{\theta+1}+\mathrm{T}^{\theta}}{2}$
is the Crank-Nicolson discretization, $\mathrm{T}_{s}$ is the solid
temperature, $E c$ is the Eckert number that accounts for heat
generation due to viscous dissipation of the fluid, $C_{p r}$ and
$k_{r}$ is the ratio of heat capacity and thermal conductivity of solid
to fluid and $\lambda_{T}$ is the temperature based Lagrange multiplier.
Note that in Figure 2, the flow temperature problem is computed in step
4 and the solid temperature Lagrange multiplier problem is computed in
step 7. Before ending the time loop, the solver outputs visualization
files (xdmf format), restart files (h5 format), post processing txt
files and also updates [(ALE in
FEniCS)](https://bitbucket.org/fenics-project/dolfin/src/master/dolfin/ale/ALE.cpp)
and smoothens the solid Lagrangian mesh.
:::

::: {.cell .markdown cell_id="a3a3b2d5eeb74f89b32db8c95d55d2d4" deepnote_app_coordinates="{\"h\":5,\"w\":12,\"x\":0,\"y\":43}" deepnote_cell_type="markdown"}
## common

This module is the heart of `<b>`{=html}vanDANA_solver`</b>`{=html} and
mainly comprises of different variational problems.

Step 1, 2 & 3 \| `<b>`{=html}flow_variational_problem.py `</b>`{=html}
Step 4 \|
`<b>`{=html}flow_temperature_variational_problem.py`</b>`{=html} Step 5
\| `<b>`{=html}solid_variational_problem.py`</b>`{=html} Step 6 & 7 \|
`<b>`{=html}lagrange_variational_problem.py`</b>`{=html}

Each of these files contain separate class definitions and subroutines
for pre-assembly, run-time assembly, solvers, and post-processing for
their respective variational problems. Here, we have provided some basic
post-processing subroutines for calculation of drag, lift, Nusset
number, Jacobian and vorticity (for a 2D case problem). Again, the user
need not require making any changes to any of the classes unless they
desire to implement a custom user-defined subroutine.
:::

::: {.cell .markdown cell_id="e082b0b9157e4684bd5612b4f52b8bd8" deepnote_app_coordinates="{\"h\":5,\"w\":12,\"x\":0,\"y\":0}" deepnote_cell_type="markdown" tags="[]"}
### constitutive_eq.py {#constitutive_eqpy}

For the constitutive nature of the material field, we define the fluid
as Newtonian and the solid as neo-Hookean and compressible or
incompressible:

Newtonian fluid:
$\quad \boldsymbol{\sigma}(\boldsymbol{u}, p)=-p \boldsymbol{I}+\frac{2}{R e} \boldsymbol{E}(\boldsymbol{u})$

```{=html}
<div style="text-align: right"> (10a) </div>
```
Compressible Neo-Hookean solid:
$\quad \sigma_{s}(\boldsymbol{x}, \boldsymbol{X})=\lambda_{0} \ln J \boldsymbol{I}+G(\boldsymbol{B}-\boldsymbol{I})$

```{=html}
<div style="text-align: right"> (10b) </div>
```
Incompressible Neo-Hookean solid:
$\boldsymbol{\sigma}_{s}(\boldsymbol{x}, \boldsymbol{X})=-p_{s} \boldsymbol{I}+G(\boldsymbol{B}-\boldsymbol{I})$

```{=html}
<div style="text-align: right"> (10c) </div>
```
where,
$\boldsymbol{E}(\boldsymbol{u})=\frac{1}{2}\left(\nabla \boldsymbol{u}+(\nabla \boldsymbol{u})^{T}\right)$
is the strain rate tensor, $\lambda_{0}$ is the non-dimensional
compressibility, $G$ is the non-dimenionsional shear modulus,
$\boldsymbol{B}=\boldsymbol{F} \boldsymbol{F}^{T}$ is the finger tensor,
$\boldsymbol{F}=\frac{\partial \boldsymbol{x}}{\partial \boldsymbol{X}}$
is the deformation gradient tensor and
$J=\operatorname{det}(\boldsymbol{F})$. In case of an incompressible
solid we also force the Jacobian to be 1 by adding the term
$\langle J-1, \boldsymbol{\omega}\rangle_{P_{0}}$ to the solid problem
(step 5 in Figure 2).
:::

::: {.cell .markdown cell_id="ad827a818e86459cb08037b7053050ea" deepnote_app_coordinates="{\"h\":5,\"w\":12,\"x\":0,\"y\":0}" deepnote_cell_type="markdown" tags="[]"}
### delta_interpolation.py {#delta_interpolationpy}

The smeared delta functions are implemented using a custom C++ code
(which runs in MPI) and the FEniCS backend allows the integration of
custom C++ strings by exposing them as a python object with pybind11.
This is done using the
[compile_cpp_code](https://bitbucket.org/fenics-project/dolfin/src/master/python/dolfin/jit/pybind11jit.py)
functionality in FEniCS.

In such an IB-FSI solver, the delta functions essentially determine the
order of accuracy of the IB method. Due to the inherent nature of the
delta functions, it is obligatory that the fluid grid surrounding the
solid is structured ${ }^{[3]}$ and for good accuracy, we recommend that
the solid mesh size stays between
`<b>`{=html}0.8`<i>`{=html}h`</i>`{=html} -
1.5`<i>`{=html}h`</i>`{=html}`</b>`{=html} at all times where
`<b>`{=html}`<i>`{=html}h`</i>`{=html}`</b>`{=html} is the uniform fluid
grid size ${ }^{[1]}$. Here, note that here the solid grid can be
unstructured and non-uniform. Also, the delta functions utilize
piecewise functions Φ with different support areas and continuity
properties ${ }^{[4]}$ which directly dictate their behavior in terms of
numerical stability and accuracy. We implement nine different piecewise
functions ${ }^{[4]}$ and after extensive testing we find that the
4-point piecewise function (Eq. 7) provides a good balance between
stability and accuracy, hence it is set as our default choice in
user_parameters.py.
:::

::: {.cell .markdown cell_id="3bff20c7ed28406bab9a55f8cd19d292" deepnote_app_coordinates="{\"h\":5,\"w\":12,\"x\":0,\"y\":0}" deepnote_cell_type="markdown" tags="[]"}
### fem_stabilizations.py {#fem_stabilizationspy}

In the Finite-Element realm, for highly advective flows and species
transport (i.e., at high $R e$ and $P e$ ) it is known that the standard
Galerkin formulation can lead to spurious oscillations in the solution.
Hence to prevent this phenomenon, it is necessary to incorporate
appropriate stabilization schemes to ensure numerical stability. In the
`<b>`{=html}vanDANA_solver`</b>`{=html} we use the well-established
Streamline-Upwind Petrov Galerkin (SUPG) stabilization by introducing
the following term to Eqs. (1) and (8):

$$
\int_{\Omega} \gamma_{S U P G}(\boldsymbol{u}) P(\boldsymbol{u}, \boldsymbol{m}) \cdot \boldsymbol{R}^{\boldsymbol{\theta}} d \Omega
$$

```{=html}
<div style="text-align: right"> (11a) </div>
```
where,

$$
\gamma_{S U P G}(\boldsymbol{u})=\alpha\left[\left(\frac{2}{\Delta \tau}\right)^{2}+\left(\frac{2|| \boldsymbol{u} \|}{h}\right)^{2}+9\left(\frac{4}{\mathbb{C} h^{2}}\right)^{2}\right]^{-\frac{1}{2}},
$$

```{=html}
<div style="text-align: right">(11b) </div>
```
$\quad \quad$ In the above equation , $\boldsymbol{R}$ is the residual,
$P(\boldsymbol{u}, \boldsymbol{m})=\boldsymbol{u} \cdot \nabla \boldsymbol{m}$
is the SUPG operator, $h$ denotes the cell diameter, $\mathbb{C}$
denotes the respective control parameter (i.e., $R e$ or $P e$ ) and the
coefficient $\alpha$ is set as $0.85$. Next, we also introduce another
residual based pressure-stabilization term to Eq. (2), namely the
pressure stabilizing Petrov Galerkin (PSPG) scheme that is given as:

$$
\int_{\Omega} \gamma_{P S P G}(\boldsymbol{u}) \nabla q \cdot \boldsymbol{R}^{\boldsymbol{\theta}} d \Omega,
$$

```{=html}
<div style="text-align: right"> (12) </div>
```
where $\gamma_{P S P G}=\gamma_{S U P G}$. The effect of such a PSPG
stabilization scheme is witnessed when circumventing pressure
instabilities that arise due to unstable LBB pairs, for e.g., in case of
equal order discretization (P1-P1) for velocity and pressure. Note that
since we already use the IPCS fractional stepping scheme to decouple the
velocity and pressure, we can still use P1-P1 elements without PSPG
stabilization. However, if using equal order discretization for velocity
and pressure in case of highly advective flows, it is good practice to
use PSPG stabilization since it helps to reduce spurious oscillations in
the pressure field. In our code, both SUPG and the PSPG stabilizations
are implemented using an explicit formulation.

At times in presence of sharp gradients/discontinuities in the flow
field, the SUPG stabilization alone is not sufficient to damp numerical
instabilities. This is mainly because the SUPG operator
$P(\boldsymbol{u}, \boldsymbol{m})$ adds numerical diffusion only in the
streamwise direction. To overcome this, we also introduce a crosswind
stabilizing operator $\Lambda(\boldsymbol{u}, \boldsymbol{\phi})$ to
prevent the localized undershooting/overshooting of the field variables
in the presence of sharp gradients/discontinuities. As such the
following term taken from R. Codina ${ }^{[5]}$, is added to Eqs. (1)
and (8):

$$
\int_{\Omega}\gamma_{C W}(\boldsymbol{\phi}) \Lambda(\boldsymbol{u}, \boldsymbol{\phi}): \boldsymbol{\nabla} \boldsymbol{m} d \Omega
$$

```{=html}
<div style="text-align: right"> (13a) </div>
```
$\quad$where,

$$
\begin{gathered}
\gamma_{C W}(\boldsymbol{\phi})=\frac{1}{2} \beta_{c} h \frac{|| \boldsymbol{R}^{\boldsymbol{\phi}} \|}{\|\nabla \boldsymbol{\phi}\|}
\end{gathered}
$$

```{=html}
<div style="text-align: right"> (13b) </div>
```
In the above Eq\'s,
$\Lambda(\boldsymbol{u}, \boldsymbol{\phi})=\left\{\begin{array}{rr}\boldsymbol{I}-\frac{\boldsymbol{u} \otimes \boldsymbol{u}}{|\boldsymbol{u}|^{2}} & \text { if } \boldsymbol{u} \neq \mathbf{0} \\ \mathbf{0} & \text { if } \boldsymbol{u}=\mathbf{0}\end{array}\right\}\cdot \nabla \boldsymbol{\phi}$,
and $\beta_{c}=$
$\max \left(0, 0.7-\frac{2|| \nabla \phi \|}{\mathbb{C} h|| \boldsymbol{R}^{\phi} \|}\right)$
and $\boldsymbol{\phi}$ is the unknown field variable. In our code, we
implement an implicit formulation for the crosswind stabilization term.
:::

::: {.cell .markdown cell_id="7bbccd4315bf408ea56e915e2c0fcbf5" deepnote_app_coordinates="{\"h\":5,\"w\":12,\"x\":0,\"y\":0}" deepnote_cell_type="markdown" tags="[]"}
### solver_options.py {#solver_optionspy}

The FEniCS language uses several linear algebra backends and the default
choice is set as PETSc. For iterative solvers with the PETSc backend,
FEniCS provides a library of [Krylov solvers and
preconditioners](https://fenicsproject.org/pub/tutorial/html/._ftut1017.html)
to choose from.

In vanDANA, the user control to modify iterative solvers,
preconditioners, absolute/relative tolerances, and monitoring
convergence for all the physical variational problems is provided in
solver_options.py. Here for the flow solver, we use python dictionaries
to initialize solver parameters and the convergence is based on an
absolute tolerance value of 10${ }^{-8}$. For the non-linear solid
problem, we use the newtons method wherein the convergence is determined
based on a relative tolerance value of 10${ }^{-6}$ (see Figure 4). The
$\Delta \boldsymbol{x}$ from the previous time step is provided as
intital guess to the newtons solver and we note that it normally takes
about 2-3 iterations to converge.
:::

::: {.cell .code cell_id="0dade15d21654aaebe12652983e6458a" deepnote_app_coordinates="{\"h\":5,\"w\":12,\"x\":0,\"y\":0}" deepnote_cell_type="code" output_cleared="true" tags="[]"}
``` python
# Solver parameters
krylov_solvers=dict(
    monitor_convergence=False,
    report=False,
    error_on_nonconvergence=True,
    nonzero_initial_guess=True,
    maximum_iterations=300,
    absolute_tolerance=1e-8)

# Solver dictionaries
tentative_velocity_solver=dict(
    solver_type='bicgstab',
    preconditioner_type='jacobi')

solid_displacement_parameters = {"newton_solver":{"linear_solver":solid_momentum_solver['solver_type'], "preconditioner":'hypre_amg', "report":True, \
                                                  "error_on_nonconvergence":True, "absolute_tolerance":1e-15, "relative_tolerance":1e-6, "maximum_iterations":50}}
```
:::

::: {.cell .markdown cell_id="788933fd1f7d459096e9cdd023539f04" deepnote_app_coordinates="{\"h\":5,\"w\":12,\"x\":0,\"y\":0}" deepnote_cell_type="markdown" tags="[]"}
`<b>`{=html}Figure 4: $\quad$ User control to append the choice of
iterative solvers, tolerances and monitoring convergence in
solver_options.py.`</b>`{=html}
:::

::: {.cell .markdown cell_id="fda19b56aed44f02977f3f9e145416d4" deepnote_app_coordinates="{\"h\":5,\"w\":12,\"x\":0,\"y\":0}" deepnote_cell_type="markdown" tags="[]"}
## user_inputs

This python module provides input to the
`<b>`{=html}vanDANA_solver`</b>`{=html} mainly in terms of control
parameters, time step, meshes, initial and boundary conditions for any
given physical problem. Here we allow user control from two files
namely - user_parameters.py, boundary_initial_conditions.py.

In case of any custom inputs to the solver, the user may add special
expressions or subroutines in problem_specific.py.
:::

::: {.cell .markdown cell_id="006bdfb34ed54ce28e655f59b6eb30c5" deepnote_app_coordinates="{\"h\":5,\"w\":12,\"x\":0,\"y\":0}" deepnote_cell_type="markdown" tags="[]"}
### import mesh

The accepted file format for meshes in vanDANA_solver is HDF5. This is
mainly because in FEniCS,
[HDF5file](https://fenicsproject.org/olddocs/dolfin/1.3.0/python/programmers-reference/cpp/io/HDF5File.html)
format is readily compatible for parallel input/output using MPI.

We encourage the user to prepare a tetrahedral mesh on any open-source
platform (for eg. Gmsh, Salome) and export it to an HDF5 format for
input to the vanDANA_solver. The user needs to rename the fluid mesh
file as `<b>`{=html}file_f.h5`</b>`{=html} and the solid mesh file as
`<b>`{=html}file_s.h5`</b>`{=html} before it can be imported by the
vanDANA_solver. Both these files are read-in from the user_inputs
module.
:::

::: {.cell .markdown cell_id="4aef829bcb0e445ab23ec6fb6c88c8ff" deepnote_app_coordinates="{\"h\":5,\"w\":12,\"x\":0,\"y\":37}" deepnote_cell_type="markdown"}
### user_parameters.py {#user_parameterspy}
:::

::: {.cell .code cell_id="8fb496e6d3da4f779dee04013477ff2b" deepnote_app_coordinates="{\"h\":5,\"w\":12,\"x\":0,\"y\":0}" deepnote_cell_type="code" deepnote_to_be_reexecuted="true" output_cleared="true" source_hash="4e046ed6" tags="[]"}
``` python
restart = False									# Restart parameter

# Physics of the problem
# ---------------------------------------------------------------------
problem_physics = dict(
				  solve_temperature = True,		# enter "True" if you want to solve for temperature

				  solve_FSI = True,				# enter "True" if you want to solve for fluid-structure interaction
				  
				  compressible_solid = True,	# enter "True" if compressible: Also remember to specify compressibility (Ld)
				  								# enter "False" if incompressible

				  viscous_dissipation = False,	# Heat release due to viscous gradients 

				  body_force = False,      		# Gravitational force (uniform volumetric force)								 
				)
	
def f_dir(dim):									# Body force direction : -ve y direction (by default)
	
	n = -1*tensors.unit_vector(1, dim) 
	return n

interpolation_fx = 'phi4'						# Delta-function interpolation for FSI problems

# FEM stabilization and constants
# ---------------------------------------------------------------------
stabilization_parameters = dict(	

	# Navier-stokes
	SUPG_NS = False,							# explicit
	PSPG_NS = False,							# explicit		
	crosswind_NS = False,						# implicit
	backflow_NS = False,					

	# Energy-equation
	SUPG_HT = False,							# explicit
	crosswind_HT = False						# implicit
)

alpha = Constant(0.85)                   	  	# SUPG/PSPG stabilization constant 
C_cw = Constant(0.7)                       		# Crosswind stabilization constant (As per R Codina : quadratic elements: 0.35, for linear elements: 0.7)

# Physical parameters    
# ---------------------------------------------------------------------
physical_parameters = dict(
 
	g = 9.81,									# Gravity (m/s2)						  

	# Fluid 
	rho_f = 1,									# Density (kg/m3)
	nu = 1,										# Dynamic viscosity (kg/m.s)
	Spht_f = 1,         						# Specific heat (J/kg.C)
	K_f = 1,									# Thermal conductivity (W/m.C)

	# Solid
	rho_s = 10,									# Density (kg/m3)
	Sm = 0,										# Shear modulus (N/m2)
	Ld = 0,										# Compressibility (N/m2)
	Spht_s = 0.11,								# Specific heat (J/kg.C)
	K_s = 1.2 									# Thermal conductivity (W/m.C)
)

def calc_non_dimensional_solid_properties(g, rho_f, nu, Spht_f, K_f, rho_s, Sm, Ld, Spht_s, K_s, Lsc, Vsc, T0, Tm, Tsc):

	rho = rho_s/rho_f
	Spht = Spht_s/Spht_f
	K = K_s/K_f
	Ld = Ld/(rho_f*Vsc*Vsc)
	Sm = Sm/(rho_f*Vsc*Vsc)
	
	return rho, Spht, K, Ld, Sm

# Characteristic scales
# ---------------------------------------------------------------------
characteristic_scales = dict(
	
	Lsc = 1,			            # m          
	Vsc = 1,	         		    # m/s
	T0 = -1*52,						# lower_temp (C)
	Tm = 37							# higher_temp (c)
)

# Temporal control
# ---------------------------------------------------------------------
time_control = dict(
				 C_no = 0.35, 					# Maximum possible Courant number
			   	 dt = 0.0025,  					# Time-step: constant throughout runtime if adjustable-timestep is "False"
			   	 T = 100,						# Total runtime
			   	 adjustable_timestep = True 	# Calculate time-step using max Courant no. during runtime: used to accelerate temporal solution
			   )

# FEM degree of variables
# ---------------------------------------------------------------------
fem_degree = dict(
				velocity_degree = 2,
				pressure_degree = 1,
				temperature_degree = 2, 
				displacement_degree = 2,
				lagrange_degree = 1
			   )

# Non-dimensional numbers
# ---------------------------------------------------------------------
def calc_non_dimensional_numbers(g, rho_f, nu, Spht_f, K_f, rho_s, Sm, Ld, Spht_s, K_s, Lsc, Vsc, T0, Tm, Tsc):

	Re = rho_f*(Vsc*Lsc)/nu            
	Pr = (Spht_f*nu)/K_f 
	Ec = (Vsc*Vsc)/(Spht_f*(Tm-T0))
	Fr = Vsc/sqrt(g*Lsc) 

	return Re, Pr, Ec, Fr

# Enter "True" if you want to post-process data
# ---------------------------------------------------------------------
post_process = True

# File printing / solid-remeshing control
# ---------------------------------------------------------------------
print_control = dict(
                  a = 40,   # for printing variables and restart files
                  b = 50,  	# for post processing data
                  c = 20, 	# for simulation_wall_time text file
                  d = 5,   	# for remeshing solid current-configuration mesh		
                  e = 20    # for runtime_tsp_courant_no_stats text file	
                )

# If 2D problem?: Do u want to calculate stream function and vorticity! # Note to self: streamfunction is not defined for 3D.
# --------------------------------------------------------------------- 
calc_stream_function = True
```
:::

::: {.cell .markdown cell_id="69ca2527c30c4f048973a2d659a55539" deepnote_app_coordinates="{\"h\":5,\"w\":12,\"x\":0,\"y\":0}" deepnote_cell_type="markdown" tags="[]"}
`<b>`{=html}Figure 5: $\quad$All control parameters are provided as
input to the vanDANA_solver (args) from user_parameters.py. The ones
most frequently used, can also be provided as input in the form of
keywords to the executable (vanDANA.py) using (args).
:::

::: {.cell .markdown cell_id="845d4e33ba2a4f399f977fe822f71e55" deepnote_app_coordinates="{\"h\":5,\"w\":12,\"x\":0,\"y\":0}" deepnote_cell_type="markdown" tags="[]"}
This file defines all the major control parameters (except
boundary/initial conditions) for any generalized problem setup. In
Figure 5, most parameters are grouped using python dictionaries and the
comments are self-sufficient for readers initialization. Note that
dictionaries like `problem_physics` and `stabilization_parameters` can
be initialized easily using booleans. For example, if
`problem_physics['solve_FSI'] == False`, the vanDANA_solver switches to
a standard CFD flow solver with heat transfer.

Another important feature is the `<b>`{=html}restart`</b>`{=html}
parameter, which if \"True\" will flawlessly restart the entire
simulation from the last saved time step. In user_parameters.py, we
allow for the user to input dimensional parametric values (and
corresponding characteristic scales), during runtime the code
automatically calculates the non-dimensional variables using functions
like `calc_non_dimensional_numbers` and
`calc_non_dimenional_solid_properties`.

We also provide a few more important features :

-   The order of finite-element basis functions for all variables can be
    specified using the `fem_degree` dictionary.
-   In the `time_control` dictionary if adjustable_timestep is set as
    True, the code will automatically calculate the time step during
    run-time by limiting the maximum Courant number to
    `time_control['C_no']`. Otherwise if adjustable_timestep is set as
    False, it will use `time_control['dt']` as a constant time step.
-   The `print_control` dictionary permits control to specify different
    counters used for printing, mesh smoothing etc. Note that if the
    `<b>`{=html}post_process`</b>`{=html} flag is set as False, the code
    will override the `print_control['b']` value and disable the output
    of post_processing text files.
-   The `calc_stream_function` flag if True, enables calculation of
    vorticity and stream function only for a 2D problem setup. This is
    mainly because the steam function is defined only in 2D and we use
    the following equation:

$$
\boldsymbol{\nabla} ^2 \psi = -\omega
$$ `<div style="text-align: right">`{=html} (14) `</div>`{=html}

-   The flag `problem_physics['body_force'] == True`, the function
    `<b>`{=html}f_dir(dim)`</b>`{=html} and the value of
    `physical_paramters['g']` are used to set the activation, direction
    and magnitude of volumetric gravitational force in the overall FSI
    domain.
-   The effect of stabilization schemes can also be varied using the two
    stabilization constants `<b>`{=html}alpha`</b>`{=html} and
    `<b>`{=html}C_cw`</b>`{=html}.

`<b>`{=html}Important:`</b>`{=html} One of the most central advantage of
our code is the ease of definition in user_parameters.py and the
flexiblity to update dictionaries using keywords from the terminal,
which allows for easy restart and testing of the vanDANA_solver.
:::

::: {.cell .markdown cell_id="cc7883967b5e41b3b4090a7f7ee608b8" deepnote_app_coordinates="{\"h\":5,\"w\":12,\"x\":0,\"y\":19}" deepnote_cell_type="markdown"}
### boundary_initial_conditions.py {#boundary_initial_conditionspy}
:::

::: {.cell .code cell_id="a6b6604c770a4a8b842fb782d8ff89f1" deepnote_app_coordinates="{\"h\":5,\"w\":12,\"x\":0,\"y\":25}" deepnote_cell_type="code" deepnote_to_be_reexecuted="false" execution_millis="370" execution_start="1673387555049" output_cleared="true" source_hash="47c312c7"}
``` python
class PeriodicDomain(SubDomain):

    def inside(self, x, on_boundary):
        return bool(x[2] < DOLFIN_EPS and x[2] > -DOLFIN_EPS and on_boundary)

    def map(self, x, y):
        y[0] = x[0]
        y[1] = x[1]
        y[2] = x[2] - 1.0

constrained_domain = PeriodicDomain()    # None

parabolic_profile = Expression('6.0*x[1]*(4.1 - x[1])/(4.1*4.1)', degree=2)

class Point_pressure(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 4.2) and near(x[1], 5.) and near(x[2], 2.)

# Boundary conditions
def fluid_create_boundary_conditions(fluid_mesh, inflow, **V):

	boundaries = fluid_mesh.get_mesh_boundaries()

	# velocity
	bcu_left_x = DirichletBC(V['fluid'][0], parabolic_profile, boundaries, 1)
	bcu_bottom_x = DirichletBC(V['fluid'][0], Constant(0), boundaries, 2)
	bcu_top_x = DirichletBC(V['fluid'][0], Constant(0), boundaries, 4)
	bcu_x = [bcu_left_x, bcu_bottom_x, bcu_top_x]

	bcu_left_y = DirichletBC(V['fluid'][0], Constant(0), boundaries, 1)
	bcu_bottom_y = DirichletBC(V['fluid'][0], Constant(0), boundaries, 2)
	bcu_top_y = DirichletBC(V['fluid'][0], Constant(0), boundaries, 4)
	bcu_y = [bcu_left_y, bcu_bottom_y, bcu_top_y]

    bcu_left_z = DirichletBC(V['fluid'][0], Constant(0), boundaries, 1)
	bcu_bottom_z = DirichletBC(V['fluid'][0], Constant(0), boundaries, 2)
	bcu_z = [bcu_left_z, bcu_bottom_z]

	bcu = [bcu_x, bcu_y, bcu_z]

	# pressure
	bcp_right = DirichletBC(V['fluid'][1], Constant(0), boundaries, 3)
	bcp = [bcp_right]

	# Streamfunction
	wall  = 'on_boundary'
	bcPSI = DirichletBC(V['fluid'][1], 0, wall)

	bcs = dict(velocity = bcu, pressure = bcp, streamfunction = bcPSI)

	if problem_physics['solve_temperature'] == True:
		# temperature
		bcT_left = DirichletBC(V['fluid_temp'][0], Constant(1), boundaries, 1)
		bcT_top = DirichletBC(V['fluid_temp'][0], Constant(0), boundaries, 4)
		bcT = [bcT_left, bcT_top]
		
		bcs.update(temperature = bcT)
			
	return bcs


def solid_create_boundary_conditions(solid_mesh_R, compressible_solid, dt, **V):

	boundaries = solid_mesh_R.get_mesh_boundaries()
	
	# Note to self: Boundary conditions are for incremental displacement (delta D)

	# Solid
	if compressible_solid == False:
		bcx_cylinder = DirichletBC(V['solid'][1].sub(0), Constant((0, 0, 0)), boundaries, 1)
	elif compressible_solid == True:
	    bcx_cylinder = DirichletBC(V['solid'][0], Constant((0, 0, 0)), boundaries, 1)

	bcx = [bcx_cylinder]  
	return bcx    


# Initial conditions
def fluid_create_initial_conditions(u_, p_, T_):

	# Velocity / pressure
	for i in range(3):
		u_[i][0].vector()[:] = 0.0
		u_[i][1].vector()[:] = 0.0
		u_[i][2].vector()[:] = 0.0
		p_[i].vector()[:] = 0.0

	# Temperature
	for i in range(3):
		T_[i].vector()[:] = 0.0
	

def solid_create_initial_conditions(Dp_, mix, dt):
	
	# Solid pressure (only defined for incompressible solid)
	assign(mix.sub(1), interpolate(Constant(0), mix.sub(1).function_space().collapse()))

	# Cumulative displacement
	Dp_[0].vector()[:] = 0.0 

	# Incremental displacement (delta D)
	Dp_[1].vector()[:] = 0.0 # V_init*dt
	Dp_[2].vector()[:] = 0.0 # V_init*dt
	assign(mix.sub(0), interpolate(Expression(('0.0', '0.0', '0.0'), degree = 2), mix.sub(0).function_space().collapse()))
```
:::

::: {.cell .markdown cell_id="773f978bfdca408c8e492b0bb4decf5d" deepnote_app_coordinates="{\"h\":5,\"w\":12,\"x\":0,\"y\":0}" deepnote_cell_type="markdown" tags="[]"}
`<b>`{=html}Figure 6:$\quad$The boundary and initial conditons are
provided as input to the vanDANA_solver from
boundary_initial_conditions.py.
:::

::: {.cell .markdown cell_id="0fdb0ce3518a4b268198c745a428bf20" deepnote_app_coordinates="{\"h\":5,\"w\":12,\"x\":0,\"y\":0}" deepnote_cell_type="markdown" tags="[]"}
To provide boundary conditions for any CFD setup, we need to flag the
boundaries and regions/subdomains in the physical domain or meshfile.
The flaging is generally done in advance during mesh preparation and
hereafter can used in boundary_intial_conditions.py. In Figure 6,
separate functions are provided to setup boundary/initial conditions for
the fluid and solid domain.

For the solid, the boundary conditions are provided for
$\Delta \boldsymbol{x}$, instead of $\boldsymbol{x}$, where
$\Delta \boldsymbol{x}$ is defined as a vector on the appropriate
function space - for eg:
`DirichletBC(V['solid'][1].sub(0), Constant((0, 0, 0)), boundaries, 1)`.
Otherwise for the fluid, boundary condtions are provided component wise
as scalers for each direction - for eg:
`bcu_bottom_x = DirichletBC(V['fluid'][0], Constant(0), boundaries, 2)`.
The same applies to setting up initial conditions. This is mainly
because the flow solver
([common/flow_variational_problem.py](https://github.com/patelte8/vanDANA/blob/IB-FSI/common/flow_variational_problem.py))
is split to solve component wise, so as to achieve better efficiency. If
you have multiple boundary conditions, note that you still have to
combine them into one list, for eg:
`bcu_x = [bcu_left_x, bcu_bottom_x, bcu_top_x]`. Also, if your physics
requires the need for `<b>`{=html}periodic boundary
conditions`</b>`{=html}, this can be setup using the class
[PeriodicDomain(Subdomain)](https://fenicsproject.org/olddocs/dolfin/1.4.0/python/demo/documented/periodic/python/documentation.html)
in FEniCS.
:::

::: {.cell .markdown cell_id="4636318c93554e4c97e5375f3800ed8c" deepnote_app_coordinates="{\"h\":5,\"w\":12,\"x\":0,\"y\":0}" deepnote_cell_type="markdown" tags="[]"}
## utilities

The utilities module deals with the handling of background functions for
the code and comprises of three files : utils.py, read.py and write.py.
Almost all of the functionality for read, write, restart, MPI, counters,
memory usage and timing of different modules is initialized and carried
out using this module. For regular use of the vanDANA solver, the user
need not make any changes to this module.
:::

::: {.cell .markdown cell_id="575e21946744443bac2cb7ae6e76130c" deepnote_app_coordinates="{\"h\":5,\"w\":12,\"x\":0,\"y\":0}" deepnote_cell_type="markdown" tags="[]"}
# Case Setup

To setup your own physical problem,

1.  Generate the fluid and immersed solid mesh as HDF5 format and move
    it inside the user_inputs module. Remember to flag the essential
    boundaries and domains to apply boundary conditions. We recommend
    using [Gmsh](https://gmsh.info) to generate the mesh files.

2.  Setup your desired control parameters in user_parameters.py
    including total run time and time-step.

3.  Setup your required boundary and initial conditions in
    boundary_initial_conditions.py. If need arises, add custom
    expressions/subroutines in problem_specific.py or else leave it
    empty.

4.  If you wish to control the iterative solver settings (for eg:
    tolerances, monitoring convergence) for any particular variational
    problem, change the default values in solver_options.py.

5.  To calculate post-processing quantities at any flagged boundary,
    navigate to `<b>`{=html}post_process_data`</b>`{=html} subroutine in
    the set of common/\... variational_problem.py files and write
    desired code to export data from the variational problem. If you do
    not wish to output any post processing data, set the post_process
    flag in user_parameters.py as False.

6.  Navigate to the root folder and run the executable (vanDANA.py).

Alternatively, one can also append important control parameters from the
terminal. You can do so by adding arguments to the executable. This will
allow efficient testing of the solver, for eg. :
`mpirun.mpich -n 64 python3 vanDANA.py -T=0.05 -a=1000 -velocity_degree=1`
:::

::: {.cell .markdown cell_id="9f4b55f2d01b43709b87682035dff5c2" deepnote_app_coordinates="{\"h\":5,\"w\":12,\"x\":0,\"y\":91}" deepnote_cell_type="markdown"}
# Running the code

A singularity build file is provided that will install necessary
libraries & setup the environment to run the code.

1.  Install singularity (version 3.5 or newer) by following the
    instruction in
    [here](https://docs.sylabs.io/guides/3.6/admin-guide/installation.html).

2.  Build a singularity container using the build file
    [(src/fenics_2019_dev)](https://github.com/patelte8/vanDANA/blob/IB-FSI/src/fenics_2019_dev)
    by

    `sudo singularity build <container_name>.img fenics_2019_dev`

3.  Once the container is built, you can launch the singularity
    container by

    `singularity run <container_name>.img`

4.  Go to the root folder and run `vanDANA.py` using

    `mpirun.mpich -n <#processors> python3 vanDANA.py`
:::

::: {.cell .markdown cell_id="5e95bb9c9b7c477b992bdf809a40cfff" deepnote_app_coordinates="{\"h\":5,\"w\":12,\"x\":0,\"y\":49}" deepnote_cell_type="markdown"}
# Turek, Hron benchmark (FSI2)

We will now illustrate the implemention of a 2D classical benchmark
(Turek and Hron ${ }^{[6]}$) using our vanDANA_solver. The fluid
domain - $\Omega = [0, 11] \times[0, 4.1]$ comprises of a uniform mesh
$700\times260$
([src/mesh.py](https://github.com/patelte8/vanDANA/blob/IB-FSI/src/mesh.py))
and the solid is an elastic slender flag attached to the back of a rigid
cylinder with center as $(2,2)$ and diameter as $D=1$. At beginning -
reference configuration $P_0$, the flag is symmetric with respect to the
cylinder center and its length (along the centerline) and thinkness are
$l=3.5$ and $h=0.2$ respectively. The solid mesh is constructed such
that is uniform inside the flag and has 12 mesh cells along its
thickness
([src/flag.geo](https://github.com/patelte8/vanDANA/blob/IB-FSI/src/flag.geo)).
The boundary conditions for the top and bottom boundary are no-slip,
left boundary is a constant inlet velocity profile of \$U=\\frac{6U_0
y(H-y)}{H\^2}; U_0=1 \$ and the right boundary is outlet with pressure
set as $0$ and $\nabla\boldsymbol{u}\cdot\boldsymbol{n}=0$. The
characteristic length and velocity scales are $D$ and $U_0$ respectively
and the $Re$ is set as $100$. For the elastic flag - the density ratio
is $ρ_r=10$, non-dimensional shear modulus is $G = 500$, and the solid
is considered compressible neo-Hookean with $λ_0=2000$.
:::

::: {.cell .markdown cell_id="e8125c50b536447e9ea4dcb61a67f54f" deepnote_app_coordinates="{\"h\":5,\"w\":12,\"x\":0,\"y\":0}" deepnote_cell_type="markdown" tags="[]"}
## Results and Post Processing
:::

::: {.cell .markdown cell_id="27c8f1b341d040baa1cee5ce7b88a775" deepnote_app_coordinates="{\"h\":5,\"w\":12,\"x\":0,\"y\":0}" deepnote_cell_type="markdown" tags="[]"}
```{=html}
<center><img src="/work/turek_benchmark.gif" width="500" height="500"></center>
```
:::

::: {.cell .markdown cell_id="ff0faaa63bbd4f07bc4a9a95b616430e" deepnote_app_coordinates="{\"h\":5,\"w\":12,\"x\":0,\"y\":0}" deepnote_cell_type="text-cell-h1" formattedRanges="[]" is_collapsed="false" tags="[]"}
# References
:::

::: {.cell .markdown cell_id="4cf201b045e241cc8f1bad1795bd0399" deepnote_app_coordinates="{\"h\":5,\"w\":12,\"x\":0,\"y\":0}" deepnote_cell_type="markdown" tags="[]"}
\[1\] Yu, Zhaosheng. \"A DLM/FD method for fluid/flexible-body
interactions.\" Journal of computational physics 207, no. 1 (2005): 1-27

\[2\] Yu, Zhaosheng, Xueming Shao, and Anthony Wachs. \"A fictitious
domain method for particulate flows with heat transfer.\" Journal of
Computational Physics 217, no. 2 (2006): 424-452

\[3\] Wang, Xingshi, and Lucy T. Zhang. \"Interpolation functions in the
immersed boundary and finite element methods.\" Computational Mechanics
45, no. 4 (2010): 321-334

\[4\] Yang, Xiaolei, Xing Zhang, Zhilin Li, and Guo-Wei He. \"A
smoothing technique for discrete delta functions with application to
immersed boundary method in moving boundary simulations.\" Journal of
Computational Physics 228, no. 20 (2009): 7821-7836

\[5\] R. Codina, A discontinuity-capturing crosswind-dissipation for the
finite element solution of the convection-diffusion equation, Computer
Methods in Applied Mechanics and Engineering. 110 (1993) 325--342

\[6\] Turek, Stefan, and Jaroslav Hron. Proposal for numerical
benchmarking of fluid-structure interaction between an elastic object
and laminar incompressible flow. Springer Berlin Heidelberg, 2006
:::

::: {.cell .markdown created_in_deepnote_cell="true" deepnote_cell_type="markdown" tags="[]"}
`<a style='text-decoration:none;line-height:16px;display:flex;color:#5B5B62;padding:10px;justify-content:end;' href='https://deepnote.com?utm_source=created-in-deepnote-cell&projectId=97fcbefa-f470-4c41-abd7-d08be048be4a' target="_blank">`{=html}
`<img alt='Created in deepnote.com' style='display:inline;max-height:16px;margin:0px;margin-right:7.5px;' src='data:image/svg+xml;base64,PD94bWwgdmVyc2lvbj0iMS4wIiBlbmNvZGluZz0iVVRGLTgiPz4KPHN2ZyB3aWR0aD0iODBweCIgaGVpZ2h0PSI4MHB4IiB2aWV3Qm94PSIwIDAgODAgODAiIHZlcnNpb249IjEuMSIgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzIwMDAvc3ZnIiB4bWxuczp4bGluaz0iaHR0cDovL3d3dy53My5vcmcvMTk5OS94bGluayI+CiAgICA8IS0tIEdlbmVyYXRvcjogU2tldGNoIDU0LjEgKDc2NDkwKSAtIGh0dHBzOi8vc2tldGNoYXBwLmNvbSAtLT4KICAgIDx0aXRsZT5Hcm91cCAzPC90aXRsZT4KICAgIDxkZXNjPkNyZWF0ZWQgd2l0aCBTa2V0Y2guPC9kZXNjPgogICAgPGcgaWQ9IkxhbmRpbmciIHN0cm9rZT0ibm9uZSIgc3Ryb2tlLXdpZHRoPSIxIiBmaWxsPSJub25lIiBmaWxsLXJ1bGU9ImV2ZW5vZGQiPgogICAgICAgIDxnIGlkPSJBcnRib2FyZCIgdHJhbnNmb3JtPSJ0cmFuc2xhdGUoLTEyMzUuMDAwMDAwLCAtNzkuMDAwMDAwKSI+CiAgICAgICAgICAgIDxnIGlkPSJHcm91cC0zIiB0cmFuc2Zvcm09InRyYW5zbGF0ZSgxMjM1LjAwMDAwMCwgNzkuMDAwMDAwKSI+CiAgICAgICAgICAgICAgICA8cG9seWdvbiBpZD0iUGF0aC0yMCIgZmlsbD0iIzAyNjVCNCIgcG9pbnRzPSIyLjM3NjIzNzYyIDgwIDM4LjA0NzY2NjcgODAgNTcuODIxNzgyMiA3My44MDU3NTkyIDU3LjgyMTc4MjIgMzIuNzU5MjczOSAzOS4xNDAyMjc4IDMxLjY4MzE2ODMiPjwvcG9seWdvbj4KICAgICAgICAgICAgICAgIDxwYXRoIGQ9Ik0zNS4wMDc3MTgsODAgQzQyLjkwNjIwMDcsNzYuNDU0OTM1OCA0Ny41NjQ5MTY3LDcxLjU0MjI2NzEgNDguOTgzODY2LDY1LjI2MTk5MzkgQzUxLjExMjI4OTksNTUuODQxNTg0MiA0MS42NzcxNzk1LDQ5LjIxMjIyODQgMjUuNjIzOTg0Niw0OS4yMTIyMjg0IEMyNS40ODQ5Mjg5LDQ5LjEyNjg0NDggMjkuODI2MTI5Niw0My4yODM4MjQ4IDM4LjY0NzU4NjksMzEuNjgzMTY4MyBMNzIuODcxMjg3MSwzMi41NTQ0MjUgTDY1LjI4MDk3Myw2Ny42NzYzNDIxIEw1MS4xMTIyODk5LDc3LjM3NjE0NCBMMzUuMDA3NzE4LDgwIFoiIGlkPSJQYXRoLTIyIiBmaWxsPSIjMDAyODY4Ij48L3BhdGg+CiAgICAgICAgICAgICAgICA8cGF0aCBkPSJNMCwzNy43MzA0NDA1IEwyNy4xMTQ1MzcsMC4yNTcxMTE0MzYgQzYyLjM3MTUxMjMsLTEuOTkwNzE3MDEgODAsMTAuNTAwMzkyNyA4MCwzNy43MzA0NDA1IEM4MCw2NC45NjA0ODgyIDY0Ljc3NjUwMzgsNzkuMDUwMzQxNCAzNC4zMjk1MTEzLDgwIEM0Ny4wNTUzNDg5LDc3LjU2NzA4MDggNTMuNDE4MjY3Nyw3MC4zMTM2MTAzIDUzLjQxODI2NzcsNTguMjM5NTg4NSBDNTMuNDE4MjY3Nyw0MC4xMjg1NTU3IDM2LjMwMzk1NDQsMzcuNzMwNDQwNSAyNS4yMjc0MTcsMzcuNzMwNDQwNSBDMTcuODQzMDU4NiwzNy43MzA0NDA1IDkuNDMzOTE5NjYsMzcuNzMwNDQwNSAwLDM3LjczMDQ0MDUgWiIgaWQ9IlBhdGgtMTkiIGZpbGw9IiMzNzkzRUYiPjwvcGF0aD4KICAgICAgICAgICAgPC9nPgogICAgICAgIDwvZz4KICAgIDwvZz4KPC9zdmc+' >`{=html}
`</img>`{=html} Created in
`<span style='font-weight:600;margin-left:4px;'>`{=html}Deepnote`</span>`{=html}`</a>`{=html}
:::
