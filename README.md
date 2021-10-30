# A CFD Tutorial in Julia: Compressible Blasius
A CFD Tutorial in Julia: Introduction to Compressible Laminar Boundary-Layer Theory journal's code repository. The codes are developed in the Computational Hypersonics and Aerodynamics Laboratory.

## **Abstract**
A boundary-layer is a thin fluid layer near a solid surface, and viscous effects dominate it. The laminar boundary-layer calculations appear in many aerodynamics problems including, skin friction drag, flow separation, and aerodynamic heating. A student must understand the flow physics and the numerical implementation to conduct successful simulations in advanced undergraduate- and graduate-level fluid dynamics/aerodynamics courses. Numerical simulations require writing computer codes. Therefore, choosing a fast and user-friendly programming language is essential to reduce code development and simulation times. Julia is a new programming language that combines performance and productivity. The present study derived the compressible Blasius equations from Navier-Stokes' equations and numerically solved the resulting equations using the Julia programming language. The fourth-order Runge-Kutta method is used for the numerical discretization, and Newton's iteration method is employed to calculate the missing boundary condition. In addition, Burgers', heat, and compressible Blasius equations are solved both in Julia and MATLAB. The runtime comparison showed that Julia with $for$ loops is 2.8 to 50 times faster than MATLAB in 4 out of 5 cases. We also released the Julia codes on our GitHub page to shorten the learning curve for interested readers.


#### **Oz, F.; Kara, K. A CFD Tutorial in Julia: Introduction to Compressible Laminar Boundary-Layer Theory. Fluids, 2021, **_Under review_**.**

#### **Instructions**

Julia setup files can be downloaded from their website (https://julialang.org/downloads/). The website also includes instructions on how to install Julia on Windows, Linux, and mac operating systems. It is common to use external packages for Julia. In order to do that, Pkg, which is Julia's built-in package manager, can be used. Once Julia is opened, Pkg can be activated with the "]" button in Windows. In Linux, calling "julia" in the terminal will open it. After that "Pkg.add("Pluto")" will trigger the setup process for that package. In here, we used Pluto as an example because, in GitHub, our codes are developed in the Pluto environment. After Pluto is installed. Pluto can be run with "Pluto.run()". This command will open a new tab in the browser which you can run your Julia codes. After that, the "using Pluto" line must be placed to the top of the file. For "Plots" package, the commands will be "Pkg.add("Plots")" and "using Plots". Since the Plots package does not have a GUI, there is not a command called "Plots.run()".

## **Compressible Blasius Equations**
Boundary-layer velocity and temperature profiles on the flat plate can be projected onto single profile wich is self-similar profile. It can be represented using the ordinary differential equations (ODEs) below (Equations are not visible in dark theme. Please check the equations from the main solver or from the paper if necessary.):

<a href="https://www.codecogs.com/eqnedit.php?latex=(cf'')'&plus;ff''&space;=0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?(cf'')'&plus;ff''&space;=0" title="(cf'')'+ff'' =0" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=(a_1g'&plus;a_2f'f'')'&plus;fg'=0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?(a_1g'&plus;a_2f'f'')'&plus;fg'=0" title="(a_1g'+a_2f'f'')'+fg'=0" /></a>
    
where 

<a href="https://www.codecogs.com/eqnedit.php?latex=f'=\frac{u}{u_e}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?f'=\frac{u}{u_e}" title="f'=\frac{u}{u_e}" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=c=\frac{\rho&space;\mu}{\rho_e&space;\mu_e}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?c=\frac{\rho&space;\mu}{\rho_e&space;\mu_e}" title="c=\frac{\rho \mu}{\rho_e \mu_e}" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=g=\frac{H}{H_e}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?g=\frac{H}{H_e}" title="g=\frac{H}{H_e}" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=a_1=\frac{c}{\sigma}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?a_1=\frac{c}{\sigma}" title="a_1=\frac{c}{\sigma}" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=a_2=\frac{(\gamma-1)M^2}{1&plus;(\frac{\gamma-1}{2})M^2}\left(1-\frac{1}{\sigma}\right)c" target="_blank"><img src="https://latex.codecogs.com/gif.latex?a_2=\frac{(\gamma-1)M^2}{1&plus;(\frac{\gamma-1}{2})M^2}\left(1-\frac{1}{\sigma}\right)c" title="a_2=\frac{(\gamma-1)M^2}{1+(\frac{\gamma-1}{2})M^2}\left(1-\frac{1}{\sigma}\right)c" /></a>

and H is the enthalpy, γ is the ratio of specific heats, M is the edge Mach number, and σ is the Prandtl number. σ and M can be defined as

<a href="https://www.codecogs.com/eqnedit.php?latex=M=\frac{u_e}{\sqrt{\gamma&space;\mathfrak{R}T_e}}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?M=\frac{u_e}{\sqrt{\gamma&space;\mathfrak{R}T_e}}" title="M=\frac{u_e}{\sqrt{\gamma \mathfrak{R}T_e}}" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=\sigma=\frac{\mu&space;c_p}{k}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\sigma=\frac{\mu&space;c_p}{k}" title="\sigma=\frac{\mu c_p}{k}" /></a>
    
In this code, σ is assumed as 0.72. The viscosity μ is a function of T and it is calculated as

<a href="https://www.codecogs.com/eqnedit.php?latex=\mu&space;=&space;c_1\frac{T^{3/2}}{(T&plus;c_2)}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mu&space;=&space;c_1\frac{T^{3/2}}{(T&plus;c_2)}" title="\mu = c_1\frac{T^{3/2}}{(T+c_2)}" /></a>

c₂ is 110.4 Kelvin. c₁ is disappearing on the nondimensionalizing process. The boundary conditions for the system of ODEs are
    
<a href="https://www.codecogs.com/eqnedit.php?latex=y=0;\hspace{12pt}f=f'=0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?y=0;\hspace{12pt}f=f'=0" title="y=0;\hspace{12pt}f=f'=0" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=y\rightarrow&space;\infty;\hspace{12pt}f',g&space;\rightarrow&space;0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?y\rightarrow&space;\infty;\hspace{12pt}f',g&space;\rightarrow&space;0" title="y\rightarrow \infty;\hspace{12pt}f',g \rightarrow 0" /></a>

The resultant equations along with the boundary conditions are solved with the Runge-Kutta scheme with Newton's iteration method for the missing boundary condition.

The solution vector is ploted and compared with results from Schlichting's Boundary-layer theory book.

<img src="./001-Compressible_Blasius/M-4_5-T-61-Vel.png" width="50%" height="50%">
<img src="./001-Compressible_Blasius/M-4_5-T-61-Temp.png" width="50%" height="50%">

Details of RK:
Numerical Recipes, Cambridge

Details of Similarity solution formulation:
Boundary-Layer Theory, 7ᵗʰ edition, Schlichting

Feel free to ask questions!

## **Please read the full paper for further details.**

Feel free to ask questions!

*Furkan Oz*  
[foz@okstate.edu](foz@okstate.edu)  
  
*Kursat Kara*  
[kursat.kara@okstate.edu](kursat.kara@okstate.edu)  

