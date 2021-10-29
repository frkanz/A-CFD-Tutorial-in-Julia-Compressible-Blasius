# A CFD Tutorial in Julia Compressible Blasius
# A CFD Tutorial in Julia
A CFD Tutorial in Julia: Introduction to Compressible Laminar Boundary-Layer Theory journal's code repository. The codes are developed in the Computational Hypersonics and Aerodynamics Laboratory.

## **Abstract**
A boundary-layer is a thin fluid layer near a solid surface, and viscous effects dominate it. The laminar boundary-layer calculations appear in many aerodynamics problems including, skin friction drag, flow separation, and aerodynamic heating. A student must understand the flow physics and the numerical implementation to conduct successful simulations in advanced undergraduate- and graduate-level fluid dynamics/aerodynamics courses. Numerical simulations require writing computer codes. Therefore, choosing a fast and user-friendly programming language is essential to reduce code development and simulation times. Julia is a new programming language that combines performance and productivity. The present study derived the compressible Blasius equations from Navier-Stokes' equations and numerically solved the resulting equations using the Julia programming language. The fourth-order Runge-Kutta method is used for the numerical discretization, and Newton's iteration method is employed to calculate the missing boundary condition. In addition, Burgers', heat, and compressible Blasius equations are solved both in Julia and MATLAB. The runtime comparison showed that Julia with $for$ loops is 2.8 to 50 times faster than MATLAB in 4 out of 5 cases. We also released the Julia codes on our GitHub page to shorten the learning curve for interested readers.


#### **Oz, F.; Kara, K. A CFD Tutorial in Julia: Introduction to Compressible Laminar Boundary-Layer Theory. Fluids, 2021, **_Under review_**.**

#### **Instructions**

Julia setup files can be downloaded from their website (https://julialang.org/downloads/). The website also includes instructions on how to install Julia on Windows, Linux, and mac operating systems. It is common to use external packages for Julia. In order to do that, Pkg, which is Julia's built-in package manager, can be used. Once Julia is opened, Pkg can be activated with the "]" button in Windows. In Linux, calling "julia" in the terminal will open it. After that "Pkg.add("Pluto")" will trigger the setup process for that package. In here, we used Pluto as an example because, in GitHub, our codes are developed in the Pluto environment. After Pluto is installed. Pluto can be run with "Pluto.run()". This command will open a new tab in the browser which you can run your Julia codes. After that, the "using Pluto" line must be placed to the top of the file. For "Plots" package, the commands will be "Pkg.add("Plots")" and "using Plots". Since the Plots package does not have a GUI, there is not a command called "Plots.run()".

## **Compressible Blasius Equations**
Compressible flow over an infinite flat plate can be represented using the ordinary differential equation (ODE) below:

f''' + 0.5 f f'' = 0  

Boundary Conditions:

f(0) = 0  
f'(0) = 0  
f'(∞) = 1  

#### **Solution Method:**

Let's assume:

p = f(η)  
h = f'(η)  
		
Then, we obtain two ODEs equations and boundary conditions
ODEs:

h'' + p h' = 0  
p' - h = 0  

BCs:  

p(0) = 0  
h(0) = 0  
h(∞) = 1  

The equations given above are solved using Thomas algorithm. The core algorithm is:

A = [ 2/Δη² + p[i]/(2 * Δη) for i =1 : N]  
B = [-4/Δη² for i=1 : N]  
C = [ 2/Δη² - p[i]/(2 * Δη) for i =1 : N]  
D = [ 0 for i=1 : N]  

for i=2 : N-1  
G[i] = - ( C[i] * G[i-1] + D[i] )/(B[i] + C[i] * H[i-1])  
H[i] = - A[i] /(B[i] + C[i] * H[i-1])  
end  

for i=N-1 : -1 : 2  
h[i] = G[i] + H[i] * h[i+1]  
end  

The solution vector is ploted and compared with results from Schlichting's Boundary-layer theory book.

<img src="./001-Blasius_Flow/BlasiusProfile.png" width="50%" height="50%">

## **Please read the full paper for further details.**

Feel free to ask questions!

*Furkan Oz*  
[foz@okstate.edu](foz@okstate.edu)  
  
*Kursat Kara*  
[kursat.kara@okstate.edu](kursat.kara@okstate.edu)  

