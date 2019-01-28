# 1D Heat Diffusion
This simulation of diffusion of heat in one dimension is done using MPI. Optimal utilization of the CPU cores is ensured by dividing the work equally among the different processor ranks as much as possible.

The heat diffusion equation which for 1 dimension is :


<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{du}{dt}=\alpha\frac{d^2u}{dx^2}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{du}{dt}=\alpha\frac{d^2u}{dx^2}" title="\frac{du}{dt}=\alpha\frac{d^2u}{dx^2}" /></a>

For this project, the α, heat constant, is one. Another assumpution is that: u is constant
at both endpoints. Using a stepsize of h and k

<a href="https://www.codecogs.com/eqnedit.php?latex=x_i=x_0&plus;ih" target="_blank"><img src="https://latex.codecogs.com/gif.latex?x_i=x_0&plus;ih" title="x_i=x_0+ih" /></a>

<a href="https://www.codecogs.com/eqnedit.php?latex=t_i=t_0&plus;ih" target="_blank"><img src="https://latex.codecogs.com/gif.latex?t_i=t_0&plus;ih" title="t_i=t_0+ih" /></a>

discritizing and solving this equation for u at time n+1 and position j yields :

<a href="https://www.codecogs.com/eqnedit.php?latex=u_j^{n&plus;1}=(1-2r)u_j^n&plus;ru_{j-1}^n&plus;ru_{j&plus;1}^n" target="_blank"><img src="https://latex.codecogs.com/gif.latex?u_j^{n&plus;1}=(1-2r)u_j^n&plus;ru_{j-1}^n&plus;ru_{j&plus;1}^n" title="u_j^{n+1}=(1-2r)u_j^n+ru_{j-1}^n+ru_{j+1}^n" /></a>   
with 

<a href="https://www.codecogs.com/eqnedit.php?latex=r=\frac{k}{h^2}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?r=\frac{k}{h^2}" title="r=\frac{k}{h^2}" /></a>

doubles are used to calculate and store the final values.
