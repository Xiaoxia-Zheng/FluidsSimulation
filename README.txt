----------------------
3D Fluids Simulation
----------------------
This project is a gridbase 3d fluids simulation. The particles will move based on their initial velocity and interactions with other particles through the background grid. There are five basic functions, one that transfers velocity from the particles to the grid, one that transfers velocity updates from the grid to the particles, then add gravity int the z-direction, enforce slip boundary conditions, compute a pressure field. And subtract the gradient of the pressure to enforce the divergence free condition. Finally advect the particles through the grid's velocity field. Each timestep you should advect the particles, transfer velocities from the particles to the grid, and then transfer updates from the grid to the particle. Using simple Euler integration to advect the particles. A staggered grid for storing the velocity field (i.e. the x-component of the velocity on x-faces, the y-component on y-faces, and the z-component on z-faces). Using trilinear weight functions for transfers between the grid and the particles (and vice versa). The boundary cells have zero velocity and particles never enter these cells.




------------------
Input File Format
------------------
The main input file will be in .json format. This file will specify the main simulation parameters and one particle file. Particle fils
The first line will be the number of particles. Each subsequent line will have six floating point values giving the initial position (x, y, z) and velocity (u, v, w) of the particle.





----------------------
Command Line Arguments
----------------------
Example: ./main input.json output-%05d.part

	"input.json" is the file that we need to input
	"output-%05d.part" is the format that we output the files. 





