----------------------
CMSC691 Assignment 6
----------------------
Xiaoxia Zheng / CE83376


----------------------
Command Line Arguments
----------------------
Example: ./main input.json output-%05d.part

	"input.json" is the file that we need to input
	"output-%05d.part" is the format that we output the files. 

-----------------
Project sturcture
-----------------
1. Read all input files.
2. Do advection of all particles and update their new position.
3. Convert particle velocity to the grid.
4. Add gravity in the z-direction, enforce slip boundary conditions.
5. Enforce the divergence free condition. I set the threshole as 1e-6.
6. Convert grid velocity back to the particle.
7. Output files.


------------------------------------    
Problems I met and resource I use
------------------------------------
1. Materials that the professor provided help a lot.
2. I also read a CMU slice about grid base fluid animation, it helps me too.
3. And there is a video helps me , the link is: https://www.youtube.com/watch?v=6C0MYfjmtoo&t=3106s


