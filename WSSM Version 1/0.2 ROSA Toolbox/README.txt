Code Map:

Start with A2_run_ROSA:
	- ROSAshape gives you initial fire line data.
	- ROSAmain is the main code.
	- ROSAplot plots the outputted fire spread.

ROSAmain:
	- ROSAtstep computes the next time step.
	- ROSAmerger checks for merging and merges fire lines.
	- ROSAfimage outputs image of fire spread, if selected.

ROSAtstep:
	- firestep computes delta z
	- timestepAAA computes the pyrowind effect using multiply connected AAA-LS
	- ROSApimage outputs image of phi contour plot, if selected.
	- RKfun computes Runge Kutta timestepping e.g. z_{t+1} = z_t + delta z * delta t

----

External functions needed:

LineCurvature2D: calculate curvature of a 2D curve.
Ref: Dirk-Jan Kroon (2024). 2D Line Curvature and Normals (https://www.mathworks.com/matlabcentral/fileexchange/32696-2d-line-curvature-and-normals), MATLAB Central File Exchange. Retrieved August 13, 2024.

Chebfun package: This contains the aaa algorithm and other necessary functions.
Ref: https://www.chebfun.org/ , T. A. Driscoll, N. Hale, and L. N. Trefethen, Chebfun Guide. Pafnuty Press, Oxford, 2014

Self-intersect: detect self-overlapping portions of a fire line.
Ref: Canos, A.J., 2024. Fast and robust self-intersections. MATLAB Central File Exchange. Retrieved August 1, 2024. URL:(https://www.mathworks.com/matlabcentral/fileexchange/13351-fast-and-robust-self-intersections).

interparc: interpolates a line to give desired resolution of points on the line.
Ref: D’Errico, J., 2024. interparc. MATLAB Central File Exchange. Retrieved August 1, 2024. 
URL: (https://www.mathworks.com/matlabcentral/fileexchange/34874-interparc).