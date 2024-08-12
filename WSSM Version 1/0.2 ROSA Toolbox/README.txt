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