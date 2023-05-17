README.txt

There are 3 modules in my code listing: 
	1. RayTracer.py
	2. BundleTest.py
	3. RMSComparison.py

1. RayTracer.py
	This module contains all the classes and functions necessary to model an optical ray tracer.

2. BundleTest.py
	This module provides the code for elements in 3 different scenarios:
		-single spherical refracting surface
        	-plano-convex singlet lens:
            		-with the plane surface facing the input
            		-with the convex surface facing the input
	To run the code for one of these scenarios, simply remove the """ at the part of the script which is labelled accordingly.
		Note: the code is originally set up to run the optical tracer for the single spherical refracting surface.
	The components are given the values as specified in the project script.
	To consider the effect of different elements, simply alter the values given when defining 'sphere1', 'plane1', or 'bun':
		The code will automatically define the output plane based on the input values.

3. RMSComparison.py
	This module is simply used to plot the RMS comparison graph for each lens scenario, as given in BundleTest.py.
	To run the code for one of these scenarios, simply remove the """ at the part of the script which is labelled accordingly.
		Note: the code is originally set up to produce an RMS comparison graph for the single spherical refracting surface.
	The components are given the values as specified in the project script.
	To change the lens being considered, simply alter the values given when defining 'sphere1', 'plane1', 'focal', or 'buni'.