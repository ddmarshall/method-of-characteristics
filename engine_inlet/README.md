# Air-Inlet-Method-of-Characteristics-Analysis-Tool (AIMCAT)

This python-based tool contains modules used to model axisymmetric and two-dimensional supersonic air inlets. The primary method used 
is the method of characteristics which models the most complex regions of the supersonic flow.

# Requirements

The following packages are required: 
numpy, scipy, matplotlib, json, ... 

# How to Use

The intended workflow for this tool is as follows: 

1. configure input files 
    a. run file (.JSON) - contains parameters including freestream properties, gas constants, and mesh settings (see template)
    b. geometry file (.py) - class containing functions for y(x) and dy/dx(x) for centerbody and cowl surfaces (see template)
    c. plot profile (.JSON) - contains settings for generating figures after a solution has been run 
    d. saved solution (.pickle) - serialized python object from a previous run

2. run the main.py file by calling python main.py. You will be greeted a UI which can explored at your will 

3. start playing around! generate a mesh, save a solution, export, make a plot, etc. 


