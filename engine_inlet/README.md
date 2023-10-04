# Air-Inlet-Method-of-Characteristics-Analysis-Tool (AIMCAT)

This python-based tool contains modules used to model axisymmetric and 
two-dimensional supersonic air inlets. The primary method used is the method of 
characteristics which models the most complex regions of the supersonic flow.

# Requirements

The following packages are required: 
numpy, scipy, matplotlib, json, pandas

I have been exclusively working with Python 3.9 for the development of this code, 
so I can't promise compatibility with previous versions.

# How to Use

The intended workflow for this tool is as follows: 

1. configure input files 
    a. run file (.json) - contains parameters including freestream properties, 
        gas constants, and mesh settings (see Input folder)
    b. geometry file (.json) - class containing functions for y(x) and dy/dx(x) 
        for centerbody and cowl surfaces (see geometry folder)
    c. plot profile (.json) - contains settings for generating figures after a 
        solution has been run (see post_processing folder)

2. run the main class by importing aimcat, calling aimcat.main() with your .json
input files. See the example AIMCAT_Control_Script.py for an example on running the 
main class. 
