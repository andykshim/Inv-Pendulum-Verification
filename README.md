# Inverted Pendulum Verification

### Authors:
Andrew Shim (akshim2)  
Richard Engel (reengle2)

Project for ECE 584 @ UIUC  
Prof. Sayan Mitra (mitras)




The project uses the tool Flow*.
Flow* is available at https://github.com/chenxin415/flowstar.git  
The *flowstar* folder contains the baseline model and examples from Flow*.  
The *flowstar_models* folder contains our single and double pendulum controller implemented using Flow*, as well as matlab files used to determine LQR parameters

## Directions
All source code can be compiled using **make**. Upon execution of the compiled code, the appropriate *output.plt* files will be generated in the *outputs* folder.  
Use **"gnuplot output.plt"** to generate a *.eps* file  
You can convert *.eps* to *.png* file via  **"convert output.eps output.png"**. 

The shell script *outputs.sh* will automatically convert to *.eps* and *.png*.


For certain plots, a frame-by-frame analysis may be useful.
Use our script via **"visualizer.py outputs/output.plt"**
