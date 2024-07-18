This quadfitter is an algorithm to reconstruct the time and position of light emitting particle events through 3d-multilateration. A theoretical overview is provided at: https://www.overleaf.com/project/64de99d86667cbc867a5bdaf.

To test the math of the quad-fitter on completely artificial data, run example.py after changing hpmts to an array of your choosing.

To run the Python version of the quadfitter, call "python3 quad_fitter.py [name of root file with events to be reconstructed]". For accurate graphs, you must manually change the variables "x", "y", and "z" to be the expected values of the reconstruction. The generated graphs will appear in the "graphs" folder. Example data for 3MeV electrons in different locations around the detector can be found in the "data" folder.

Finally, summary_graph.py provides a way to generate summary graphs like the one found in the theoretial overview. You will need to input your own data points.
