# Selective-MumfordShah
Implementation of model 1 from "On Two Convex Variational Models and Their Iterative Solutions for Selective Segmentation of Images with Intensity Inhomogeneity"

See main.m file.

In order for successful segmentation, you will probably need to tweak the parameter "eta" in the optimisation part. Dependant on noise level, lambda may also need to be tweaked

Additionally, the threshold parameter "th" will need to be tweaked.



## Example input:
![Example output](input.PNG "Example input").

## Model output:
![Model output](exoutput.PNG "Example output").

## Segmentation by thresholding:
![Segmentation](contour.PNG "Segmentation").
