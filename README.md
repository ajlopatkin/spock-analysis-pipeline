# spock-analysis-pipeline

Automated Image Analysis Pipeline to Detect and Quantify Bacterial Colony Fluorescence

## Contents

- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Demo](#demo)
- [License](./LICENSE)
- [Citation](#citation)

# Overview

This is a custom image analysis pipeline for analysis of the Solid Media Portable Cell Killing (SPOCK) Assay. The image analysis is meant to identify cell colonies within a .TIFF image of a plate and collect various statistical properties of these colonies so that screening for adjuvants can be performed quickly and at broad scale.

Method:
Imported .tif files were first converted to greyscale and then binarized using an Otsu global image threshold to locate the placement of the nitrocellulose membrane. Since the location and angle of the membrane varied slightly between each image, rotational and spatial transformations were normalized to a rectangular grid aligned with the membrane. A Hough transform was used to detect all colonies within this region, and simple morphological operations were performed to eliminate small particles and image noise. To filter remaining technical artifacts such as slight membrane creasing, a 16 x 24 grid was placed within the determined rectangular region to identify coordinates with one, greater than one, or zero colonies present. Exactly one colony per grid coordinate was determined based on relative size and placement of all detected circles within an individual square (nearest to the center). Since the outline of killed colonies could be detected using fluorescence imaging, if no colony satisfied our defined criteria, the grid coordinate was designated as empty due to technical error. A colony mask was then generated within the nitrocellulose coordinates, and the average intensity of the complement of this mask (e.g., nitrocellulose membrane without the colonies) was used to adjust the background of each image. Finally, background-adjusted colony intensities were normalized by row, column, and plate average to account for edge-specific and image-specific fluorescence effects.

# Repo Contents

- [spock_code](./spock_code): `MATLAB` code for analysis and plotting.
- [docs](./docs): SPOCK paper.
- [data](./data): Sample data to run `MATLAB` code.

# System Requirements

## Hardware Requirements

The `spock-analysis-pipeline` code requires only a standard computer with enough RAM to support the manipulation of 5-8MB images. For minimal performance, this will be a computer with about 2 GB of RAM. Image manipulation can be taxing on the computer, so for optimal performance, we recommend a computer with the following specs:

RAM: 16+ GB  
CPU: 4+ cores, 3.3+ GHz/core

## Software Requirements

### OS and MATLAB Requirements

The code development version is tested on *Mac OSX* operating system, with MATLAB R2017b. The following toolboxes/modules are required for this code to execute:
- MATLAB 2017 or more recent
- Image Processing Toolbox

Earlier versions of MATLAB may have different versions of commands (such as `dir`) or may be missing commands altogether (such as `imbinarize`).

Linux, Windows and Mac OSX are supported with a valid MATLAB installation.

# Installation Guide

To install, simply place copy the spock_code folder to your MATLAB folder structure, and add the functions to your path if desired.

# Demo

The following data has been included for the purpose of reproducing results and demonstrating code functionality:
- Analysis: 2 raw plate images each from Replicates 2 and 3
- Plotting/Comparison: Processed and compiled results from all plates and all replicates

To run the code, follow the guidelines below:
## Analysis

To run the image analysis pipeline, copy the contents of either the [data/Replicate2](./data/Replicate2) or [data/Replicate3](./data/Replicate3) folders into the spock_code [Data](./spock_code/Data) folder. Then, run img_analysis.m; set im_debug=1 to see intermediate steps while the code is running. After the run is complete, you will see the results written into the spock_code [Matfiles](./spock_code/Matfiles) folder; this is a structure called `experiment` that further contains a 16x24 structure that captures all details of the plate, such as mean fluorescence, colony size, etc. Note that user interaction is required during the run to draw the outline of the plate; to do this, when the image is presented, click on the bottom-left of the plate, and then double-click on the top left; do the same for the bottom right and top right. This will place the grid coordinates according to your lines.

## Comparison/Plotting

To run the comparison code, copy the contents of [/data/Analysis](./data/Analysis/) into [/spock_code/Comparison](./spock_code/comparison) and run `plotting_comparison.m`. This will generate figures based on all of the data for all replicates contained in Matfiles1, Matfiles2 and Matfiles3. These figures replicate the publication figures from the Image Analysis portion.

# Citation
The citation for this project will be provided upon publication. 
