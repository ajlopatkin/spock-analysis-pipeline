# spock-analysis-pipeline

Automated Image Analysis Pipeline to Detect and Quantify Bacterial Colony Fluorescence

## Contents

- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Demo](#demo)
- [Results](#results)
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
- MATLAB Core
- Image Analysis Toolbox

Linux, Windows and Mac OSX are supported with a valid MATLAB installation.

# Installation Guide

To install, simply place copy the spock_code folder to your MATLAB folder structure, and add the functions to your path if desired.

# Demo

# Results

# Citation
