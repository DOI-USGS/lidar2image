# lidar2image

This library provides C++ and Python utilities in support of creating a ground control network using an LROC NAC image and coincident LOLA lidar shot data.

## Building the Program

This software manages dependencies using Anaconda.  If you do not have Anaconda, then you must [download](https://www.anaconda.com/download) it before proceeding.

In order to build a minimal set of programs necessary to create the ground control network, a user can follow these steps from a conda-enabled terminal:

1. Create and activate a conda environment.

        conda create -n lidar2image -f environment.yml
        conda activate lidar2image
1. Build the lidar_image_align binary

        cd lidar2image_processing
        mkdir build && cd build
        cmake ..
        make


