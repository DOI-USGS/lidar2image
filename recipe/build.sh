#!/bin/bash

$PYTHON -m pip install . --no-deps -vv

mkdir -p lidar2image_processing/build
cd lidar2image_processing/build

cmake ${CMAKE_ARGS} \
    -GNinja \
    -DCMAKE_BUILD_TYPE=Release \
    -DISISROOT="${PREFIX}" \
    ..
ninja install
cd ../..

mkdir -p lidar_processing/build
cd lidar_processing/build

cmake ${CMAKE_ARGS} \
    -GNinja \
    -DISISROOT="${PREFIX}" \
    -DCMAKE_PREFIX_PATH="${PREFIX}" \
    -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
    ..

ninja install