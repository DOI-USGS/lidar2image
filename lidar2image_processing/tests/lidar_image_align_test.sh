#!/bin/sh


cameraFilename=data/AS15-M-1134-crop.lev1.cub
pvlFilename=${cameraFilename/.cub/.pvl}

./lidar_image_align -t -l data/RDR_hadley_LOLA.csv -i $cameraFilename

exit $STATUS
