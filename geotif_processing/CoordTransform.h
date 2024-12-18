#ifndef _COORD_TRANSFORM_H_
#define _COORD_TRANSFORM_H_

#include <string>
#include <vector>
#include <gdal_priv.h>
#include <cpl_string.h>
#include <Eigen/Core>

namespace at
{
    std::vector<double> projectedToGeographic(std::vector<double> projCoord, GDALDataset* data);
    std::vector<double> geographicToProjected(std::vector<double> geoCoord, GDALDataset* data);
    std::vector<double> geographicToPixel(std::vector<double> geoCoord, GDALDataset* data);
    std::vector<double> pixelToProjected(std::vector<double> pixel, GDALDataset* data);
    std::vector<double> pixelToGeographic(std::vector<double> pixel, GDALDataset* data);
    std::vector<double> projectedToPixel(std::vector<double> projCoord, GDALDataset* data);
    std::vector<double> geographicToCartesian(std::vector<double> geoCoord, GDALDataset* data);   
    Eigen::Vector3f geographicToCartesian(Eigen::Vector3f lonLatRad);
    Eigen::Vector3f geographicToCartesian(Eigen::Vector3f lonLatRad, double radius);
}

#endif //  _COORD_TRANSFORM_H_































