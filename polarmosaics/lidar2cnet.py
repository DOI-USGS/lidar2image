#! python

import argparse
import math
import sys

import numpy as np
import pandas as pd
import plio
from plio.io import isis_serial_number
from plio.io.io_controlnetwork import to_isis

import pvl
import spiceypy
import warnings


def parse_args():
    desc = "Process a lidar2image GCP file and an associated image file."
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('gcp_file', help="The lidar2image generate gcp file.")
    parser.add_argument('cube', help='The associated cube file.')
    parser.add_argument('output', help='The output control network file')
    parser.add_argument('--lat_sigma', default=10.0, type=float, help='Latitude Sigma Value')
    parser.add_argument('--lon_sigma', default=10.0, type=float, help='Longitude Sigma Value')
    parser.add_argument('--rad_sigma', default=15.0, type=float, help='Radius Sigma Value')
    parser.add_argument('-u', '--username', default='Unknown', help='The username assigned to the control network')
    parser.add_argument('-d', '--description', default='Unknown', help='The description of the control network')
    parser.add_argument('-p', '--pointid_prefix', default=None, help='The prefix to apply to all point ids')
    parser.add_argument('-s', '--pointid_suffix', default=None, help='The suffix to apply to all point ids')
    parser.add_argument('-c', '--chooser_name', default='lidar2image', help='The default point chooser name')
    parser.add_argument('-i', '--networkid', default='None', help='The ISIS network identifier')
    parser.add_argument('-k', '--keep-duplicates', dest='duplicates', default=False, help='If passed, keep duplicate values')
    parser.add_argument('--suppress-ref-warnings', dest='suppress', default=False, help="If true, suppress warnings related to reference measure")
    return parser.parse_args()


def compute_covariance(lat, lon, rad, latsigma=10., lonsigma=10., radsigma=15., semimajor_axis=None):
        lat = math.radians(lat)
        lon = math.radians(lon)
        if semimajor_axis is None:
            semimajor_axis = rad

        # SetSphericalSigmasDistance
        scaled_lat_sigma = np.divide(latsigma, semimajor_axis, dtype= np.float64)

        # This is specific to each lon.
        scaled_lon_sigma = np.divide(lonsigma , (math.cos(lat) * semimajor_axis), dtype= np.float64)

        # Calculate the covariance matrix in latitudinal coordinates
        # assuming no correlation
        cov = np.eye(3,3, dtype='float64')
        cov[0,0] = np.power(scaled_lat_sigma, 2)
        cov[1,1] = np.power(scaled_lon_sigma, 2)
        cov[2,2] = np.power(radsigma, 2)

        # Calculate the jacobian of the transformation from latitudinal to rectangular coordinates
        j = np.zeros((3,3), dtype='float64')
        cosphi = np.cos(lat)
        sinphi = np.sin(lat)
        coslambda = np.cos(lon)
        sinlambda = np.sin(lon)
        rcosphi = np.multiply(rad, cosphi)
        rsinphi = np.multiply(rad, sinphi)
        j[0,0] = np.multiply(-rsinphi, coslambda)
        j[0,1] = np.multiply(-rcosphi,  sinlambda)
        j[0,2] = np.multiply(cosphi, coslambda)
        j[1,0] = np.multiply(-rsinphi, sinlambda)
        j[1,1] = np.multiply(rcosphi, coslambda)
        j[1,2] = np.multiply(cosphi, sinlambda)
        j[2,0] = rcosphi
        j[2,1] = 0.
        j[2,2] = sinphi

        # Conjugate the latitudinal covariance matrix by the jacobian (error propagation formula)
        mat = j.dot(cov)
        mat = mat.dot(j.T)
        # Extract the upper triangle 
        rectcov = np.zeros((2,3), dtype='float64')
        rectcov[0,0] = mat[0,0]
        rectcov[0,1] = mat[0,1]
        rectcov[0,2] = mat[0,2]
        rectcov[1,0] = mat[1,1]
        rectcov[1,1] = mat[1,2]
        rectcov[1,2] = mat[2,2]

        return rectcov

def convert_gcs_to_body_fixed(lat, lon, rad):
    return spiceypy.latrec(rad,lon,lat)

def main(args):
    np.set_printoptions(precision=15)
    # Read the GCP file into a data frame
    gcp = pd.read_csv(args.gcp_file, names=['UTC Time', 'Apriori Longitude', 'Apriori Latitude',
                                            'Apriori Radius', 'A', 'B', 'C', 'Image Name',
                                            'Sample', 'Line'])
    # Grab the input image label
    label = pvl.load(args.cube)

    # Grab the target from the label
    target = plio.utils.utils.find_in_dict(label, 'TargetName')

    # Compute the ISIS3 serial number
    serial = isis_serial_number.generate_serial_number(args.cube)

    columns = ['id', 'pointType', 'serialnumber', 'measureType',
               'sample', 'line',
               'aprioriX', 'aprioriY', 'aprioriZ','aprioriCovar',
               'latitudeConstrained', 'longitudeConstrained', 'radiusConstrained']

    data = []
    for i, r in gcp.iterrows():
        lat, lon, rad = r[['Apriori Latitude', 'Apriori Longitude', 'Apriori Radius']]
        cov = compute_covariance(lat, lon, rad, latsigma=args.lat_sigma, lonsigma=args.lon_sigma, radsigma=args.rad_sigma)
        rect = convert_gcs_to_body_fixed(math.radians(r['Apriori Latitude']),
                                         math.radians(r['Apriori Longitude']),
                                         r['Apriori Radius'])
        data.append([i, 3, serial, 0, r.Sample, r.Line,
               rect[0], rect[1], rect[2], cov,
              1, 1, 1])

    df = pd.DataFrame(data, columns=columns)

    # Clean out the bad points
    df = df[(df['sample'] > 0) & (df['line'] > 0)]

    # Clean out the duplicates if not k
    if args.duplicates == False:
        df = df.drop_duplicates(['sample', 'line'], keep=False)

    # we want to add choosername=args.chooser_name back to this when it's added to plio
    with warnings.catch_warnings():
        if args.suppress != False:
            warnings.filterwarnings("ignore", message="Unable to identify referenceIndex")
        to_isis(df, args.output, mode='wb', version=2,
                targetname=target, description=args.description,
                username=args.username, pointid_prefix=args.pointid_prefix,
                pointid_suffix=args.pointid_suffix, 
                networkid=args.networkid)

if __name__ == '__main__':
    args = parse_args()

    sys.exit(main(args))
