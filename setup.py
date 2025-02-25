from setuptools import setup, find_packages

setup(
        name="lidar2image",
        version='1.0.0',
        packages=find_packages(),
        author = "USGS Astrogeology",
        description = "Utilities for creating ground control networks from LIDAR data",
        license = "CC0-1.0",
        url = "https://github.com/DOI-USGS/lidar2image",
        include_package_data=True,
        install_requires=[
            'numpy',
            'pandas',
            'plio',
            'pvl',
            'spiceypy'],
        entry_points={
            'console_scripts': [ 'lidar2cnet = polarmosaics.lidar2cnet:main']
        }

)
