#!/usr/bin/env python
import pathlib
import setuptools
parent = pathlib.Path(__file__).parent.absolute()
conf = setuptools.config.read_configuration(parent / 'setup.cfg')
meta = conf['metadata']
setuptools.setup(
    name=meta['name'],
    version=meta['version'],
    author=meta['author'],
    author_email=meta['author_email'],
    description=meta['description'],
    long_description=meta['long_description'],
    long_description_content_type="text/markdown",
    url=meta['url'],
    packages=setuptools.find_packages(),
    python_requires='>=3.6',
    install_requires=[
        'matplotlib',
        'netCDF4',
        'pyproj',
        'shapely',
        'fiona',
        'f90nml',
        'ordered_set',
        'psutil',
        'paramiko',
        'scipy',
        'wget',
        'appdirs',
    ],
    entry_points={
        'console_scripts': [
            'plot_mesh = pyschism.cmd.plot_mesh:main',
            "schrun = pyschism.__main__:main",
            'tidal_run = pyschism.cmd.tidal_run:main',
        ]
    },
    tests_require=['nose'],
    test_suite='nose.collector',
)
