#!/usr/bin/env python
import pathlib
import setuptools  # type: ignore[import]
import subprocess
import sys

subprocess.check_call(
    [sys.executable, '-m', 'pip', 'install', '--upgrade', 'pip'])

subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'wheel'])

try:
    from dunamai import Version
except ImportError:
    subprocess.check_call(
            [sys.executable, '-m', 'pip', 'install', 'dunamai']
    )
    from dunamai import Version  # type: ignore[import]

try:
    version = Version.from_any_vcs().serialize()
except RuntimeError:
    version = '0.0.0'


parent = pathlib.Path(__file__).parent.absolute()
conf = setuptools.config.read_configuration(parent / 'setup.cfg')
meta = conf['metadata']
setuptools.setup(
    name=meta['name'],
    version=version,
    author=meta['author'],
    author_email=meta['author_email'],
    description=meta['description'],
    long_description=meta['long_description'],
    long_description_content_type="text/markdown",
    url=meta['url'],
    packages=setuptools.find_packages(),
    python_requires='>=3.6, <3.9',
    setup_requires=['wheel', 'setuptools_scm', 'setuptools>=41.2'],
    include_package_data=True,
    extras_require={'dev': ['coverage', 'flake8', 'nose']},
    install_requires=[
        'matplotlib',
        'netcdf-flattener>=1.2.0',
        'netCDF4',
        'pyproj',
        'shapely',
        'fiona',
        'f90nml',
        'psutil',
        'scipy',
        'wget',
        'appdirs',
        'cf-python',
        'sqlalchemy',
        'geopandas',
        'pyugrid',
        'pytz',
        'boto3',
        'rtree'
    ],
    entry_points={'console_scripts': ['pyschism = pyschism.__main__:main']},
    tests_require=['nose'],
    test_suite='nose.collector',
)
