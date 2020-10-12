#!/usr/bin/env python
import pathlib
# from dunamai import Version  # type: ignore[import]
import setuptools  # type: ignore[import]
parent = pathlib.Path(__file__).parent.absolute()
conf = setuptools.config.read_configuration(parent / 'setup.cfg')
meta = conf['metadata']
setuptools.setup(
    name=meta['name'],
    # version=Version.from_any_vcs().serialize(),
    author=meta['author'],
    author_email=meta['author_email'],
    description=meta['description'],
    long_description=meta['long_description'],
    long_description_content_type="text/markdown",
    url=meta['url'],
    packages=setuptools.find_packages(),
    python_requires='>=3.6',
    setup_requires=['setuptools_scm',
                    # 'dunamai',
                    'setuptools>=41.2'
                    ],
    include_package_data=True,
    extras_require={'dev': ['coverage',
                            # 'dunamai',
                            'flake8', 'nose']},
    install_requires=[
        # 'dunamai',
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
            # 'plot_mesh = pyschism.cmd.plot_mesh:main',
            # "schrun = pyschism.__main__:main",
            # 'tidal_run = pyschism.cmd.tidal_run:main',
            'plot = pyschism.__main__:plot'
        ]
    },
    tests_require=['nose'],
    test_suite='nose.collector',
)
