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
        'matplotlib',
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
        'cf-plot',
        'sqlalchemy',
        'pyugrid'
    ],
    entry_points={
        'console_scripts': [
            'pyschism = pyschism.__main__:main'
        ]
    },
    tests_require=['nose'],
    test_suite='nose.collector',
)
