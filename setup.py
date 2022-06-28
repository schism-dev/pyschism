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
except ValueError as e:
    if "time data '%cI' does not match format '%Y-%m-%dT%H:%M:%S%z'" in str(e):
        version = '0.0.0'
    else:
        raise

class BuildSchism(setuptools.Command):

    description = "build external SCHISM dependencies"

    user_options = [
        ('url=', None, 'Path for git clone of SCHISM source.'),
        ('branch=', None, 'Branch to use for install'),
        ('hydro=', None, 'Branch to use for install')
    ]

    def initialize_options(self):
        self.url = None
        self.branch = None
        self.hydro = None

    def finalize_options(self):
        self.url = 'https://github.com/schism-dev/schism' if self.url is None \
            else self.url
        self.branch = 'master' if self.branch is None else self.branch
        self.hydro = True if self.hydro is None else bool(self.hydro)

    def run(self):
        print(self.url, self.branch, self.hydro)
        # subprocess.check_call(["git", "clone", "submodules/jigsaw-python"])
        # # install jigsawpy
        # os.chdir(PARENT / 'submodules/jigsaw-python')
        # subprocess.check_call(["git", "checkout", "master"])
        # self.announce('INSTALLING JIGSAWPY', level=3)
        # subprocess.check_call(["python", "setup.py", "install"])
        # # install jigsaw
        # self.announce(
        #     'INSTALLING JIGSAW LIBRARY AND BINARIES FROM '
        #     'https://github.com/dengwirda/jigsaw-python', level=3)
        # os.chdir("external/jigsaw")
        # os.makedirs("build", exist_ok=True)
        # os.chdir("build")
        # gcc, cpp = self._check_gcc_version()
        # subprocess.check_call(
        #     ["cmake", "..",
        #      "-DCMAKE_BUILD_TYPE=Release",
        #      f"-DCMAKE_INSTALL_PREFIX={PYENV_PREFIX}",
        #      f"-DCMAKE_C_COMPILER={gcc}",
        #      f"-DCMAKE_CXX_COMPILER={cpp}",
        #      ])
        # subprocess.check_call(["make", f"-j{cpu_count()}", "install"])
        # libsaw_prefix = list(PYENV_PREFIX.glob("**/*jigsawpy*")).pop() / '_lib'
        # os.makedirs(libsaw_prefix, exist_ok=True)
        # envlib = PYENV_PREFIX / 'lib' / SYSLIB[platform.system()]
        # os.symlink(envlib, libsaw_prefix / envlib.name)
        # os.chdir(PARENT)
        # subprocess.check_call(
        #   ["git", "submodule", "deinit", "-f", "submodules/jigsaw-python"])

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
    packages=setuptools.find_packages(exclude=['tests', 'examples', 'docs', 'docker']),
    python_requires='>3.6',
    setup_requires=['wheel', 'setuptools_scm', 'setuptools>=41.2',
                    'netcdf-flattener>=1.2.0'],
    include_package_data=True,
    extras_require={'dev': ['coverage', 'flake8', 'nose']},
    cmdclass={
        "build_schism": BuildSchism
    },
    install_requires=[
        'pygeos',
        'geopandas',
        'netcdf-flattener>=1.2.0',
        'netCDF4',
        'f90nml',
        'psutil',
        'scipy',
        'wget',
        'cf-python',
        'metpy',
        'sqlalchemy',
        'pyugrid',
        'boto3',
        'rtree',
        'numba',
        'tqdm',
        'tqdm-logging-wrapper',
        'xmltodict',
        'cdsapi',
        'seawater',
        'xarray',
        'cfgrib',
        'zarr',
        'fsspec',
        'stormevents',
        'utm',
        # 'dask_geopandas @ git+git://github.com/geopandas/dask-geopandas.git@master',
    ],
    entry_points={'console_scripts': ['pyschism = pyschism.__main__:main']},
    tests_require=['nose'],
    test_suite='nose.collector',
)
