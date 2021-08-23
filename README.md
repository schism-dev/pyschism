# PySCHISM
[![codecov](https://codecov.io/gh/schism-dev/pyschism/branch/main/graph/badge.svg?token=VE9PHEACBZ)](https://codecov.io/gh/schism-dev/pyschism)

## A Python interface for SCHISM model runs.

### Installation:

#### Pre-requisites
It is highly recommended that you run this software using a [Python virtual environment](https://gist.github.com/jreniel/c2dd4f2f68f9d8172355461b5337f236), and that you use Python>=3.6 (preferrably, using the latest available Python version is encouraged). You may use conda or venv to satisfy this dependency.
You should also have the cdunits library installed. In ubuntu systems this is achieved by:
```bash
apt-get install udunits-bin
```


#### Install option 1: pip
```bash
pip install pyschism
```

#### Install option 2: clone repo
To install, clone this repository, and navigate into it:
``` bash
git clone https://github.com/schism-dev/pyschism
cd pyschism
```
Then make sure to activate the target Python environment (this step is not necessary if you chose not to use a virtual environment).
After making sure your target environment is active, you can install the package using pip:

```bash
pip install .
```

#### If you are a developer
If you are a developer, it is recommended that you clone the repo.
After you add the `-e` flag to the pip install command in order to install in developer mode.

```bash
pip install -e .
```
---
### Usage examples:

#### Using the Library

##### Example 1: Full domain Hgrid plot:
``` python
from pyschism.mesh import Hgrid
hgrid = Hgrid.open('hgrid.gr3')
hgrid.make_plot(show=True)
```

##### Example 2: Write mesh to QGIS friendly format
```python
# NOTE: 2dm files can be read by QGIS > 3.0
from pyschism.mesh import Hgrid
hgrid = Hgrid.open('hgrid.gr3')
hgrid.write("/path/to/output/file.2dm", fmt='2dm')
```
---
Questions, comments and suggestions are welcome. Please follow the instructions on the `CONTRIBUTING.md` file for contributions. For bug reports and feature requests, please open an issue using the issue tracker.
Author contact: jrcalzada@wm.edu


---
    “Marconi is a good fellow. Let him continue. He is using seventeen of my patents.”
    Nikola Tesla
---