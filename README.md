# PySCHISM
[![codecov](https://codecov.io/gh/schism-dev/pyschism/branch/main/graph/badge.svg?token=VE9PHEACBZ)](https://codecov.io/gh/schism-dev/pyschism)
[![DOI](https://zenodo.org/badge/233075737.svg)](https://zenodo.org/badge/latestdoi/233075737)

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
## Online manual

https://schism-dev.github.io/schism/master/getting-started/pre-processing-with-pyschism/overview.html

## References

If you used this software as part of your work, please use the following citation format.

Jaime R Calzada, Linlin Cui, & Joseph Zhang. (2023). schism-dev/pyschism: v0.1.5 (v0.1.5). Zenodo. https://doi.org/10.5281/zenodo.7623122

---



Questions, comments and suggestions are welcome. Please follow the instructions on the `CONTRIBUTING.md` file for contributions. For bug reports and feature requests, please open an issue using the issue tracker.
Main author name: Jaime R Calzada
Author contact: jrcalzada@vims.edu


---
    “Marconi is a good fellow. Let him continue. He is using seventeen of my patents.”
    Nikola Tesla
---
