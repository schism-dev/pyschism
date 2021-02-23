# PySCHISM
![coverage](tests/coverage.svg)

## A Python interface for SCHISM model runs.

### Installation:

It is highly recommended that you run this software using a [Python virtual environment](https://gist.github.com/jreniel/c2dd4f2f68f9d8172355461b5337f236), and that you use Python>=3.6 (preferrably, using the latest available Python version is encouraged). You may use conda or venv to satisfy this dependency.

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
If you are a developer, you  may install in developer mode instead:
```bash
pip install -e .
```
---
### Usage examples:

#### Using the CLI
##### Example 1: Full domain hgrid plot from the terminal.
``` bash
plot_mesh /path/to/hgrid.gr3 --plot-boundaries --plot-elements
```
![example_1_hgrid](https://raw.githubusercontent.com/schism-dev/pyschism/master/examples/example_1/hgrid.png)

#### Using the Library
Hint: You can test the library functions from the command line (without having to write a .py file) by using `python -c` and wrapping the commands between a pair of quotes, for example:
```bash
python -c "
from pyschism.mesh import Hgrid
print(Hgrid.open('/path/to/hgrid.gr3'))
"
```
##### Example 1: Full domain Hgrid plot:
``` python
from pyschism.mesh import Hgrid
hgrid = Hgrid.open('hgrid.gr3')
hgrid.make_plot(show=True)
```
##### Example 2: Write boundaries to shapefile
```python
# open mesh as example above
from pyschism.mesh import Hgrid
hgrid = Hgrid.open('hgrid.gr3', crs="EPSG:3395")  # For shapefile output, coordinate reference system should be specified.
hgrid.write_boundaries("/path/to/output/dir", overwrite=True)
```

##### Example 3: Write mesh to QGIS friendly format
```python
# NOTE: 2dm files can be read by QGIS > 3.0
from pyschism.mesh import Hgrid
hgrid = Hgrid.open('hgrid.gr3')
hgrid.write("/path/to/output/file.2dm", fmt='2dm')
```
---
Questions, comments and suggestions are welcome. Please follow the instructions on the `CONTRIBUTING.md` file for contributions. For bug reports and feature requests, please open an issue using the issue tracker.
Author contact: jrcalzada@wm.edu

