# PySCHISM

## A Python interface for SCHISM model runs.

### Installation:

It is highly recommended that you run this software using a Python virtual environment, and that you use Python >= 3.6. To install, clone this repository, navigate into it and install using pip:
``` bash
pip install .
```



### Usage examples:

#### Using the CLI
##### Example 1: Full domain Hgrid plot from the terminal.
``` bash
plot_mesh /path/to/hgrid.gr3 --plot-boundaries --plot-elements
```
![example_1_hgrid](https://raw.githubusercontent.com/schism-dev/pyschism/master/examples/example_1/hgrid.png)

#### Using the Library
##### Example 1: Full domain Hgrid plot:
``` python
from pyschism.mesh import Hgrid
hgrid = Hgrid.open('hgrid.gr3')
hgrid.make_plot(show=True)
```
