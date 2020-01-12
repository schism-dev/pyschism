# PySCHISM

## A Python interface for SCHISM model runs.

### Installation:

It is highly recommended that you run this software using a [Python virtual environment](https://gist.github.com/jreniel/c2dd4f2f68f9d8172355461b5337f236), and that you use Python==3.8.

To install clone this repository, and navigate into it:
``` bash
git clone https://github.com/schism-dev/pyschism && cd pyschism
```

### Usage examples:

#### Using the CLI
##### Example 1: Full domain hgrid plot from the terminal.
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
