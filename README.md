# PySchism

## A Python interface for SCHISM model runs.

### Installation:

To install, make sure you have you have activated a Python>=3.6 virtual environment. Clone this repository, cd into it and install using pip:
``` bash
pip install .
```
It is highly recommended that you run this software using a Python virtual environment, and that you use Python >= 3.6.


### Usage examples:

##### Full Hgrid domain plot
``` python
from pyschism.mesh import Hgrid
hgrid = Hgrid.open('hgrid.gr3')
hgrid.make_plot(show=True)
```



