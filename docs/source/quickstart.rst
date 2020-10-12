Quickstart guide
================

Perhaps the simplest thing one can try to do is read-in and hgrid file from a path on the disk.
There is more than one way of achieving this, but perhaps the simplest approach it to use the Hgrid class:

.. code::

    from pyschism import Hgrid
    hgrid = Hgrid.open('hgrid.gr3')

Now that you have instantiated the Hgrid class with a grid file, you may see a quick-plot of the mesh by calling the make_plot() method, and by specifying show=True:

.. code::

    hgrid.make_plot(show=True)



A more advanced workflow
^^^^^^^^^^^^^^^^^^^^^^^^

Suppose we want to setup a full tidal-only SCHISM run.
