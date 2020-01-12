import importlib
from pyschism.mesh import Mesh

__all__ = [
    "Mesh"
]

if importlib.util.find_spec("colored_traceback") is not None:
    import colored_traceback
    colored_traceback.add_hook(always=True)
