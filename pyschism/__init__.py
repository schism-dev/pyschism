from pyschism.mesh import Mesh

__all__ = [
    "Mesh"
]

try:
    import colored_traceback
    colored_traceback.add_hook(always=True)
except ModuleNotFoundError:
    pass
