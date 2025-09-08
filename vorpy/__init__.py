"""
VORPY - Voronoi analysis of molecular structures

Lightweight package init: importing `vorpy` should NOT import heavy, optional
dependencies (GUI, plotting, etc.) or reach into internal modules with side effects.
"""

from __future__ import annotations

from importlib import import_module


# Version (prefer local module, fall back to importlib.metadata if you move it later)
try:
    from .src.version import __version__
except Exception:  # pragma: no cover
    try:
        from importlib.metadata import version as _pkg_version  # py3.8+
        __version__ = _pkg_version("vorpy")
    except Exception:  # pragma: no cover
        __version__ = "0.0.0"


def run(*, file=None, load_files=None, settings=None, groups=None, exports=None):
    """
    Lazy entry point for the GUI / orchestrator.
    Imports only when called so that `import vorpy` stays lightweight.
    """
    Run = _lazy_import("vorpy.src.run", "Run")
    app = Run(
        file=file,
        load_files=load_files,
        settings=settings,
        groups=groups,
        exports=exports,
    )
    return app


def VorPyGUI(*args, **kwargs):
    """
    Backward-compatible accessor to the GUI class without importing GUI at module import time.
    Usage: gui_cls = vorpy.VorPyGUI(); app = gui_cls(...)
    """
    return _lazy_import("vorpy.src.GUI.vorpy_gui", "VorPyGUI")


def _lazy_import(module_path: str, attr: str):
    """Import `attr` from `module_path` on demand."""
    mod = import_module(module_path)
    return getattr(mod, attr)


__all__ = [
    "__version__",
    "run",
    "VorPyGUI",
]
