from __future__ import annotations

from typing import Any


# --- Optional compiled modules ------------------------------------------------

try:
    from . import _calc as _calc  # compiled extension (.pyd/.so)
except Exception:
    _calc = None  # type: ignore[assignment]


try:
    from . import _vert_cpp as _vert_cpp  # compiled extension (.pyd/.so)
except Exception:
    _vert_cpp = None  # type: ignore[assignment]


# --- Re-export selected symbols ----------------------------------------------

if _vert_cpp is not None:
    calc_vert = _vert_cpp.calc_vert
    calc_flat_vert = _vert_cpp.calc_flat_vert
    verify_aw = _vert_cpp.verify_aw
    verify_prm = _vert_cpp.verify_prm
    verify_pow = _vert_cpp.verify_pow


# Optionally re-export anything from _calc you want at package level.
# If you don't want "import *" behavior, remove this block.
if _calc is not None:
    for _name in getattr(_calc, "__all__", dir(_calc)):
        if _name.startswith("_"):
            continue
        globals()[_name] = getattr(_calc, _name)


__all__ = [
    "_calc",
    "_vert_cpp",
    "calc_vert",
    "calc_flat_vert",
    "verify_aw",
    "verify_prm",
    "verify_pow",
]
