"""Read-only analysis layers built on solved VorPy networks."""

from vorpy.src.analyze.apollonius import (
    AlphaComplex,
    ApolloniusComplex,
    ApolloniusSimplex,
    PrimalFeatureRef,
    UnsupportedFeature,
)
from vorpy.src.analyze.aw_interface_curvature import (
    AWInterfaceCurvatureResult,
    InterfaceEdgeRecord,
    InterfaceSurfaceRecord,
    analyze_aw_interface_curvature,
)

__all__ = [
    "AlphaComplex",
    "ApolloniusComplex",
    "ApolloniusSimplex",
    "PrimalFeatureRef",
    "UnsupportedFeature",
    "AWInterfaceCurvatureResult",
    "InterfaceEdgeRecord",
    "InterfaceSurfaceRecord",
    "analyze_aw_interface_curvature",
]
