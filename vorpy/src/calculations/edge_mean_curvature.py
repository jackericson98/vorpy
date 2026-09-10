"""Network-level edge contributions to integrated mean curvature."""

from vorpy.src.calculations.edge_curvature import calculate_aw_network_edge_curvatures


def calculate_aw_network_edge_mean_curvatures(net, quadrature_order=32, tolerance=1e-6):
    """Return edge-M from the fused production edge-curvature traversal."""
    mean_values, _ = calculate_aw_network_edge_curvatures(
        net,
        quadrature_order=quadrature_order,
        tolerance=tolerance,
    )
    return mean_values
