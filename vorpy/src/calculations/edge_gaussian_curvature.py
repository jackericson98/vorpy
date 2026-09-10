"""Network-level edge contribution to integrated Gaussian curvature."""

from vorpy.src.calculations.edge_curvature import calculate_aw_network_edge_curvatures


def calculate_aw_network_edge_gaussian_curvatures(net, quadrature_order=32, tolerance=1e-6):
    """Return edge-G from the fused production edge-curvature traversal."""
    _, gaussian_values = calculate_aw_network_edge_curvatures(
        net,
        quadrature_order=quadrature_order,
        tolerance=tolerance,
    )
    return gaussian_values
