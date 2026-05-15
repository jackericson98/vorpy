import pandas as pd

from vorpy.src.interface.interface import Interface


class FakeNet:
    def __init__(self):
        self.surfs = pd.DataFrame([
            {"balls": [0, 1], "Surface Area": 2.0, "Mean Curvature": 0.5, "Gauss Curvature": 0.1},
            {"balls": [1, 2], "Surface Area": 4.0, "Mean Curvature": 1.0, "Gauss Curvature": 0.2},
            {"balls": [2, 3], "Surface Area": 8.0, "Mean Curvature": 2.0, "Gauss Curvature": 0.3},
            {"balls": [0, 3], "Surface Area": 16.0, "Mean Curvature": 4.0, "Gauss Curvature": 0.4},
        ])

        self.edges = pd.DataFrame([
            {"balls": [0, 1, 4], "surfs": [0]},
            {"balls": [1, 2, 4], "surfs": [1]},
            {"balls": [2, 3, 4], "surfs": [2]},
            {"balls": [0, 3, 4], "surfs": [3]},
        ])

        self.verts = pd.DataFrame([
            {"balls": [0, 1, 4, 5], "surfs": [0]},
            {"balls": [1, 2, 4, 5], "surfs": [1]},
            {"balls": [2, 3, 4, 5], "surfs": [2]},
            {"balls": [0, 3, 4, 5], "surfs": [3]},
        ])


def test_basic_interface():
    net = FakeNet()

    iface = Interface(
        net=net,
        balls1={0, 1},
        balls2={2, 3},
        name="A_B_iface"
    )

    print("\n=== BASIC INTERFACE TEST ===")
    print("Interface name:", iface.name)
    print("Interface balls:", iface.balls)

    print("\nSurfaces:")
    print(iface.surfs)

    print("\nEdges:")
    print(iface.edges)

    print("\nVerts:")
    print(iface.verts)

    print("\nSummary:")
    print("Surface area:", iface.surface_area)
    print("Mean curvature:", iface.mean_curvature)
    print("Gauss curvature:", iface.gauss_curvature)

    assert len(iface.surfs) == 2
    assert set(iface.surfs.index) == {1, 3}
    assert iface.balls == [0, 1, 2, 3]
    assert len(iface.edges) == 2
    assert len(iface.verts) == 2
    assert iface.surface_area == 20.0

    expected_mean = ((4.0 * 1.0) + (16.0 * 4.0)) / 20.0
    expected_gauss = ((4.0 * 0.2) + (16.0 * 0.4)) / 20.0

    assert iface.mean_curvature == expected_mean
    assert iface.gauss_curvature == expected_gauss

    print("\nBASIC INTERFACE TEST PASSED")


def test_no_interface():
    net = FakeNet()

    iface = Interface(
        net=net,
        balls1={0},
        balls2={2},
        name="no_iface"
    )

    print("\n=== NO INTERFACE TEST ===")
    print("Interface surfaces:", len(iface.surfs))
    print("Interface edges:", len(iface.edges))
    print("Interface verts:", len(iface.verts))
    print("Surface area:", iface.surface_area)
    print("Mean curvature:", iface.mean_curvature)

    assert len(iface.surfs) == 0
    assert len(iface.edges) == 0
    assert len(iface.verts) == 0
    assert iface.surface_area == 0.0
    assert iface.mean_curvature is None
    assert iface.gauss_curvature is None

    print("\nNO INTERFACE TEST PASSED")


def test_overlap_error():
    net = FakeNet()

    print("\n=== OVERLAP ERROR TEST ===")

    try:
        Interface(
            net=net,
            balls1={0, 1},
            balls2={1, 2},
            name="bad_iface"
        )

    except ValueError as error:
        print("Caught expected error:")
        print(error)
        print("\nOVERLAP ERROR TEST PASSED")
        return

    raise AssertionError("Expected overlap ValueError was not raised.")


def test_export():
    net = FakeNet()

    iface = Interface(
        net=net,
        balls1={0, 1},
        balls2={2, 3},
        name="A_B_iface"
    )

    out_dir = "interface_test_outputs"
    iface.export(directory=out_dir)

    print("\n=== EXPORT TEST ===")
    print(f"Exported interface files to: {out_dir}")
    print("EXPORT TEST PASSED")


def main():
    test_basic_interface()
    test_no_interface()
    test_overlap_error()
    test_export()

    print("\nALL INTERFACE TESTS PASSED")


if __name__ == "__main__":
    main()
