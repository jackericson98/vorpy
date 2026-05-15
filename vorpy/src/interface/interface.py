import os
from dataclasses import dataclass, field

import numpy as np
import pandas as pd


@dataclass
class Interface:
    net: object
    balls1: set[int]
    balls2: set[int]
    name: str = "interface"

    balls: list[int] = field(default_factory=list)
    verts: pd.DataFrame | None = None
    edges: pd.DataFrame | None = None
    surfs: pd.DataFrame | None = None

    surface_area: float = 0.0
    mean_curvature: float | None = None
    gauss_curvature: float | None = None

    def __post_init__(self):
        self.balls1 = set(self.balls1)
        self.balls2 = set(self.balls2)

        self._validate()
        self.build()

    def _validate(self):
        if self.net is None:
            raise ValueError("Interface requires a solved network.")

        if getattr(self.net, "surfs", None) is None:
            raise ValueError("Interface requires net.surfs.")

        if getattr(self.net, "edges", None) is None:
            raise ValueError("Interface requires net.edges.")

        if getattr(self.net, "verts", None) is None:
            raise ValueError("Interface requires net.verts.")

        overlap = self.balls1.intersection(self.balls2)
        if overlap:
            raise ValueError(
                f"Interface selections overlap. "
                f"{len(overlap)} balls are present in both selections. "
                f"Example overlap indices: {sorted(overlap)[:10]}"
            )

    def build(self):
        self.surfs = self._get_interface_surfaces()
        self.balls = self._get_interface_balls()
        self.edges = self._get_interface_edges()
        self.verts = self._get_interface_verts()

        self._calculate_summary()

    def _get_interface_surfaces(self):
        surfs = self.net.surfs.copy()

        def is_interface_surface(balls):
            if balls is None or len(balls) != 2:
                return False

            b0, b1 = balls

            return (
                (b0 in self.balls1 and b1 in self.balls2)
                or
                (b0 in self.balls2 and b1 in self.balls1)
            )

        return surfs[surfs["balls"].apply(is_interface_surface)].copy()

    def _get_interface_balls(self):
        if self.surfs is None or len(self.surfs) == 0:
            return []

        balls = set()

        for surf_balls in self.surfs["balls"]:
            balls.update(surf_balls)

        return sorted(balls)

    def _get_interface_edges(self):
        if self.surfs is None or len(self.surfs) == 0:
            return self.net.edges.iloc[0:0].copy()

        iface_surf_indices = set(self.surfs.index)

        if "surfs" in self.net.edges.columns:
            def touches_iface_surface(surf_indices):
                if surf_indices is None:
                    return False

                return any(surf in iface_surf_indices for surf in surf_indices)

            return self.net.edges[self.net.edges["surfs"].apply(touches_iface_surface)].copy()

        iface_balls = set(self.balls)

        def touches_iface_ball(edge_balls):
            if edge_balls is None:
                return False

            return any(ball in iface_balls for ball in edge_balls)

        return self.net.edges[self.net.edges["balls"].apply(touches_iface_ball)].copy()

    def _get_interface_verts(self):
        if self.surfs is None or len(self.surfs) == 0:
            return self.net.verts.iloc[0:0].copy()

        iface_surf_indices = set(self.surfs.index)

        if "surfs" in self.net.verts.columns:
            def touches_iface_surface(surf_indices):
                if surf_indices is None:
                    return False

                return any(surf in iface_surf_indices for surf in surf_indices)

            return self.net.verts[self.net.verts["surfs"].apply(touches_iface_surface)].copy()

        iface_balls = set(self.balls)

        def touches_iface_ball(vert_balls):
            if vert_balls is None:
                return False

            return any(ball in iface_balls for ball in vert_balls)

        return self.net.verts[self.net.verts["balls"].apply(touches_iface_ball)].copy()

    def _calculate_summary(self):
        if self.surfs is None or len(self.surfs) == 0:
            self.surface_area = 0.0
            self.mean_curvature = None
            self.gauss_curvature = None
            return

        area_col = self._find_col(self.surfs, ["Surface Area", "surface_area", "sa"])
        mean_col = self._find_col(self.surfs, ["Mean Curvature", "mean_curv", "mean_curvature"])
        gauss_col = self._find_col(self.surfs, ["Gauss Curvature", "gauss_curv", "gauss_curvature"])

        if area_col is None:
            raise ValueError(
                "Could not calculate interface surface area because no surface-area column was found."
            )

        areas = self.surfs[area_col].astype(float).to_numpy()
        self.surface_area = float(np.sum(areas))

        if mean_col is not None and self.surface_area > 0:
            vals = self.surfs[mean_col].astype(float).to_numpy()
            self.mean_curvature = float(np.sum(vals * areas) / self.surface_area)

        if gauss_col is not None and self.surface_area > 0:
            vals = self.surfs[gauss_col].astype(float).to_numpy()
            self.gauss_curvature = float(np.sum(vals * areas) / self.surface_area)

    def export_surfs(self, directory=None):
        self._ensure_directory(directory)

        path = os.path.join(directory or os.getcwd(), f"{self.name}_surfs.csv")
        self.surfs.to_csv(path, index=True)

    def export_edges(self, directory=None):
        self._ensure_directory(directory)

        path = os.path.join(directory or os.getcwd(), f"{self.name}_edges.csv")
        self.edges.to_csv(path, index=True)

    def export_verts(self, directory=None):
        self._ensure_directory(directory)

        path = os.path.join(directory or os.getcwd(), f"{self.name}_verts.csv")
        self.verts.to_csv(path, index=True)

    def export_info(self, directory=None):
        self._ensure_directory(directory)

        path = os.path.join(directory or os.getcwd(), f"{self.name}_info.txt")

        with open(path, "w", encoding="utf-8") as info:
            info.write(f"{self.name}\n\n")
            info.write("Interface Summary\n\n")
            info.write(f"  Selection 1 Balls: {len(self.balls1)}\n")
            info.write(f"  Selection 2 Balls: {len(self.balls2)}\n")
            info.write(f"  Interface Balls: {len(self.balls)}\n")
            info.write(f"  Interface Surfaces: {len(self.surfs)}\n")
            info.write(f"  Interface Edges: {len(self.edges)}\n")
            info.write(f"  Interface Vertices: {len(self.verts)}\n\n")
            info.write(f"  Surface Area: {self.surface_area:.6f} Å²\n")
            info.write(f"  Mean Curvature: {self.mean_curvature}\n")
            info.write(f"  Gaussian Curvature: {self.gauss_curvature}\n")

    def export(self, directory=None, surfs=True, edges=True, verts=True, info=True):
        if surfs:
            self.export_surfs(directory=directory)

        if edges:
            self.export_edges(directory=directory)

        if verts:
            self.export_verts(directory=directory)

        if info:
            self.export_info(directory=directory)

    @staticmethod
    def _find_col(df, possible_names):
        for name in possible_names:
            if name in df.columns:
                return name

        return None

    @staticmethod
    def _ensure_directory(directory):
        if directory is not None:
            os.makedirs(directory, exist_ok=True)