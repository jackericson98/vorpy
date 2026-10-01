"""Render topology_diagnostic's actual pocket triangles with optional VTK.

Usage: python scripts/render_pocket_preview.py PATH_TO_DIAGNOSTIC_OUTPUT
The cutaway is explicitly illustrative; topology is measured on the intact mesh.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import vtk


def render(directory):
    meshes = json.loads((directory / "presentation_mesh.json").read_text())
    colors = {"outer": (0.24, 0.57, 0.68), "wall": (0.88, 0.56, 0.28),
              "floor": (1.0, 0.80, 0.38)}
    for cutaway in (False, True):
        renderer = vtk.vtkRenderer()
        renderer.SetBackground(1, 1, 1)
        for region, triangles in meshes.items():
            points, cells = vtk.vtkPoints(), vtk.vtkCellArray()
            for triangle in triangles:
                ids = [points.InsertNextPoint(*v) for v in triangle["vertices"]]
                cells.InsertNextCell(3, ids)
            poly = vtk.vtkPolyData()
            poly.SetPoints(points)
            poly.SetPolys(cells)
            mapper = vtk.vtkPolyDataMapper()
            if cutaway:
                plane = vtk.vtkPlane()
                plane.SetNormal(0, 1, 0)
                clip = vtk.vtkClipPolyData()
                clip.SetInputData(poly)
                clip.SetClipFunction(plane)
                clip.Update()
                mapper.SetInputData(clip.GetOutput())
            else:
                mapper.SetInputData(poly)
            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            actor.GetProperty().SetColor(*colors[region])
            actor.GetProperty().SetAmbient(0.42)
            actor.GetProperty().SetDiffuse(0.58)
            actor.GetProperty().SetSpecular(0.12)
            renderer.AddActor(actor)
        label = vtk.vtkTextActor()
        label.SetInput("DEEP POCKET | cutaway of measured boundary" if cutaway else
                       "DEEP POCKET | intact measured boundary")
        label.SetPosition(25, 850)
        label.GetTextProperty().SetFontSize(23)
        label.GetTextProperty().SetColor(0.12, 0.18, 0.23)
        renderer.AddViewProp(label)
        subtitle = vtk.vtkTextActor()
        subtitle.SetInput("Blue: outer body    Orange: pocket walls    Gold: closed floor")
        subtitle.SetPosition(25, 22)
        subtitle.GetTextProperty().SetFontSize(20)
        subtitle.GetTextProperty().SetColor(0.12, 0.18, 0.23)
        renderer.AddViewProp(subtitle)
        camera = renderer.GetActiveCamera()
        camera.SetFocalPoint(0, 0, 0)
        camera.SetViewUp(0, 1, 0) if not cutaway else camera.SetViewUp(0, 0, 1)
        camera.SetPosition(0, -50*np.sin(np.deg2rad(8)), 50*np.cos(np.deg2rad(8))) if not cutaway else camera.SetPosition(24, -48, 25)
        camera.ParallelProjectionOn()
        camera.SetParallelScale(13)
        renderer.ResetCameraClippingRange()
        window = vtk.vtkRenderWindow()
        window.SetOffScreenRendering(1)
        window.SetSize(1100, 900)
        window.SetMultiSamples(8)
        window.AddRenderer(renderer)
        window.Render()
        capture = vtk.vtkWindowToImageFilter()
        capture.SetInput(window)
        capture.Update()
        writer = vtk.vtkPNGWriter()
        writer.SetFileName(str(directory / ("pocket_cutaway.png" if cutaway else "pocket_intact.png")))
        writer.SetInputConnection(capture.GetOutputPort())
        writer.Write()
        window.Finalize()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    render(parser.parse_args().directory.resolve())
