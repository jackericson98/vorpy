import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QWidget
from PyQt5.QtCore import Qt
from OpenGL.GL import *
from OpenGL.GLU import *
from PyQt5.QtOpenGL import QGLWidget
import numpy as np


class GLWidget(QGLWidget):
    def __init__(self, parent=None):
        super(GLWidget, self).__init__(parent)
        self.setMinimumSize(400, 300)

    def initializeGL(self):
        glClearColor(0, 0, 0, 1)
        glEnable(GL_DEPTH_TEST)
        glEnable(GL_LIGHT0)
        glEnable(GL_LIGHTING)
        glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE)
        glEnable(GL_COLOR_MATERIAL)

    def resizeGL(self, width, height):
        glViewport(0, 0, width, height)
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        gluPerspective(45, width / float(height), 0.1, 200.0)
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()

    def paintGL(self):
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glLoadIdentity()
        glTranslatef(0.0, 0.0, -5.0)
        self.draw_sphere()

    def draw_sphere(self):
        glPushMatrix()
        glColor3f(1.0, 0, 0)  # Red color
        quadric = gluNewQuadric()
        gluSphere(quadric, 1.0, 32, 16)  # Draw sphere with radius 1.0
        glPopMatrix()


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('3D Viewer with PyQt5 and OpenGL')
        self.setCentralWidget(GLWidget())
        self.show()


if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = MainWindow()
    sys.exit(app.exec_())
