#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <stdexcept>
#include <cstddef>

#include "core/calcs.hpp"


namespace py = pybind11;

static double calc_dist(py::array_t<double, py::array::c_style | py::array::forcecast> a,
                        py::array_t<double, py::array::c_style | py::array::forcecast> b)
{
    if (a.ndim() != 1 || b.ndim() != 1) throw std::runtime_error("calc_dist expects 1D arrays");
    if (a.shape(0) != b.shape(0)) throw std::runtime_error("calc_dist expects same length");

    const auto n = static_cast<std::size_t>(a.shape(0));
    py::gil_scoped_release release;
    return calc_dist_ptr(a.data(), b.data(), n);
}

static double calc_tri_3pts(
    py::array_t<double, py::array::c_style | py::array::forcecast> p0,
    py::array_t<double, py::array::c_style | py::array::forcecast> p1,
    py::array_t<double, py::array::c_style | py::array::forcecast> p2)
{
    if (p0.ndim()!=1 || p1.ndim()!=1 || p2.ndim()!=1) throw std::runtime_error("calc_tri: expected 1D vectors");
    if (p0.shape(0)!=3 || p1.shape(0)!=3 || p2.shape(0)!=3) throw std::runtime_error("calc_tri: expected length-3 vectors");

    py::gil_scoped_release release;
    return calc_tri_pts3(p0.data(), p1.data(), p2.data());
}


static double calc_tri_tri(
    py::array_t<double, py::array::c_style | py::array::forcecast> tri)
{
    if (tri.ndim()!=2 || tri.shape(0)!=3 || tri.shape(1)!=3) throw std::runtime_error("calc_tri: expected shape (3,3)");

    const double* t = tri.data();
    const double* p0 = t + 0;
    const double* p1 = t + 3;
    const double* p2 = t + 6;

    py::gil_scoped_release release;
    return calc_tri_pts3(p0, p1, p2);
}


PYBIND11_MODULE(_calc, m)
{
    m.doc() = "VorPy native kernels (pybind11)";

    m.def("calc_dist", &calc_dist, "Euclidean distance between two 1D vectors");
    m.def("calc_tri",
      py::overload_cast<
          py::array_t<double, py::array::c_style | py::array::forcecast>,
          py::array_t<double, py::array::c_style | py::array::forcecast>,
          py::array_t<double, py::array::c_style | py::array::forcecast>
      >(&calc_tri_3pts),
      "Triangle area from 3 points");

    m.def("calc_tri",
      py::overload_cast<
          py::array_t<double, py::array::c_style | py::array::forcecast>
      >(&calc_tri_tri),
      "Triangle area from a (3,3) array");

}

