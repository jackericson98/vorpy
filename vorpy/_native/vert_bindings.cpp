#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "core/vert_core.hpp"


namespace py = pybind11;

static std::array<std::array<double, 3>, 4> as_locs_4x3(const py::array_t<double, py::array::c_style | py::array::forcecast>& a)
{
    if (a.ndim() != 2 || a.shape(0) != 4 || a.shape(1) != 3)
    {
        throw std::invalid_argument("locs must have shape (4,3)");
    }

    auto r = a.unchecked<2>();
    std::array<std::array<double, 3>, 4> out{};
    for (ssize_t i = 0; i < 4; ++i)
    {
        out[i] = {r(i, 0), r(i, 1), r(i, 2)};
    }
    return out;
}

static std::array<double, 4> as_rads_4(const py::array_t<double, py::array::c_style | py::array::forcecast>& a)
{
    if (a.ndim() != 1 || a.shape(0) != 4)
    {
        throw std::invalid_argument("rads must have shape (4,)");
    }

    auto r = a.unchecked<1>();
    return {r(0), r(1), r(2), r(3)};
}

static std::vector<std::array<double, 3>> as_locs_nx3(const py::array_t<double, py::array::c_style | py::array::forcecast>& a)
{
    if (a.ndim() != 2 || a.shape(1) != 3)
    {
        throw std::invalid_argument("test_locs must have shape (N,3)");
    }

    auto r = a.unchecked<2>();
    std::vector<std::array<double, 3>> out;
    out.reserve(static_cast<std::size_t>(a.shape(0)));

    for (ssize_t i = 0; i < a.shape(0); ++i)
    {
        out.push_back({r(i, 0), r(i, 1), r(i, 2)});
    }

    return out;
}

static std::vector<double> as_rads_n(const py::array_t<double, py::array::c_style | py::array::forcecast>& a)
{
    if (a.ndim() != 1)
    {
        throw std::invalid_argument("test_rads must have shape (N,)");
    }

    auto r = a.unchecked<1>();
    std::vector<double> out;
    out.reserve(static_cast<std::size_t>(a.shape(0)));

    for (ssize_t i = 0; i < a.shape(0); ++i)
    {
        out.push_back(r(i));
    }

    return out;
}

PYBIND11_MODULE(_vert_cpp, m)
{
    m.doc() = "C++ port of vorpy.src.calculations.vert";

    m.def(
        "calc_vert",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> locs,
           py::array_t<double, py::array::c_style | py::array::forcecast> rads) -> py::tuple
        {
            const auto L = as_locs_4x3(locs);
            const auto R = as_rads_4(rads);

            const auto res = vorpy::calc_vert(L, R);

            py::object loc1 = res.loc1 ? py::cast(*res.loc1) : py::none();
            py::object rad1 = res.rad1 ? py::cast(*res.rad1) : py::none();

            py::object loc2 = res.loc2 ? py::cast(*res.loc2) : py::none();
            py::object rad2 = res.rad2 ? py::cast(*res.rad2) : py::none();

            return py::make_tuple(loc1, rad1, loc2, rad2);
        },
        py::arg("locs"),
        py::arg("rads")
    );

    m.def(
        "calc_flat_vert",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> locs,
           py::array_t<double, py::array::c_style | py::array::forcecast> rads,
           bool power) -> py::tuple
        {
            const auto L = as_locs_4x3(locs);
            const auto R = as_rads_4(rads);

            const auto [p, rad] = vorpy::calc_flat_vert(L, R, power);

            if (!p.has_value() || !rad.has_value())
            {
                return py::make_tuple(py::none(), py::none());
            }

            return py::make_tuple(py::cast(*p), py::cast(*rad));
        },
        py::arg("locs"),
        py::arg("rads"),
        py::arg("power") = false
    );


    m.def(
        "verify_aw",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> loc,
           double rad,
           py::array_t<double, py::array::c_style | py::array::forcecast> test_locs,
           py::array_t<double, py::array::c_style | py::array::forcecast> test_rads)
        {
            if (loc.ndim() != 1 || loc.shape(0) != 3)
            {
                throw std::invalid_argument("loc must have shape (3,)");
            }

            auto lr = loc.unchecked<1>();
            const std::array<double, 3> L{lr(0), lr(1), lr(2)};

            const auto TL = as_locs_nx3(test_locs);
            const auto TR = as_rads_n(test_rads);

            if (TL.size() != TR.size())
            {
                throw std::invalid_argument("test_locs and test_rads must have matching N");
            }

            return vorpy::verify_aw(L, rad, TL, TR);
        },
        py::arg("loc"),
        py::arg("rad"),
        py::arg("test_locs"),
        py::arg("test_rads"));

    m.def(
        "verify_prm",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> loc,
           double rad,
           py::array_t<double, py::array::c_style | py::array::forcecast> test_locs)
        {
            if (loc.ndim() != 1 || loc.shape(0) != 3)
            {
                throw std::invalid_argument("loc must have shape (3,)");
            }

            auto lr = loc.unchecked<1>();
            const std::array<double, 3> L{lr(0), lr(1), lr(2)};

            const auto TL = as_locs_nx3(test_locs);

            return vorpy::verify_prm(L, rad, TL);
        },
        py::arg("loc"),
        py::arg("rad"),
        py::arg("test_locs"));

    m.def(
        "verify_pow",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> loc,
           double rad,
           py::array_t<double, py::array::c_style | py::array::forcecast> test_locs,
           py::array_t<double, py::array::c_style | py::array::forcecast> test_rads)
        {
            if (loc.ndim() != 1 || loc.shape(0) != 3)
            {
                throw std::invalid_argument("loc must have shape (3,)");
            }

            auto lr = loc.unchecked<1>();
            const std::array<double, 3> L{lr(0), lr(1), lr(2)};

            const auto TL = as_locs_nx3(test_locs);
            const auto TR = as_rads_n(test_rads);

            if (TL.size() != TR.size())
            {
                throw std::invalid_argument("test_locs and test_rads must have matching N");
            }

            return vorpy::verify_pow(L, rad, TL, TR);
        },
        py::arg("loc"),
        py::arg("rad"),
        py::arg("test_locs"),
        py::arg("test_rads"));
}
