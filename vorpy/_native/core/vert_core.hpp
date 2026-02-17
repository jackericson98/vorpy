#pragma once

#include <array>
#include <optional>
#include <utility>
#include <vector>

namespace vorpy
{
    struct VertResult
    {
        std::optional<std::array<double, 3>> loc1;
        std::optional<double>               rad1;

        std::optional<std::array<double, 3>> loc2;
        std::optional<double>               rad2;
    };


    // Core entry points --------------------------------------------------------

    // Equivalent of Python calc_vert(locs, rads)
    // locs: 4 points in 3D, rads: 4 radii
    // Returns up to two vertex solutions.
    VertResult calc_vert(const std::array<std::array<double, 3>, 4>& locs,
                         const std::array<double, 4>& rads);


    // Equivalent of Python calc_flat_vert(locs, rads, power=False)
    // Returns (loc, rad) if solvable, else (nullopt, nullopt).
    std::pair<std::optional<std::array<double, 3>>, std::optional<double>>
    calc_flat_vert(const std::array<std::array<double, 3>, 4>& locs,
                   const std::array<double, 4>& rads,
                   bool power);


    // Verification helpers -----------------------------------------------------

    bool verify_aw(const std::array<double, 3>& loc,
                   double rad,
                   const std::vector<std::array<double, 3>>& test_locs,
                   const std::vector<double>& test_rads);


    bool verify_prm(const std::array<double, 3>& loc,
                    double rad,
                    const std::vector<std::array<double, 3>>& test_locs);


    bool verify_pow(const std::array<double, 3>& loc,
                    double rad,
                    const std::vector<std::array<double, 3>>& test_locs,
                    const std::vector<double>& test_rads);

} // namespace vorpy
