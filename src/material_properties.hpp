#pragma once

namespace amrex::weldsim
{
    /// @brief Density
    double density(double T);

    /// @brief Thermal conductivity
    double conductivity(double T);

    /// @brief Specific heat
    double specific_heat(double T);
}