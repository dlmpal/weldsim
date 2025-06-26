#pragma once

#include <AMReX_BCRec.H>

namespace amrex::weldsim
{
    struct WeldParams
    {
        /// @brief Simulation time while the electrode is active
        Real weld_time;

        ///@brief Simulation time after welding has ended
        Real cooling_time;

        /// @brief Electrode velocity
        Real vel;

        /// @brief Electrode power
        Real Q;

        /// @brief Electrode position along y-axis
        Real center_y;

        /// @brief Electrode position along z-axis
        Real center_z;

        /// @brief Gaussian "radius"
        Real sigma;
    };

    struct HeatTransferParams
    {
        /// @brief Ambient temperature
        Real T_ambient;

        ///@brief Initial temperature
        Real T_init;

        /// @brief Convection heat transfer coefficient
        Real htc;

        /// @brief Emissivity
        Real eps;
    };

    enum class SolverType
    {
        Explicit = 0,
        Implicit = 1
    };

    struct ProbParams
    {
        /// @brief Simulation start time
        Real time_start;

        /// @brief Simulation stop time
        Real time_stop;

        /// @brief Whether to use explicit
        /// or implicit timestepping
        SolverType solver_type;

        /// @brief Level 0 timestep
        Real coarse_dt;

        /// @brief Points file name
        std::string points_file;

        /// @brief Coordinates of points for which to record
        /// temperature at each (level 0) timestep
        Vector<Array<Real, AMREX_SPACEDIM>> points;

        /// @brief Boundary condition record
        Vector<BCRec> bcr;

        /// @brief Parameters related to heat transfer
        HeatTransferParams ht_params;

        /// @brief Parameters related to the welding configuration
        WeldParams weld_params;
    };

    /// @brief Read problem parameters from the input file
    ProbParams read_params();
}