#pragma once

#include "params.hpp"
#include <AMReX_MultiFab.H>

namespace amrex::weldsim
{
    class SingleLevelSolver
    {
    public:
        SingleLevelSolver();
        void init();
        void run();

    private:
        /// @brief Problem parameters
        ProbParams params;

        /// @brief Verbosity
        int verbose;

        /// @brief Root plot file name
        std::string plot_file;

        /// @brief Plot interval
        int plot_int;

        static const int n_comp = 1;
        static const int n_ghost = 1;

        /// @brief Old timestep values
        MultiFab T_old;

        /// @brief New timestep values 
        MultiFab T_new;

        /// @brief Geometry
        Geometry gm;
    };
}