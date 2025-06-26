#pragma once

#include "thermal_solver.hpp"
#include <AMReX_LevelBld.H>

namespace amrex::weldsim
{
    class ThermalSolverLevelBld : public LevelBld
    {
    public:
        virtual void variableSetUp() override;
        virtual void variableCleanUp() override;
        virtual AmrLevel *operator()() override;
        virtual AmrLevel *operator()(Amr &amr, int level,
                                     const Geometry &gm,
                                     const BoxArray &ba,
                                     const DistributionMapping &dm,
                                     Real time) override;
    };

    LevelBld *get_level_build();
}