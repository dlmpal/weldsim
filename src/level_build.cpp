#include "level_build.hpp"

namespace amrex::weldsim
{
    //=============================================================================
    void ThermalSolverLevelBld::variableSetUp()
    {
        ThermalSolverLevel::variable_setup();
    }
    //=============================================================================
    void ThermalSolverLevelBld::variableCleanUp()
    {
        ThermalSolverLevel::variable_cleanup();
    }
    //=============================================================================
    AmrLevel *ThermalSolverLevelBld::operator()()
    {
        return new ThermalSolverLevel;
    }
    //=============================================================================
    AmrLevel *ThermalSolverLevelBld::operator()(Amr &amr, int level,
                                                const Geometry &gm,
                                                const BoxArray &ba,
                                                const DistributionMapping &dm,
                                                Real time)
    {
        return new ThermalSolverLevel(amr, level, gm, ba, dm, time);
    }
    //=============================================================================
    // Global instance
    ThermalSolverLevelBld lvlbld;
    LevelBld *get_level_build()
    {
        return &lvlbld;
    }
}