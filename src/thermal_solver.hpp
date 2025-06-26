#pragma once

#include "params.hpp"
#include <AMReX_AmrLevel.H>
#include <AMReX_Amr.H>

namespace amrex::weldsim
{
    class ThermalSolverLevel : public AmrLevel
    {
    public:
        ThermalSolverLevel() = default;

        ThermalSolverLevel(Amr &amr, int lev,
                           const Geometry &gm,
                           const BoxArray &ba,
                           const DistributionMapping &dm,
                           Real time);

        virtual ~ThermalSolverLevel();

        /// @brief Define data descriptors
        static void variable_setup();

        /// @brief Cleanup data descriptors
        static void variable_cleanup();

        /// @brief Initialize data at problem start-up
        virtual void initData() override;

        /// @brief Initialize data on the current level from its equivalent before regridding
        virtual void init(AmrLevel &old) override;

        /// @brief Initialize data on the current level from the coarse level
        virtual void init() override;

        /// @brief Compute initial timestep size for all levels
        /// @note Level 0 does this for all levels
        virtual void computeInitialDt(int finest_level, int sub_cycle,
                                      Vector<int> &n_cycle,
                                      const Vector<IntVect> &ref_ratio,
                                      Vector<Real> &dt_level,
                                      Real stop_time) override;

        /// @brief Compute timestep size for all levels after regridding
        virtual void computeNewDt(int finest_level, int sub_cycle,
                                  Vector<int> &n_cycle,
                                  const Vector<IntVect> &ref_ratio,
                                  Vector<Real> &dt_min,
                                  Vector<Real> &dt_level,
                                  Real stop_time, int post_regrid_flag) override;

        /// @brief Advance current level for one timestep
        virtual Real advance(Real time, Real dt,
                             int iteration, int ncycle) override;

        //! Do work after each time step
        virtual void post_timestep(int iteration) override;

        /// @brief Do work after regrid()
        virtual void post_regrid(int, int) override {}

        /// @brief Do work after init()
        virtual void post_init(Real) override {}

        /// @brief Error estimation for regridding.
        virtual void errorEst(TagBoxArray &tb,
                              int clearval, int tagval,
                              Real time, int n_error_buf = 0,
                              int ngrow = 0) override;

    private:
        ThermalSolverLevel &get_level(int level) const
        {
            return static_cast<ThermalSolverLevel &>(parent->getLevel(level));
        }

        static const int state_idx = 0;
        static const int n_comp = 1;
        static const int n_ghost = 1;

    public:
        static ProbParams params;
        static Real threshold;
        static std::ofstream points_file;
    };
}