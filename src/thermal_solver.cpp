#include "thermal_solver.hpp"
#include "advance_implicit.hpp"
#include "bc.hpp"
#include <AMReX_ParmParse.H>
#include <filesystem>

namespace amrex::weldsim
{
    //=============================================================================
    ProbParams ThermalSolverLevel::params{};
    Real ThermalSolverLevel::threshold{};
    std::ofstream ThermalSolverLevel::points_file{};
    //=============================================================================
    ThermalSolverLevel::ThermalSolverLevel(Amr &amr, int lev,
                                           const Geometry &gm,
                                           const BoxArray &ba,
                                           const DistributionMapping &dm,
                                           Real time)
        : AmrLevel(amr, lev, gm, ba, dm, time)
    {
    }
    //=============================================================================
    ThermalSolverLevel::~ThermalSolverLevel()
    {
    }
    //=============================================================================
    void ThermalSolverLevel::variable_setup()
    {
        // Only read the parameters from file during the first call
        static bool init = false;
        if (init == false)
        {
            init = true;

            // Read the refinement temperature threshold
            ParmParse parser("amr");
            parser.get("threshold", threshold);

            // Read problem parameters
            params = read_params();

            // Initialize the points file
            if (ParallelDescriptor::IOProcessor())
            {
                // Create the points file
                if (params.points.size() > 0)
                {
                    std::filesystem::path path(params.points_file);
                    std::filesystem::create_directories(path.parent_path());
                    points_file.open(params.points_file);
                    points_file << "Time,";
                    for (std::size_t i = 0; i < params.points.size() - 1; i++)
                    {
                        points_file << "Point" + std::to_string(i + 1) + ",";
                    }
                    points_file << "Point" + std::to_string(params.points.size()) + "\n";
                }
            }
        }

        // State descriptor for state (i.e. temperature)
        desc_lst.addDescriptor(state_idx, IndexType::TheCellType(),
                               StateDescriptor::TimeCenter::Point,
                               n_ghost, n_comp, &cell_cons_interp);
        desc_lst.setComponent(state_idx, state_idx, {"T"}, params.bcr,
                              StateDescriptor::BndryFunc(bc_fill_null));
    }
    //=============================================================================
    void ThermalSolverLevel::variable_cleanup()
    {
        points_file.close();
        desc_lst.clear();
    }
    //=============================================================================
    void ThermalSolverLevel::initData()
    {
        auto &T_new = get_new_data(state_idx);
        T_new.setVal(params.ht_params.T_init);
    }
    //=============================================================================
    void ThermalSolverLevel::init(AmrLevel &old)
    {
        Real dt_new = parent->dtLevel(Level());
        Real time_cur = old.get_state_data(state_idx).curTime();
        Real time_prev = old.get_state_data(state_idx).prevTime();
        Real dt_old = time_cur - time_prev;
        setTimeLevel(time_cur, dt_old, dt_new);

        // Interpolate temperature values from the old to the new/current level
        auto &T_new = get_new_data(state_idx);
        FillPatch(old, T_new, n_ghost, time_cur, state_idx, 0, n_comp);
    }
    //=============================================================================
    void ThermalSolverLevel::init()
    {
        Real dt = parent->dtLevel(Level());
        Real time_cur = get_level(Level() - 1).state[state_idx].curTime();
        Real time_prev = get_level(Level() - 1).state[state_idx].prevTime();
        Real dt_old = (time_cur - time_prev) / static_cast<Real>(parent->MaxRefRatio(Level() - 1));
        setTimeLevel(time_cur, dt_old, dt);

        // Interpolate temperature values from the coarse to the fine (current) level
        auto &T_new = get_new_data(state_idx);
        FillCoarsePatch(T_new, 0, time_cur, state_idx, 0, n_comp, n_ghost);
    }
    //=============================================================================
    void ThermalSolverLevel::computeInitialDt(int finest_level, int sub_cycle,
                                              Vector<int> &n_cycle,
                                              const Vector<IntVect> &ref_ratio,
                                              Vector<Real> &dt_level,
                                              Real stop_time)
    {
        // The coarsest level does this for all levels
        if (Level() > 0)
        {
            return;
        }

        // Number of fine level steps for each level 0 step:
        // For the levels i and i+1, with a refinement ratio of r,
        // the following is true: dt_(i+1) = 1/r * dt_i
        // If the refinement ratio is constant for all levels,
        // this recursively leads to the following:
        // dt_i = 1/r^i * dt_0
        // n_steps stores the value of (1/r^i) for each level
        Vector<int> n_steps(n_cycle.size());
        std::partial_sum(n_cycle.begin(),
                         n_cycle.end(),
                         n_steps.begin(),
                         std::multiplies<int>());

        // Compute the value of the level 0 timestep
        // This is computed from CFL considerations, either on
        // level 0, or finer levels
        Real dt_0 = params.coarse_dt;
        if (state[state_idx].curTime() > params.weld_params.weld_time + 10)
        {
            dt_0 *= 10;
        }

        // Make sure that advancing by dt_0 does not exceed
        // the simulation's stop time
        if (stop_time > 0)
        {
            const Real eps = 0.001 * dt_0;
            const Real cur_time = state[state_idx].curTime();
            if ((cur_time + dt_0) > (stop_time - eps))
            {
                dt_0 = stop_time - cur_time;
            }
        }

        // Compute the timestep values for finer levels
        for (int ilev = 0; ilev <= finest_level; ++ilev)
        {
            dt_level[ilev] = dt_0 / Real(n_steps[ilev]);
        }
    }
    //=============================================================================
    void ThermalSolverLevel::computeNewDt(int finest_level, int sub_cycle,
                                          Vector<int> &n_cycle,
                                          const Vector<IntVect> &ref_ratio,
                                          Vector<Real> &dt_min,
                                          Vector<Real> &dt_level,
                                          Real stop_time, int post_regrid_flag)
    {
        // Nothing new to be done here...
        computeInitialDt(finest_level, sub_cycle, n_cycle,
                         ref_ratio, dt_level, stop_time);
    }
    //=============================================================================
    void ThermalSolverLevel::post_timestep(int iteration)
    {
        // Correct current level values by a weighted average
        // from values on the immediately finer level (average down)
        if (Level() < parent->finestLevel())
        {
            const auto &fine_level = get_level(Level() + 1);
            const auto &T_fine = fine_level.get_new_data(state_idx);
            auto &T_crse = get_new_data(state_idx);
            average_down(T_fine, T_crse, fine_level.Geom(), Geom(),
                         0, n_comp, parent->refRatio(Level()));
        }

        // Update points file
        if (Level() == 0)
        {
            // Update points file
            if (ParallelDescriptor::IOProcessor() and not params.points.empty())
            {
                points_file << state[state_idx].curTime() << ",";
                for (std::size_t i = 0; i < params.points.size(); i++)
                {
                    const auto idx = Geom().CellIndex(params.points[i].data());
                    for (MFIter mfi(get_new_data(state_idx)); mfi.isValid(); ++mfi)
                    {
                        const auto box = mfi.validbox();
                        const auto arr = get_new_data(state_idx).array(mfi);
                        if (arr.contains(idx))
                        {
                            points_file << arr(idx);
                            if (i < params.points.size() - 1)
                            {
                                points_file << ",";
                            }
                        }
                    }
                }
                points_file << "\n";
            }
        }
    }
    //=============================================================================
    Real ThermalSolverLevel::advance(Real time, Real dt, int iteration, int ncycle)
    {
        state[state_idx].allocOldData();
        state[state_idx].swapTimeLevels(dt);

        const auto &T_old = get_old_data(state_idx);
        auto &T_new = get_new_data(state_idx);
        FillPatch(*this, T_new, n_ghost, time, state_idx, 0, n_comp);

        if (Level() == 0)
        {
            advance_implicit(params, Geom(),
                             T_old, T_new, time, dt,
                             Level());
        }
        else
        {
            // const auto &T_crse = get_level(Level() - 1).get_new_data(state_idx);
            auto &crse_level = get_level(Level() - 1);
            MultiFab T_crse(crse_level.boxArray(), crse_level.DistributionMap(), n_comp, n_ghost);
            FillPatch(crse_level, T_crse, n_ghost, state[state_idx].curTime(), state_idx, 0, n_comp);
            advance_implicit(params, Geom(),
                             T_old, T_new, time, dt,
                             Level(), &T_crse, parent);
        }

        return dt;
    }
    //=============================================================================
    void ThermalSolverLevel::errorEst(TagBoxArray &tb, int clearval, int tagval,
                                      Real time, int n_error_buf, int ngrow)
    {
        const auto &T = get_new_data(state_idx).const_arrays();
        const auto &a = tb.arrays();
        Real T_norm = get_new_data(state_idx).norm2(0);
        ParallelFor(tb, [=] AMREX_GPU_DEVICE(int bi, int i, int j, int k)
                    {
                        if (T[bi](i, j, k) > threshold)
                        {
                            a[bi](i,j,k) =  TagBox::SET;
                        } });
    }
}