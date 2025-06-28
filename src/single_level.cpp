#include "single_level.hpp"
#include "advance_explicit.hpp"
#include "advance_implicit.hpp"
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <filesystem>

namespace amrex::weldsim
{
    //=============================================================================
    SingleLevelSolver::SingleLevelSolver()
    {
        // Read problem parameters
        params = read_params();

        // Setup domain
        Array<int, AMREX_SPACEDIM> n_cells;
        {
            ParmParse parser("amr");
            parser.get("n_cell", n_cells);
            parser.get("verbose", verbose);
            parser.get("plot_files_output", plot_files_output);
            parser.get("plot_file", plot_file);
            parser.get("plot_int", plot_int);
        }
        IntVect domain_low(AMREX_D_DECL(0, 0, 0));
        IntVect domain_high(AMREX_D_DECL(n_cells[0] - 1,
                                         n_cells[1] - 1,
                                         n_cells[2] - 1));
        Box domain(domain_low, domain_high);
        BoxArray ba(domain);
        DistributionMapping dm(ba);

        // Setup state multifabs
        T_old.define(ba, dm, n_comp, n_ghost);
        T_new.define(ba, dm, n_comp, n_ghost);

        // Setup geometry
        Array<Real, AMREX_SPACEDIM> x_low;
        Array<Real, AMREX_SPACEDIM> x_high;
        {
            ParmParse parser("geometry");
            parser.get("prob_lo", x_low);
            parser.get("prob_hi", x_high);
        }
        RealBox rb(x_low, x_high);
        IntArray periodicity = {AMREX_D_DECL(0, 0, 0)}; ///< Non-periodic
        gm.define(domain, rb, CoordSys::cartesian, periodicity);
    }
    //=============================================================================
    void SingleLevelSolver::init()
    {
        T_old.setVal(params.ht_params.T_init);
        T_new.setVal(params.ht_params.T_init);
    }
    //=============================================================================
    void SingleLevelSolver::run()
    {
        auto &weld_params = params.weld_params;

        std::ofstream visit_file, points_file;
        if (ParallelDescriptor::IOProcessor())
        {
            // Create the VisIt "movie" file
            if (plot_files_output)
            {
                std::filesystem::path path(plot_file);
                if (not path.parent_path().empty())
                {
                    std::filesystem::create_directories(std::filesystem::path(plot_file).parent_path());
                }
                visit_file.open(path.parent_path() / "plt.visit");
            }

            // Create the points file
            if (params.points.size() > 0)
            {
                std::filesystem::path path(params.points_file);
                if (not path.parent_path().empty())
                {
                    std::filesystem::create_directories(path.parent_path());
                }
                points_file.open(path);
                points_file << "Time,";
                for (std::size_t i = 0; i < params.points.size() - 1; i++)
                {
                    points_file << "Point" + std::to_string(i + 1) + ",";
                }
                points_file << "Point" + std::to_string(params.points.size()) + "\n";
            }
        }

        Real time = params.time_start;
        int step = 0;
        bool after_weld_end = false;
        while (time < params.time_stop)
        {
            // Increment time
            time += params.coarse_dt;
            step++;

            // Advance solution by one timestep
            if (params.solver_type == SolverType::Explicit)
            {
                advance_explicit(params, gm, T_old, T_new, time, params.coarse_dt);
            }
            else
            {
                advance_implicit(params, gm, T_old, T_new, time, params.coarse_dt);
            }
            MultiFab::Copy(T_old, T_new, 0, 0, n_comp, n_ghost);

            // Save timestep data to file
            if (plot_int > 0 and step % plot_int == 0)
            {
                // Print current timestep info
                if (verbose > 0)
                {
                    Print() << "Time: " << time << ", " << "maxT: " << T_new.max(0) << "\n";
                }

                // Create a plotfile of the current state
                if (plot_files_output)
                {
                    const std::string &filename = Concatenate(plot_file, step, 5);
                    WriteSingleLevelPlotfile(filename, T_new, {"T"}, gm, time, step);
                    if (ParallelDescriptor::IOProcessor())
                    {
                        visit_file << std::filesystem::path(filename).filename().string() << "/Header\n";
                    }
                }

                // Update points file
                if (ParallelDescriptor::IOProcessor() and not params.points.empty())
                {
                    points_file << time << ",";
                    for (std::size_t i = 0; i < params.points.size(); i++)
                    {
                        const auto idx = gm.CellIndex(params.points[i].data());
                        for (MFIter mfi(T_new); mfi.isValid(); ++mfi)
                        {
                            const auto box = mfi.validbox();
                            const auto arr = T_new.array(mfi);
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

            /// @todo Cleanup
            // Increase timestep size or plot interval after welding has stopped
            if (time > weld_params.weld_time + 10 and after_weld_end == false)
            {
                after_weld_end = true;
                if (params.solver_type == SolverType::Implicit)
                {
                    params.coarse_dt *= 10;
                }
                else
                {
                    plot_int = plot_int * 100;
                }
            }
        }
    }
}