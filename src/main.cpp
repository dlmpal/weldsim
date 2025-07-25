#include "single_level.hpp"
#include "level_build.hpp"
#include <AMReX_Amr.H>
#include <AMReX_ParmParse.H>
#include <filesystem>

using namespace amrex;
using namespace amrex::weldsim;

int main(int argc, char *argv[])
{
    Initialize(argc, argv);
    {
        ParmParse parser("amr");
        int max_level;
        parser.get("max_level", max_level);
        std::string plot_file;
        parser.get("plot_file", plot_file);
        bool plot_files_output;
        parser.get("plot_files_output", plot_files_output);

        // Single-level solver
        if (max_level == 0)
        {
            Print() << "Running single-level solver...\n";
            SingleLevelSolver solver;
            solver.init();
            solver.run();
        }
        // Multi-level solver
        else
        {
            Print() << "Running multi-level solver (max_level=" << max_level << ")\n";
            Amr amr(get_level_build());
            Real time_start = ThermalSolverLevel::params.time_start;
            Real time_stop = ThermalSolverLevel::params.time_stop;
            amr.init(time_start, time_stop);

            // Create the VisIt "movie" file
            std::ofstream visit_file;
            if (plot_files_output)
            {
                if (ParallelDescriptor::IOProcessor())
                {
                    std::filesystem::path path(plot_file);
                    if (not path.parent_path().empty())
                    {
                        std::filesystem::create_directories(path.parent_path());
                    }
                    visit_file.open(path.parent_path() / "plt.visit");
                }
            }

            int step = 0;
            while (amr.okToContinue() and amr.cumTime() < time_stop)
            {
                amr.coarseTimeStep(time_stop);

                if (plot_files_output)
                {
                    amr.writePlotFile();
                    if (ParallelDescriptor::IOProcessor())
                    {
                        visit_file << Concatenate(std::filesystem::path(plot_file).filename().string(), step + 1) + "/Header\n";
                    }
                }

                step++;
            }
        }
    }
    Finalize();
}