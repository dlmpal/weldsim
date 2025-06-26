#include "params.hpp"
#include <AMReX_ParmParse.H>

namespace amrex::weldsim
{
    //=============================================================================
    ProbParams read_params()
    {
        ProbParams params;

        // General parameters
        {
            ParmParse parser;

            // Solver type
            std::string type_str;
            parser.get("solver_type", type_str);
            if (type_str == "explicit")
            {
                params.solver_type = SolverType::Explicit;
            }
            else if (type_str == "implicit")
            {
                params.solver_type = SolverType::Implicit;
            }
            else
            {
                Warning("SolverType: " + type_str + " not supported. Defaulting to implicit");
                params.solver_type = SolverType::Implicit;
            }

            // Coarse dt
            parser.get("coarse_dt", params.coarse_dt);

            // Points file etc.
            parser.get("points_file", params.points_file);
            int n_points;
            parser.get("n_points", n_points);
            for (int i = 0; i < n_points; i++)
            {
                Array<Real, AMREX_SPACEDIM> point_coords;
                parser.get(("point_" + std::to_string(i + 1)).c_str(), point_coords);
                params.points.push_back(point_coords);
            }
        }

        // Geometry
        Array<Real, AMREX_SPACEDIM> x_low;
        Array<Real, AMREX_SPACEDIM> x_high;
        {
            ParmParse parser("geometry");
            parser.get("prob_lo", x_low);
            parser.get("prob_hi", x_high);
        }

        // Boundary conditions
        Array<int, AMREX_SPACEDIM> bc_lo;
        Array<int, AMREX_SPACEDIM> bc_hi;
        {
            ParmParse parser("bc");
            parser.get("lo", bc_lo);
            parser.get("hi", bc_hi);
        }
        params.bcr.emplace_back(bc_lo.data(), bc_hi.data());

        // Heat transfer parameters
        auto &ht_params = params.ht_params;
        {
            ParmParse parser("ht");
            parser.get("ambient_temp", ht_params.T_ambient);
            parser.get("initial_temp", ht_params.T_init);
            parser.get("heat_transfer_coefficient", ht_params.htc);
            parser.get("emissivity", ht_params.eps);
        }

        // Welding parameters
        WeldParams &weld_params = params.weld_params;
        {
            ParmParse parser("weld");
            parser.get("cooling_time", weld_params.cooling_time);
            parser.get("velocity", weld_params.vel);
            parser.get("power", weld_params.Q);
            parser.get("center_y", weld_params.center_y);
            parser.get("center_z", weld_params.center_z);
            parser.get("sigma", weld_params.sigma);
            weld_params.weld_time = (x_high[0] - x_low[0]) / weld_params.vel;
        }

        // Set simulation start and stop time
        params.time_start = 0.0;
        params.time_stop = weld_params.weld_time + weld_params.cooling_time;

        return params;
    }
}