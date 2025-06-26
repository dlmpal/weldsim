#include "heat_source.hpp"
#include <AMReX_Geometry.H>

namespace amrex::weldsim
{
    //=============================================================================
    Real heat_source(const WeldParams &p, const Geometry &gm, int i, int j, int k, Real time)
    {
        if (time > p.weld_time)
        {
            return 0.0;
        }

        // Distance of cell center from the electrode
        RealVect x_cell(AMREX_D_DECL(gm.CellCenter(i, 0),
                                     gm.CellCenter(j, 1),
                                     gm.CellCenter(k, 2)));
        RealVect x_src(AMREX_D_DECL(p.vel * time,
                                    p.center_y,
                                    p.center_z));
        Real r = (x_src - x_cell).vectorLength();

        // Compute the source value
        Real coeff = (p.Q / p.vel) / (2 * p.sigma * p.sigma * M_PI);
        return coeff * std::exp(-r * r / (2 * p.sigma * p.sigma));
    }
}
