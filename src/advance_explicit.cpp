#include "advance_explicit.hpp"
#include "material_properties.hpp"
#include "heat_source.hpp"
#include "bc.hpp"
#include <AMReX_PhysBCFunct.H>

namespace amrex::weldsim
{
    //=============================================================================
    void advance_explicit(const ProbParams &params, const Geometry &gm,
                          MultiFab &_T, MultiFab &_T_new, Real time, Real dt)
    {
        // Update ghost cells
        GpuBndryFuncFab<BCFillExtDir> bf(BCFillExtDir{&params.ht_params});
        PhysBCFunct phys_bc(gm, params.bcr, bf);
        phys_bc(_T, 0, _T.nComp(), _T.nGrowVect(), time, 0);

#if (AMREX_SPACEDIM > 2)
        const auto dX = gm.CellSizeArray();
#else
        GpuArray<Real, 3> dX = {AMREX_D_DECL(gm.CellSize(0), gm.CellSize(1), 1.0)};
#endif
        const Real dV = dX[0] * dX[1] * dX[2];
        const GpuArray<Real, 3> dA = {dX[1] * dX[2],
                                      dX[0] * dX[2],
                                      dX[0] * dX[1]};

        for (MFIter mfi(_T_new); mfi.isValid(); ++mfi)
        {
            const auto box = mfi.validbox();
            const auto lo = box.loVect3d();
            const auto hi = box.hiVect3d();
            const auto T = _T.array(mfi);
            const auto T_new = _T_new.array(mfi);
            ParallelFor(box,
                        [=](int i, int j, int k)
                        {
                            Real flux_x = conductivity(0.5 * (T(i + 1, j, k) + T(i, j, k))) * (T(i + 1, j, k) - T(i, j, k)) / dX[0];
                            flux_x -= conductivity(0.5 * (T(i, j, k) + T(i - 1, j, k))) * (T(i, j, k) - T(i - 1, j, k)) / dX[0];

                            Real flux_y = conductivity(0.5 * (T(i, j + 1, k) + T(i, j, k))) * (T(i, j + 1, k) - T(i, j, k)) / dX[1];
                            flux_y -= conductivity(0.5 * (T(i, j, k) + T(i, j - 1, k))) * (T(i, j, k) - T(i, j - 1, k)) / dX[1];

#if (AMREX_SPACEDIM > 2)
                            Real flux_z = conductivity(0.5 * (T(i, j, k + 1) + T(i, j, k))) * (T(i, j, k + 1) - T(i, j, k)) / dX[2];
                            flux_z -= conductivity(0.5 * (T(i, j, k) + T(i, j, k - 1))) * (T(i, j, k) - T(i, j, k - 1)) / dX[2];
#endif

                            // Total diffusive flux
                            Real flux = flux_x * dA[0] + flux_y * dA[1];

#if (AMREX_SPACEDIM > 2)
                            flux += flux_z * dA[2];
#endif

                            // Electrode heat source
                            Real Q = heat_source(params.weld_params, gm, i, j, k, time);

                            Real alpha = density(T(i, j, k)) * specific_heat(T(i, j, k));
                            T_new(i, j, k) = T(i, j, k) + dt / alpha * (flux / dV + Q);
                        });
        }
    }
}