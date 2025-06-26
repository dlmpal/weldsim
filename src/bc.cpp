#include "bc.hpp"
#include "material_properties.hpp"
#include <AMReX_PhysBCFunct.H>

namespace amrex::weldsim
{
    //=============================================================================
    // Linearized heat transfer coefficient
    // for radiation, i.e:
    // eps * sig * (T^4 - T0^4) =
    // eps * sig * (T-T0)*(T+T0)*(T^2+T0^2) =
    // h_rad * (T-T0)
    Real radiation_htc(Real T, Real T0, Real eps)
    {
        constexpr double sigma = 5.670374 * 1e-8;
        return eps * sigma * (T + T0) * (T * T + T0 * T0);
    };
    //=============================================================================
    BCFillExtDir::BCFillExtDir(const HeatTransferParams *params)
        : params_(params)
    {
    }
    //=============================================================================
    void BCFillExtDir::operator()(const amrex::IntVect &iv,
                                  amrex::Array4<amrex::Real> const &T,
                                  const int dcomp,
                                  const int numcomp,
                                  amrex::GeometryData const &gm,
                                  const amrex::Real time,
                                  const amrex::BCRec *bcr,
                                  const int bcomp,
                                  const int orig_comp) const
    {
        // Quick access
        Real T0 = params_->T_ambient;

        // Compute the temperature value at the ghost cell adjacent to the boundary
        auto compute_ghost_value = [T0](Real T, Real h, Real cond, Real dx)
        {
            dx = dx * 0.5;
            Real Twall = (h * T0 + cond / dx * T) / (h + cond / dx);
            return 2 * Twall - T;
        };

        const auto lo = gm.Domain().loVect();
        const auto hi = gm.Domain().hiVect();
        const auto dX = gm.CellSize();
        for (int idir = 0; idir < AMREX_SPACEDIM; idir++)
        {
            if (iv[idir] < lo[idir])
            {
                // First interior cell adjacent to the boundary
                IntVect loc = iv;
                loc[idir] = lo[idir];

                if (bcr->lo(idir) == BCType::ext_dir)
                {

                    Real h = params_->htc + radiation_htc(T(loc), T0, params_->eps);
                    Real cond = conductivity(T(loc));

                    T(iv) = compute_ghost_value(T(loc), h, cond, dX[idir]);
                }
                else if (bcr->lo(idir) == BCType::foextrap)
                {
                    T(iv) = -T(loc);
                }
            }
            else if (iv[idir] > hi[idir])
            {
                // First interior cell adjacent to the boundary
                IntVect loc = iv;
                loc[idir] = hi[idir];

                if (bcr->hi(idir) == BCType::ext_dir)
                {
                    Real h = params_->htc + radiation_htc(T(loc), T0, params_->eps);
                    Real cond = conductivity(T(loc));

                    T(iv) = compute_ghost_value(T(loc), h, cond, dX[idir]);
                }
                else if (bcr->hi(idir) == BCType::foextrap)
                {
                    T(iv) = -T(loc);
                }
            }
        }
    }
    //=============================================================================
    void bf_fill_robin(const Geometry &gm, const ProbParams &params,
                       const MultiFab &_T, MultiFab &robin_a,
                       MultiFab &robin_b, MultiFab &robin_f)
    {
        // Quick access
        const auto &ht_params = params.ht_params;
        Real T0 = ht_params.T_ambient;
        Real htc = ht_params.htc;
        Real eps = ht_params.eps;

        // Robin B.C. for convection + radiation:
        // h * T + k * dT/dx = h * T0 ->
        // a * T + b * dT/dx = f
        for (MFIter mfi(_T); mfi.isValid(); ++mfi)
        {
            for (int idir = 0; idir < AMREX_SPACEDIM; idir++)
            {
                const auto T = _T.array(mfi);
                const auto a = robin_a.array(mfi);
                const auto b = robin_b.array(mfi);
                const auto f = robin_f.array(mfi);

                // Low
                if (params.bcr[0].lo(idir) == BCType::ext_dir)
                {
                    const auto bndry_box = adjCellLo(mfi.validbox(), idir);
                    ParallelFor(bndry_box,
                                [=](int i, int j, int k)
                                {
                                    // Get the adjacent interior cell index
                                    IntVect loc(i, j, k);
                                    loc[idir] += 1;

                                    // Total heat transfer coefficient
                                    Real h = htc + radiation_htc(T(loc), T0, eps);

                                    // Conductivity
                                    Real cond = conductivity(T(loc));

                                    a(i, j, k) = h;
                                    b(i, j, k) = cond;
                                    f(i, j, k) = h * T0;
                                });
                }

                // High
                if (params.bcr[0].hi(idir) == BCType::ext_dir)
                {
                    const auto bndry_box = adjCellHi(mfi.validbox(), idir);
                    ParallelFor(bndry_box,
                                [=](int i, int j, int k)
                                {
                                    // Get the adjacent interior cell index
                                    IntVect loc(i, j, k);
                                    loc[idir] -= 1;

                                    // Total heat transfer coefficient
                                    Real h = htc + radiation_htc(T(loc), T0, eps);

                                    // Conductivity
                                    Real cond = conductivity(T(loc));

                                    a(i, j, k) = h;
                                    b(i, j, k) = cond;
                                    f(i, j, k) = h * T0;
                                });
                }
            }
        }
    }
}