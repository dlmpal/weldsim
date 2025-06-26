#include "advance_implicit.hpp"
#include "heat_source.hpp"
#include "material_properties.hpp"
#include "bc.hpp"
#include <AMReX_PhysBCFunct.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLABecLaplacian.H>

namespace amrex::weldsim
{
    //=============================================================================
    Array<Real, 2> setup_mlmg(MLMG &mlmg)
    {
        ParmParse parser("mlmg");

        int n_iter_max;
        parser.get("n_iter_max", n_iter_max);
        mlmg.setMaxIter(n_iter_max);

        int n_fmg_iter_max;
        parser.get("n_fmg_iter_max", n_fmg_iter_max);
        mlmg.setMaxFmgIter(n_fmg_iter_max);

        int verbose;
        parser.get("verbose", verbose);
        mlmg.setVerbose(verbose);

        int bottom_verbose = 0;
        parser.get("bottom_verbose", bottom_verbose);
        mlmg.setBottomVerbose(bottom_verbose);

        Real abs_tol;
        parser.get("abs_tol", abs_tol);

        Real rel_tol;
        parser.get("rel_tol", rel_tol);

        return {abs_tol, rel_tol};
    }
    //=============================================================================
    void advance_implicit(const ProbParams &params, const Geometry &gm,
                          const MultiFab &_T_old, MultiFab &_T_new, Real time, Real dt,
                          int lev, const MultiFab *_T_crse, const Amr *amr)
    {
        // Quick access to related BoxArray and DistributionMapping
        const auto &ba = _T_old.boxArray();
        const auto &dm = _T_old.DistributionMap();

        // Setup operator
        MLABecLaplacian linop({gm}, {ba}, {dm}, LPInfo{});
        linop.setMaxOrder(2);

        // Set physical/domain B.C.
        Array<LinOpBCType, AMREX_SPACEDIM> bc_lo;
        Array<LinOpBCType, AMREX_SPACEDIM> bc_hi;
        for (int idir = 0; idir < AMREX_SPACEDIM; idir++)
        {
            // Low
            if (lev == 0 or gm.Domain().smallEnd(idir) == amr->Geom(0).Domain().smallEnd(idir))
            {
                if (params.bcr[0].lo(idir) == BCType::ext_dir)
                {
                    bc_lo[idir] = LinOpBCType::Robin;
                }
                else
                {
                    bc_lo[idir] = LinOpBCType::Neumann;
                }
            }
            else
            {
                bc_lo[idir] = LinOpBCType::Dirichlet;
            }

            // High
            if (lev == 0 or gm.Domain().bigEnd(idir) == amr->Geom(0).Domain().bigEnd(idir))
            {
                if (params.bcr[0].hi(idir) == BCType::ext_dir)
                {
                    bc_hi[idir] = LinOpBCType::Robin;
                }
                else
                {
                    bc_hi[idir] = LinOpBCType::Neumann;
                }
            }
            else
            {
                bc_hi[idir] = LinOpBCType::Dirichlet;
            }
        }
        linop.setDomainBC({bc_lo}, {bc_hi});

        // Set coarse/fine interface B.C.
        // Note: This needs to be called before setLevelBC()
        if (lev > 0)
        {
            linop.setCoarseFineBC(_T_crse, amr->refRatio(lev));
        }

        // Fill the values for the Robin B.C. multifabs
        // For Neumann/adiabatic boundaries, no work is done
        MultiFab robin_a(ba, dm, 1, 1);
        MultiFab robin_b(ba, dm, 1, 1);
        MultiFab robin_f(ba, dm, 1, 1);
        bf_fill_robin(gm, params, _T_old, robin_a, robin_b, robin_f);
        linop.setLevelBC(0, nullptr, &robin_a, &robin_b, &robin_f);

        // Set scalar coefficients
        linop.setScalars(1.0, dt);

        // "a" coefficients (per cell), i.e. rho * cp
        MultiFab a_coeff(ba, dm, 1, 0);

        // "b" coefficients (per face), i.e. kappa
        Array<MultiFab, AMREX_SPACEDIM> b_coeff;
        for (int idim = 0; idim < AMREX_SPACEDIM; idim++)
        {
            b_coeff[idim].define(convert(ba, IntVect::TheDimensionVector(idim)), dm, 1, 0);
        }

        // Fill field coefficients (a, b)
        for (MFIter mfi(_T_old); mfi.isValid(); ++mfi)
        {
            const auto T = _T_old.array(mfi);

            auto a = a_coeff.array(mfi);
            const auto box = mfi.validbox();
            ParallelFor(box, [=](int i, int j, int k)
                        { a(i, j, k) = density(T(i, j, k)) * specific_heat(T(i, j, k)); });

            for (int idir = 0; idir < AMREX_SPACEDIM; idir++)
            {
                auto b = b_coeff[idir].array(mfi);
                const auto nodal_box = mfi.nodaltilebox(idir);
                ParallelFor(nodal_box, [=](int i, int j, int k)
                            { 
                              IntVect idx_right(AMREX_D_DECL(i,j,k));
                              if(idx_right[idir] == gm.ProbDomain().lo(idir))
                              {
                                b(i, j, k) = conductivity(T(idx_right));
                              }
                              else
                              {
                                IntVect idx_left = idx_right - IntVect::TheDimensionVector(idir);
                                b(i, j, k) = conductivity(0.5 * (T(idx_left) + T(idx_right)));
                              } });
            }
        }
        linop.setACoeffs(0, a_coeff);
        linop.setBCoeffs(0, GetArrOfConstPtrs<MultiFab>(b_coeff));

        // RHS
        MultiFab rhs(ba, dm, 1, 0);
        for (MFIter mfi(_T_old); mfi.isValid(); ++mfi)
        {
            const auto &box = mfi.validbox();
            const auto T = _T_old.array(mfi);
            auto rhs_arr = rhs.array(mfi);
            ParallelFor(box, [=](int i, int j, int k)
                        {   Real alpha = density(T(i, j, k)) * specific_heat(T(i, j, k));
                            rhs_arr(i, j, k) = T(i,j,k) * alpha + dt * heat_source(params.weld_params, gm, i, j, k, time); });
        }

        // Linear solver
        MLMG mlmg(linop);
        auto [abs_tol, rel_tol] = setup_mlmg(mlmg);
        mlmg.solve({&_T_new}, {&rhs}, abs_tol, rel_tol);
    }
}