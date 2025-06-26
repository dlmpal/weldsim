#pragma once

#include "params.hpp"
#include <AMReX_Amr.H>

namespace amrex::weldsim
{
    void advance_implicit(const ProbParams &params, const Geometry &gm,
                          const MultiFab &T_old, MultiFab &T_new, Real time, Real dt,
                          int lev = 0, const MultiFab *T_crse = nullptr, const Amr *amr = nullptr);
}