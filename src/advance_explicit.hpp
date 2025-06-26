#pragma once

#include "params.hpp"
#include <AMReX_MultiFab.H>

namespace amrex::weldsim
{
    void advance_explicit(const ProbParams &params, const Geometry &gm,
                          MultiFab &T, MultiFab &T_new, Real time, Real dt);
}