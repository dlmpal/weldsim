#pragma once

#include "params.hpp"

namespace amrex::weldsim
{
    Real heat_source(const WeldParams &params,
                     const Geometry &gm,
                     int i, int j, int k, Real time);
}