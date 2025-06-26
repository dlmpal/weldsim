#pragma once

#include "params.hpp"
#include <AMReX_BCRec.H>
#include <AMReX_MultiFab.H>

namespace amrex::weldsim
{
    class BCFillExtDir
    {
    public:
        BCFillExtDir(const HeatTransferParams *params);

        void operator()(
            const amrex::IntVect &iv,
            amrex::Array4<amrex::Real> const &dest,
            const int dcomp,
            const int numcomp,
            amrex::GeometryData const &gm,
            const amrex::Real time,
            const amrex::BCRec *bcr,
            const int bcomp,
            const int orig_comp) const;

    private:
        const HeatTransferParams *params_;
    };

    class BCFillNull
    {
    public:
        BCFillNull() {};

        void operator()(
            const amrex::IntVect &iv,
            amrex::Array4<amrex::Real> const &dest,
            const int dcomp,
            const int numcomp,
            amrex::GeometryData const &gm,
            const amrex::Real time,
            const amrex::BCRec *bcr,
            const int bcomp,
            const int orig_comp) const {};
    };

    inline void bc_fill_null(Box const &bx, FArrayBox &data,
                             const int dcomp, const int numcomp,
                             Geometry const &geom, const Real time,
                             const Vector<BCRec> &bcr, const int bcomp,
                             const int scomp) {};

    void bf_fill_robin(const Geometry &gm, const ProbParams &params,
                       const MultiFab &T, MultiFab &robin_a,
                       MultiFab &robin_b, MultiFab &robin_f);
}