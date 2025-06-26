#include "material_properties.hpp"

namespace amrex::weldsim
{
    //=============================================================================
    double density(double T)
    {
        return 8050;
    }
    //=============================================================================
    double conductivity(double T)
    {
        double cond = 0;
        if (T >= 20 and T < 800)
        {
            cond = 54 - 3.33e-2 * T;
        }
        else if (T >= 800)
        {
            cond = 27.3;
        }
        return cond;
    }
    //=============================================================================
    double specific_heat(double T)
    {
        double sh = 0;
        if (T >= 20 and T < 600)
        {
            sh = 425 + (7.73e-1 - 1.69e-3 * T + 2.22e-6 * T * T) * T;
        }
        else if (T >= 600 and T < 735)
        {
            sh = 666 + 13002 / (738 - T);
        }
        else if (T >= 735 and T < 900)
        {
            sh = 545 + 17820 / (T - 731);
        }
        else if (T >= 900)
        {
            sh = 650;
        }
        return sh;
    }
}