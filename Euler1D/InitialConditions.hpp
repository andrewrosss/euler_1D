//
//  InitialConditions.hpp
//  Euler1D
//
//  Created by Andrew Ross on 2017-04-03.
//  Copyright © 2017 Andrew Ross. All rights reserved.
//

#ifndef InitialConditions_hpp
#define InitialConditions_hpp

#include <stdio.h>
#include <iostream>
#include <vector>
#include "OneDNode.hpp"
#include "FluxFunctions.h"
#include "OneDCell.hpp"

template <unsigned int N>
void SetInitialConditions(OneDCell<N> cell[], std::string IC, double &run_time,
                          int interior_cells, int nodes_per_cell, double Gamma)
{
    double RL, UL, PL, RR, UR, PR,
    if (IC == "E1")
    {
        RL = 2.281;
        UL = 164.83;
        PL = 201170.0;
        
        RR = 1.408;
        UR = 0.0;
        PR = 101100.0;
        
        //delta_t = CFL*min(delta_x/(fabs(UL) + sqrt(Gamma*PL/RL)), delta_x/(fabs(UR) + sqrt(Gamma*PR/RR)));
        if (strcmp(argv[3], "default") == 0)
        {
            run_time = 0.012;
        }
        x_shock = 2.0;
        RIEMANN(PL, UL, sound_sp(RL, PL, Gamma), Gamma, RL,
                PR, UR, sound_sp(RR, PR, Gamma), Gamma, RR,
                PATTRN,
                PM, UM,
                AML, AMR,
                USL, UHL, UTL,
                USR, UHR, UTR);
        std::cout << "USR = " << USR << "\n";
    }
    else if (IC == "E2")
    {
        RL = 1.045;
        UL = 200.0;
        PL = 300000.0;
        
        RR = 3.483;
        UR = 200.0;
        PR = 300000.0;
        
        //delta_t = CFL*min(delta_x/(fabs(UL) + sqrt(Gamma*PL/RL)), delta_x/(fabs(UR) + sqrt(Gamma*PR/RR)));
        if (strcmp(argv[3], "default") == 0)
        {
            run_time = 0.025;
        }
        x_shock = 2.0;
        RIEMANN(PL, UL, sound_sp(RL, PL, Gamma), Gamma, RL,
                PR, UR, sound_sp(RR, PR, Gamma), Gamma, RR,
                PATTRN,
                PM, UM,
                AML, AMR,
                USL, UHL, UTL,
                USR, UHR, UTR);
        std::cout << "UM = " << UM << "\n";
    }
    else if (IC == "E3")
    {
        RL = 1.598;
        UL = -383.64;
        PL = 91880.0;
        
        RR = 2.787;
        UR = -216.97;
        PR = 200000.0;
        
        //delta_t = CFL*min(delta_x/(fabs(UL) + sqrt(Gamma*PL/RL)), delta_x/(fabs(UR) + sqrt(Gamma*PR/RR)));
        if (strcmp(argv[3], "default") == 0)
        {
            run_time = 0.035;
        }
        x_shock = 5.0;
        RIEMANN(PL, UL, sound_sp(RL, PL, Gamma), Gamma, RL,
                PR, UR, sound_sp(RR, PR, Gamma), Gamma, RR,
                PATTRN,
                PM, UM,
                AML, AMR,
                USL, UHL, UTL,
                USR, UHR, UTR);
        std::cout << "UTR = " << UTR << "\n";
        std::cout << "UHR = " << UHR << "\n";
    }
    else if (IC == "E4")
    {
        RL = 4.696;
        UL = 0.0;
        PL = 404400.0;
        
        RR = 1.408;
        UR = 0.0;
        PR = 101100.0;
        
        //delta_t = CFL*min(delta_x/(fabs(UL) + sqrt(Gamma*PL/RL)), delta_x/(fabs(UR) + sqrt(Gamma*PR/RR)));
        if (strcmp(argv[3], "default") == 0)
        {
            run_time = 0.007;
        }
        x_shock = 5.0;
        RIEMANN(PL, UL, sound_sp(RL, PL, Gamma), Gamma, RL,
                PR, UR, sound_sp(RR, PR, Gamma), Gamma, RR,
                PATTRN,
                PM, UM,
                AML, AMR,
                USL, UHL, UTL,
                USR, UHR, UTR);
        std::cout << "UHL = " << UHL << "\n";
        std::cout << "UTL = " << UTL << "\n";
        std::cout << "AML = " << AML << "\n";
        std::cout << "PM = " << PM << "\n";
        std::cout << "AMR = " << AMR << "\n";
        std::cout << "UM = " << UM << "\n";
        std::cout << "USR = " << USR << "\n";
    }
}

#endif /* InitialConditions_hpp */
