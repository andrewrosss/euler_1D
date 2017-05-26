//
//  FluxFunctions.h
//  Euler1D
//
//  Created by Andrew Ross on 2017-03-28.
//  Copyright © 2017 Andrew Ross. All rights reserved.
//

#ifndef FluxFunctions_h
#define FluxFunctions_h

#include "OneDNode.hpp"
#include "OneDCell.hpp"
#include "OneDGrid.hpp"
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <cassert>
#include <fstream>
#include <cstdlib>


/////////////////////////////////////////////////////
//
// zeroth-order interpolation to reconstruct primitive
// vars at cell interfaces
//
/////////////////////////////////////////////////////

double phi(double u_L,
           double u,
           double u_R,
           double u_interp)
{
    if (u_interp - u > 0)
    {
        double u_max = fmax(fmax(u_L, u), u_R);
        return fmin((u_max - u)/(u_interp - u), 1.0);
    }
    else if (u_interp - u < 0)
    {
        double u_min = fmin(fmin(u_L, u), u_R);
        return fmin((u_min - u)/(u_interp - u), 1.0);
    }
    else
    {
        return 1.0;
    }
}

// ADD_28


template <unsigned int C, unsigned int N>
void cell_bd_values_FV_0(OneDGrid<C, N>& grid)
{
    
    double density_left, velocity_left, pressure_left;
    double density_right, velocity_right, pressure_right;
    
    for (int i = 0; i < grid.NumberOfCells(); i++)
    {
        // store the prim vars interpolated to the left
        // boundary (0th-order)
        density_left = grid.Cell(i).Node(0).GetDensity();
        velocity_left = grid.Cell(i).Node(0).GetVelocityX();
        pressure_left = grid.Cell(i).Node(0).GetPressure();
        
        // store the prim vars interpolated to the right
        // boundary (0th-order)
        density_right = grid.Cell(i).Node(0).GetDensity();
        velocity_right = grid.Cell(i).Node(0).GetVelocityX();
        pressure_right = grid.Cell(i).Node(0).GetPressure();
        
        // store the interpolated values in the appropriate
        // interface nodes
        grid.Cell(i).LeftInterface().SetPrim(density_left,
                                             velocity_left,
                                             pressure_left);
        grid.Cell(i).RightInterface().SetPrim(density_right,
                                              velocity_right,
                                              pressure_right);
    }
    
    // load the primitive variables in right interface node
    // of the left ghost cell, and load the primitive variables
    //  in the left interface node of the right ghost cell
    //
    // these values must be populated to be able to pose/solve
    // Riemann problems at the left/right domain boundaries
    int i = grid.NumberOfCells();
    density_left = grid.Cell(i).Node(0).GetDensity();
    velocity_left = grid.Cell(i).Node(0).GetVelocityX();
    pressure_left = grid.Cell(i).Node(0).GetPressure();
    grid.Cell(i).LeftInterface().SetPrim(density_left,
                                         velocity_left,
                                         pressure_left);
    
    density_right = grid.Cell(-1).Node(0).GetDensity();
    velocity_right = grid.Cell(-1).Node(0).GetVelocityX();
    pressure_right = grid.Cell(-1).Node(0).GetPressure();
    grid.Cell(-1).RightInterface().SetPrim(density_right,
                                           velocity_right,
                                           pressure_right);
}




/////////////////////////////////////////////////////
//
// first-order interpolation to reconstruct primitive
// vars at cell interfaces
//
/////////////////////////////////////////////////////

template <unsigned int C, unsigned int N>
void cell_bd_values_FV_1(OneDGrid<C, N>& grid)
{
    double density_left, velocity_left, pressure_left;
    double density_right, velocity_right, pressure_right;
    double density_limiter_L, velocity_limiter_L, pressure_limiter_L;
    double density_limiter_R, velocity_limiter_R, pressure_limiter_R;
    double density_limiter, velocity_limiter, pressure_limiter;
    
    
    for (int i = 0; i < grid.NumberOfCells(); i++)
    {
        // store the prim vars interpolated to the left
        // boundary (1st-order)
        
        // -------------------------------
        // compute the interpolated value at the left interface
        density_left = grid.Cell(i).Node(0).GetDensity()
                       - (1.0/4.0)*(grid.Cell(i+1).Node(0).GetDensity()
                                    - grid.Cell(i-1).Node(0).GetDensity());
        // compute the left limiter value
        density_limiter_L = phi(grid.Cell(i-1).Node(0).GetDensity(),
                              grid.Cell(i).Node(0).GetDensity(),
                              grid.Cell(i+1).Node(0).GetDensity(),
                              density_left);
        // compute the interpolated value at the right interface
        density_right = grid.Cell(i).Node(0).GetDensity()
        + (1.0/4.0)*(grid.Cell(i+1).Node(0).GetDensity()
                     - grid.Cell(i-1).Node(0).GetDensity());
        // compute the right limiter value
        density_limiter_R = phi(grid.Cell(i-1).Node(0).GetDensity(),
                              grid.Cell(i).Node(0).GetDensity(),
                              grid.Cell(i+1).Node(0).GetDensity(),
                              density_right);
        // compute the limiter value
        density_limiter = fmin(density_limiter_L, density_limiter_R);
        // compute the interface value using the limiter value
        density_left = grid.Cell(i).Node(0).GetDensity()
        - (density_limiter/4.0)*(grid.Cell(i+1).Node(0).GetDensity()
                     - grid.Cell(i-1).Node(0).GetDensity());
        // compute the interface value using the limiter value
        density_right = grid.Cell(i).Node(0).GetDensity()
        + (density_limiter/4.0)*(grid.Cell(i+1).Node(0).GetDensity()
                                 - grid.Cell(i-1).Node(0).GetDensity());
        
        
        // -------------------------------
        // compute the interpolated value at the left interface
        velocity_left = grid.Cell(i).Node(0).GetVelocityX()
                        - (1.0/4.0)*(grid.Cell(i+1).Node(0).GetVelocityX()
                                     - grid.Cell(i-1).Node(0).GetVelocityX());
        // compute the limiter value
        velocity_limiter_L = phi(grid.Cell(i-1).Node(0).GetVelocityX(),
                              grid.Cell(i).Node(0).GetVelocityX(),
                              grid.Cell(i+1).Node(0).GetVelocityX(),
                              velocity_left);
        // compute the interpolated value at the right interface
        velocity_right = grid.Cell(i).Node(0).GetVelocityX()
        + (1.0/4.0)*(grid.Cell(i+1).Node(0).GetVelocityX()
                     - grid.Cell(i-1).Node(0).GetVelocityX());
        // compute the limiter value
        velocity_limiter_R = phi(grid.Cell(i-1).Node(0).GetVelocityX(),
                               grid.Cell(i).Node(0).GetVelocityX(),
                               grid.Cell(i+1).Node(0).GetVelocityX(),
                               velocity_right);
        // compute the limiter value
        velocity_limiter = fmin(velocity_limiter_L, velocity_limiter_R);
        // compute the interface value using the limiter value
        velocity_left = grid.Cell(i).Node(0).GetVelocityX()
        - (velocity_limiter/4.0)*(grid.Cell(i+1).Node(0).GetVelocityX()
                     - grid.Cell(i-1).Node(0).GetVelocityX());
        // compute the interface value using the limiter value
        velocity_right = grid.Cell(i).Node(0).GetVelocityX()
        + (velocity_limiter/4.0)*(grid.Cell(i+1).Node(0).GetVelocityX()
                                  - grid.Cell(i-1).Node(0).GetVelocityX());
        
        
        // -------------------------------
        // compute the interpolated value at the left interface
        pressure_left = grid.Cell(i).Node(0).GetPressure()
                        - (1.0/4.0)*(grid.Cell(i+1).Node(0).GetPressure()
                                     - grid.Cell(i-1).Node(0).GetPressure());
        // compute the limiter value
        pressure_limiter_L = phi(grid.Cell(i-1).Node(0).GetPressure(),
                               grid.Cell(i).Node(0).GetPressure(),
                               grid.Cell(i+1).Node(0).GetPressure(),
                               pressure_left);
        // compute the interpolated value at the right interface
        pressure_right = grid.Cell(i).Node(0).GetPressure()
        + (1.0/4.0)*(grid.Cell(i+1).Node(0).GetPressure()
                     - grid.Cell(i-1).Node(0).GetPressure());
        // compute the limiter value
        pressure_limiter_R = phi(grid.Cell(i-1).Node(0).GetPressure(),
                               grid.Cell(i).Node(0).GetPressure(),
                               grid.Cell(i+1).Node(0).GetPressure(),
                               pressure_right);
        // compute the limiter value
        pressure_limiter = fmin(pressure_limiter_L, pressure_limiter_R);
        // compute the interface value using the limiter value
        pressure_left = grid.Cell(i).Node(0).GetPressure()
        - (pressure_limiter/4.0)*(grid.Cell(i+1).Node(0).GetPressure()
                     - grid.Cell(i-1).Node(0).GetPressure());
        // compute the interface value using the limiter value
        pressure_right = grid.Cell(i).Node(0).GetPressure()
        + (pressure_limiter/4.0)*(grid.Cell(i+1).Node(0).GetPressure()
                                  - grid.Cell(i-1).Node(0).GetPressure());
        
        
        // ==============================================================
        
        
        // store the interpolated values in the appropriate
        // interface nodes
        grid.Cell(i).LeftInterface().SetPrim(density_left,
                                             velocity_left,
                                             pressure_left);
        grid.Cell(i).RightInterface().SetPrim(density_right,
                                              velocity_right,
                                              pressure_right);
    }
    
    // load the primitive variables in right interface node
    // of the left ghost cell, and load the primitive variables
    // in the left interface node of the right ghost cell
    //
    // these values must be populated to be able to pose/solve
    // Riemann problems at the left/right domain boundaries
    int i = grid.NumberOfCells();
    density_left = grid.Cell(i).Node(0).GetDensity();
    velocity_left = grid.Cell(i).Node(0).GetVelocityX();
    pressure_left = grid.Cell(i).Node(0).GetPressure();
    grid.Cell(i).LeftInterface().SetPrim(density_left,
                                         velocity_left,
                                         pressure_left);
    
    density_right = grid.Cell(-1).Node(0).GetDensity();
    velocity_right = grid.Cell(-1).Node(0).GetVelocityX();
    pressure_right = grid.Cell(-1).Node(0).GetPressure();
    grid.Cell(-1).RightInterface().SetPrim(density_right,
                                           velocity_right,
                                           pressure_right);
}

// ADD_29

template <unsigned int C, unsigned int N>
void cell_bd_values_FR(OneDGrid<C, N>& grid)
{
    
    double density_left = 0.0, velocity_left = 0.0, pressure_left = 0.0;
    double density_right = 0.0, velocity_right = 0.0, pressure_right = 0.0;
    
    for (int i = 0; i < grid.NumberOfCells(); i++)
    {
        density_left = 0.0; velocity_left = 0.0; pressure_left = 0.0;
        density_right = 0.0; velocity_right = 0.0; pressure_right = 0.0;
        
        for (int j = 0; j < grid.NumberOfNodesPerCell(); j++)
        {
            density_left += grid.Cell(i).LeftInterface().Phi(j)*
                            grid.Cell(i).Node(j).GetMass();
            velocity_left += grid.Cell(i).LeftInterface().Phi(j)*
                             grid.Cell(i).Node(j).GetMomentumX();
            pressure_left += grid.Cell(i).LeftInterface().Phi(j)*
                             grid.Cell(i).Node(j).GetEnergy();
            
            density_right += grid.Cell(i).RightInterface().Phi(j)*
                             grid.Cell(i).Node(j).GetMass();
            velocity_right += grid.Cell(i).RightInterface().Phi(j)*
                              grid.Cell(i).Node(j).GetMomentumX();
            pressure_right += grid.Cell(i).RightInterface().Phi(j)*
                              grid.Cell(i).Node(j).GetEnergy();
        }
        
        //      /\
        // change mass, momentum and energy
        //
        // change setCons and remove SetPrimUsing...
        //      \/
        
        
        // store the interpolated values in the appropriate
        // interface nodes
        grid.Cell(i).LeftInterface().SetCons(density_left,
                                             velocity_left,
                                             pressure_left);
        grid.Cell(i).LeftInterface().SetPrimUsingCons();
        grid.Cell(i).RightInterface().SetCons(density_right,
                                              velocity_right,
                                              pressure_right);
        grid.Cell(i).RightInterface().SetPrimUsingCons();
    }
    
    // load the primitive variables in right interface node
    // of the left ghost cell, and load the primitive variables
    //  in the left interface node of the right ghost cell
    //
    // these values must be populated to be able to pose/solve
    // Riemann problems at the left/right domain boundaries
    
    
    int i = grid.NumberOfCells();
    density_left = grid.Cell(i).Node(0).GetDensity();
    velocity_left = grid.Cell(i).Node(0).GetVelocityX();
    pressure_left = grid.Cell(i).Node(0).GetPressure();
    grid.Cell(i).LeftInterface().SetPrim(density_left,
                                         velocity_left,
                                         pressure_left);
    
    density_right = grid.Cell(-1).Node(0).GetDensity();
    velocity_right = grid.Cell(-1).Node(0).GetVelocityX();
    pressure_right = grid.Cell(-1).Node(0).GetPressure();
    grid.Cell(-1).RightInterface().SetPrim(density_right,
                                           velocity_right,
                                           pressure_right);
    
//    int i = grid.NumberOfCells();
//    density_left = grid.Cell(i-1).RightInterface().GetDensity();
//    velocity_left = -grid.Cell(i-1).RightInterface().GetVelocityX();
//    pressure_left = grid.Cell(i-1).RightInterface().GetPressure();
//    grid.Cell(i).LeftInterface().SetPrim(density_left,
//                                         velocity_left,
//                                         pressure_left);
//    
//    density_right = grid.Cell(0).LeftInterface().GetDensity();
//    velocity_right = -grid.Cell(0).LeftInterface().GetVelocityX();
//    pressure_right = grid.Cell(0).LeftInterface().GetPressure();
//    grid.Cell(-1).RightInterface().SetPrim(density_right,
//                                           velocity_right,
//                                           pressure_right);
    
    // while we're at it we'll store the f(u_{j,k}) in
    // each of the nodes - we will need these to be there
    // when we compute the local fluxes later on
    double rho, u, p;
    for (int i = 0; i < grid.NumberOfCells(); i++)
    {
        for (int j = 0; j < grid.NumberOfNodesPerCell(); j++)
        {
            rho = grid.Cell(i).Node(j).GetDensity();
            u = grid.Cell(i).Node(j).GetVelocityX();
            p = grid.Cell(i).Node(j).GetPressure();
            grid.Cell(i).Node(j).SetLocalFluxFromPrim(rho, u, p);
        }
    }
}



/////////////////////////////////////////////////////
//
// Riemann solver
//
/////////////////////////////////////////////////////

void RIEMANN(double PL, double UL, double AL, double GL, double RL,
             double PR, double UR, double AR, double GR, double RR,
             std::string &PATTRN,
             double &PM, double &UM,
             double &AML, double &AMR,
             double &USL, double &UHL, double &UTL,
             double &USR, double &UHR, double &UTR,
             double TOLER=5.0e-10, int NMAX=20)
{
    //initialize some extra variables for the N-R method
    double MSL, PML, DPML;
    double MSR, PMR, DPMR;
    
    // Compute the ledt and right Riemann invariants
    
    double CL = UL + 2.0*AL/(GL - 1.0);
    double CR = UR - 2.0*AR/(GR - 1.0);
    
    // Check for vacuum state
    
    if ((CL - CR) <= 0.0)
    {
        PATTRN = "V";
        return;
    }
    
    // Make an initial estimate of the intermadiate state velocity
    // to bgein the Newton-Raphson iterative solution procedure
    
    double Z = 0.0;
    if (PL > PR)
    {
        Z = (GL-1.0)*AR*(pow((PL/PR),(0.50*(GL-1.0)/GL)))/((GR-1.0)*AL);
    }
    else
    {
        // PL <= PR
        Z = (GL-1.0)*AR*(pow((PL/PR),(0.50*(GR-1.0)/GR)))/((GR-1.0)*AL);
    }
    UM = (CL*Z + CR)/(1.0 + Z);
    
    // Note that in the case that two rarefaction waves are present
    // (pattern D) and the left and right specific heat ratios are equal
    // then an exact solution has been found and the iterative procedure
    // is not required. Check for this
    
    if ((GL == GR) && (UM >= UL) && (UM <= UR))
    {
        PATTRN = "D";
        AML = AL-0.50*(GL-1.0)*(UM-UL);
        AMR = AR+0.50*(GR-1.0)*(UM-UR);
        PM = PL*pow((AML/AL),(2.0*GL/(GL-1.0)));
        UHL = UL-AL;
        UTL = UM-AML;
        UHR = UR+AR;
        UTR = UM+AMR;
        return;
    }
    
    // Begin the Newton-Raphson iterative procedure and solve for
    // the velocity in the intermediate state. During this iterative
    // process the pressure in the intermediate state is also found
    
    int i = 1;
    do
    {
        if (i >= NMAX)
        {
            PATTRN = "F";
            return;
        }
        
        if (UM < UL)
        {
            MSL=0.25*(GL+1.0)*(UM-UL)/AL;
            MSL=MSL-sqrt(1.0+MSL*MSL);
            PML=PL*(1.0+GL*(UM-UL)*MSL/AL);
            DPML=2.0*GL*PL*MSL*MSL*MSL/(AL*(1.0+MSL*MSL));
        }
        else
        {
            // UM >= UL
            AML = AL - 0.5*(GL - 1.0)*(UM - UL);
            PML = PL*pow((AML/AL),(2.0*GL/(GL - 1.0)));
            DPML = -GL*PML/AML;
        }
        
        if (UM > UR)
        {
            MSR = 0.25*(GR + 1.0)*(UM - UR)/AR;
            MSR = MSR + sqrt(1.0 + MSR*MSR);
            PMR = PR*(1.0+GR*(UM - UR)*MSR/AR);
            DPMR = 2.0*GR*PR*MSR*MSR*MSR/(AR*(1.0 + MSR*MSR));
        }
        else
        {
            // UM <= UR
            AMR = AR + 0.50*(GR - 1.0)*(UM - UR);
            PMR = PR*pow((AMR/AR),(2.0*GR/(GR - 1.0)));
            DPMR = GR*PMR/AMR;
        }
        
        UM = UM - (PML - PMR)/(DPML - DPMR);
    } while (fabs(1.0 - PML/PMR) > TOLER);
    // COME BACK TO THIS
    PM = PML;
    
    // Determine the wave pattern of the Riemann Problem and assign
    // the other state variables for the intermdiate state.
    
    if ((UM < UL) && (UM > UR))
    {
        PATTRN = "A";
        AML = AL*sqrt(((GL + 1.0) + (GL - 1.0)*PM/PL)/((GL + 1.0) + (GL - 1.0)*PL/PM));
        USL = UL + MSL*AL;
        AMR = AR*sqrt(((GR + 1.0) + (GR - 1.0)*PM/PR)/((GR + 1.0) + (GR - 1.0)*PR/PM));
        USR = UR + MSR*AR;
    }
    else if ((UM < UL) && (UM <= UR))
    {
        PATTRN = "B";
        AML = AL*sqrt(((GL + 1.0) + (GL - 1.0)*PM/PL)/((GL + 1.0) + (GL - 1.0)*PL/PM));
        USL = UL + MSL*AL;
        UHR = UR + AR;
        UTR = UM + AMR;
    }
    else if ((UM >= UL) && (UM > UR))
    {
        PATTRN = "C";
        UHL = UL - AL;
        UTL = UM - AML;
        AMR = AR*sqrt(((GR + 1.0) + (GR - 1.0)*PM/PR)/((GR + 1.0) + (GR - 1.0)*PR/PM));
        USR = UR + MSR*AR;
    }
    else
    {
        // ((UM >= UL) && (UM <= UR))
        PATTRN = "D";
        UHL = UL - AL;
        UTL = UM - AML;
        UHR = UR + AR;
        UTR = UM + AMR;
    }
    //return;
}




/////////////////////////////////////////////////////
//
// compute R, U and P along (x - x_interface)/t = 0
//
/////////////////////////////////////////////////////

void interface_vals(double PL, double UL, double AL, double RL,
                    double PR, double UR, double AR, double RR,
                    double Gamma,
                    std::string PATTRN,
                    double PM, double UM,
                    double AML, double AMR,
                    double USL, double UHL, double UTL,
                    double USR, double UHR, double UTR,
                    double &R_star, double &U_star, double &P_star)
{
    if (PATTRN == "V")
    {
        std::cout << "\nV\n";
    }
    
    // shock-shock wave form
    if (PATTRN == "A")
    {
        if (0.0 < USL)
        {
            R_star = RL;
            U_star = UL;
            P_star = PL;
        }
        else if ((USL <= 0.0) && (0.0 < UM))
        {
            R_star = Gamma*PM/AML/AML;
            U_star = UM;
            P_star = PM;
        }
        else if ((UM <= 0.0) && (0.0 < USR))
        {
            R_star = Gamma*PM/AMR/AMR;
            U_star = UM;
            P_star = PM;
        }
        else
        {
            // USR <= 0.0
            R_star = RR;
            U_star = UR;
            P_star = PR;
        }
    }
    
    // shock-rarefaction wave form
    else if (PATTRN == "B")
    {
        if (0.0 < USL)
        {
            R_star = RL;
            U_star = UL;
            P_star = PL;
        }
        else if ((USL <= 0.0) && (0.0 < UM))
        {
            R_star = Gamma*PM/AML/AML;
            U_star = UM;
            P_star = PM;
        }
        else if ((UM <= 0.0) && (0.0 < UTR))
        {
            //R_star = Gamma*PM/AMR/AMR;
            R_star = RR*pow(PM/PR,1.0/Gamma);
            U_star = UM;
            P_star = PM;
        }
        else if ((UTR <= 0.0) && (0.0 <= UHR))
        {
            R_star = RR*pow(2.0/(Gamma + 1) - (Gamma - 1)/(Gamma + 1)/AR*UR, 2.0/(Gamma - 1));
            U_star = 2.0/(Gamma + 1)*((Gamma - 1)*UR/2.0 - AR);
            P_star = PR*pow(2.0/(Gamma + 1) - (Gamma - 1)/(Gamma + 1)/AR*UR, 2.0*Gamma/(Gamma - 1));
        }
        else
        {
            // UHR < 0.0
            R_star = RR;
            U_star = UR;
            P_star = PR;
        }
    }
    
    // rarefaction-shock wave form
    else if (PATTRN == "C")
    {
        if (0.0 < UHL)
        {
            R_star = RL;
            U_star = UL;
            P_star = PL;
        }
        else if ((UHL <= 0.0) && (0.0 <= UTL))
        {
            R_star = RL*pow(2.0/(Gamma + 1) + (Gamma - 1)/(Gamma + 1)/AL*UL, 2.0/(Gamma - 1));
            U_star = 2.0/(Gamma + 1)*((Gamma - 1)*UL/2.0 + AL);
            P_star = PL*pow(2.0/(Gamma + 1) + (Gamma - 1)/(Gamma + 1)/AL*UL, 2.0*Gamma/(Gamma - 1));
        }
        else if ((UTL < 0.0) && (0.0 < UM))
        {
            // R_star = Gamma*PM/AML/AML;
            R_star = RL*pow(PM/PL,1.0/Gamma);
            U_star = UM;
            P_star = PM;
        }
        else if ((UM <= 0.0) && (0.0 < USR))
        {
            R_star = Gamma*PM/AMR/AMR;
            U_star = UM;
            P_star = PM;
        }
        else
        {
            // USR < 0.0
            R_star = RR;
            U_star = UR;
            P_star = PR;
        }
    }
    
    // rarefaction-rarefaction wave form
    else
    {
        // PATTRN == "D"
        if (0.0 < UHL)
        {
            R_star = RL;
            U_star = UL;
            P_star = PL;
        }
        else if ((UHL <= 0.0) && (0.0 <= UTL))
        {
            R_star = RL*pow(2.0/(Gamma + 1) + (Gamma - 1)/(Gamma + 1)/AL*UL, 2.0/(Gamma - 1));
            U_star = 2.0/(Gamma + 1)*((Gamma - 1)*UL/2.0 + AL);
            P_star = PL*pow(2.0/(Gamma + 1) + (Gamma - 1)/(Gamma + 1)/AL*UL, 2.0*Gamma/(Gamma - 1));
        }
        else if ((UTL < 0.0) && (0.0 < UM))
        {
            //R_star = Gamma*PM/AML/AML;
            R_star = RL*pow(PM/PL,1.0/Gamma);
            U_star = UM;
            P_star = PM;
        }
        else if ((UM <= 0.0) && (0.0 < UTR))
        {
            //R_star = Gamma*PM/AMR/AMR;
            R_star = RR*pow(PM/PR,1.0/Gamma);
            U_star = UM;
            P_star = PM;
        }
        else if ((UTR <= 0.0) && (0.0 <= UHR))
        {
            R_star = RR*pow(2.0/(Gamma + 1) - (Gamma - 1)/(Gamma + 1)/AR*UR, 2.0/(Gamma - 1));
            U_star = 2.0/(Gamma + 1)*((Gamma - 1)*UR/2.0 - AR);
            P_star = PR*pow(2.0/(Gamma + 1) - (Gamma - 1)/(Gamma + 1)/AR*UR, 2.0*Gamma/(Gamma - 1));
        }
        else
        {
            // UHR < 0.0
            R_star = RR;
            U_star = UR;
            P_star = PR;
        }
    }
}


// ADD_30


/////////////////////////////////////////////////////
//
// set the interface fluxes in each of the cells
//
/////////////////////////////////////////////////////

template <unsigned int C, unsigned int N>
void compute_common_fluxes(OneDGrid<C ,N>& grid)
{
    double PL; double UL; double AL; double GL; double RL;
    double PR; double UR; double AR; double GR; double RR;
    
    std::string PATTRN;
    double PM; double UM;
    double AML; double AMR;
    double USL; double UHL; double UTL;
    double USR; double UHR; double UTR;
    
    double R_star; double U_star; double P_star;
    
    for (int i = 0; i < grid.NumberOfCells() + 1; i++)
    {
        // i = ith interface (i.e. left Bd of cell i)
        PL = grid.Cell(i-1).RightInterface().GetPressure();
        UL = grid.Cell(i-1).RightInterface().GetVelocityX();
        AL = grid.Cell(i-1).RightInterface().GetSoundSpeed();
        GL = grid.Cell(i-1).RightInterface().GetGamma();
        RL = grid.Cell(i-1).RightInterface().GetDensity();
        
        PR = grid.Cell(i).LeftInterface().GetPressure();
        UR = grid.Cell(i).LeftInterface().GetVelocityX();
        AR = grid.Cell(i).LeftInterface().GetSoundSpeed();
        GR = grid.Cell(i).LeftInterface().GetGamma();
        RR = grid.Cell(i).LeftInterface().GetDensity();
        
        RIEMANN(PL, UL, AL, GL, RL,
                PR, UR, AR, GR, RR,
                PATTRN,
                PM, UM,
                AML, AMR,
                USL, UHL, UTL,
                USR, UHR, UTR);
        
        interface_vals(PL, UL, AL, RL,
                       PR, UR, AR, RR,
                       GL,
                       PATTRN,
                       PM, UM,
                       AML, AMR,
                       USL, UHL, UTL,
                       USR, UHR, UTR,
                       R_star, U_star, P_star);
        
        grid.Cell(i-1).RightInterface().SetLocalFluxFromPrim(R_star,
                                                             U_star,
                                                             P_star);
        grid.Cell(i).LeftInterface().SetLocalFluxFromPrim(R_star,
                                                          U_star,
                                                          P_star);
        
    }
}


// CHG_31

template <unsigned int C, unsigned int N>
void compute_local_fluxes(OneDGrid<C ,N>& grid)
{
    if (grid.NumberOfNodesPerCell() == 1)
    {
        /***********************************
         *                                 *
         *          FINITE VOLUME          *
         *                                 *
         ***********************************/
        
        double local_flux_mass;
        double local_flux_momentum;
        double local_flux_energy;
        
        for (int i = 0; i < grid.NumberOfCells(); i++)
        {
            
            // compute the local mass flux
            local_flux_mass = (-1.0/grid.CellWidth())*
            (grid.Cell(i).RightInterface().GetMassFlux() -
             grid.Cell(i).LeftInterface().GetMassFlux());
            
            // compute the local momentum flux
            local_flux_momentum = (-1.0/grid.CellWidth())*
            (grid.Cell(i).RightInterface().GetMomentumFlux() -
             grid.Cell(i).LeftInterface().GetMomentumFlux());
            
            // compute the local energy flux
            local_flux_energy = (-1.0/grid.CellWidth())*
            (grid.Cell(i).RightInterface().GetEnergyFlux() -
             grid.Cell(i).LeftInterface().GetEnergyFlux());
            
            // store the fluxes in the central node
            grid.Cell(i).Node(0).SetLocalFlux(local_flux_mass,
                                              local_flux_momentum,
                                              local_flux_energy);
        }
    }
    else
    {
        /***********************************
         *                                 *
         *       FLUX-RECONSTRUCTION       *
         *                                 *
         ***********************************/
        
        // variables to aid reability
        int K = grid.NumberOfNodesPerCell();
        double dx = grid.CellWidth();
        
        // temporary variables
        double disc_mass_flux, disc_mom_flux, disc_energy_flux;
        double mass_corr_L, momentum_corr_L, energy_corr_L;
        double mass_corr_R, momentum_corr_R, energy_corr_R;
        double mass_flux[K], momentum_flux[K], energy_flux[K];
        
        for (int i = 0; i < grid.NumberOfCells(); i++)
        {
            for (int j = 0; j < K; j++)
            {
                // these variables will store the required inner
                // produsts used in computing the continuous
                // flux at (cell, node) = (i, j)
                
                // mass
                disc_mass_flux = 0.0;
                mass_corr_L = grid.Cell(i).LeftInterface().GetMassFlux();
                mass_corr_R = grid.Cell(i).RightInterface().GetMassFlux();
                
                // momentum
                disc_mom_flux = 0.0;
                momentum_corr_L = grid.Cell(i).LeftInterface().GetMomentumFlux();
                momentum_corr_R = grid.Cell(i).RightInterface().GetMomentumFlux();
                
                // energy
                disc_energy_flux = 0.0;
                energy_corr_L = grid.Cell(i).LeftInterface().GetEnergyFlux();
                energy_corr_R = grid.Cell(i).RightInterface().GetEnergyFlux();
                
                // iterate over all of the nodes and compute
                // the requred inner products
                for (int k = 0; k < K; k++)
                {
                    // mass
                    disc_mass_flux += grid.Cell(i).Node(j).Phi(k)*grid.Cell(i).Node(k).GetMassFlux();
                    mass_corr_L -= grid.Cell(i).LeftInterface().Phi(k)*grid.Cell(i).Node(k).GetMassFlux();
                    mass_corr_R -= grid.Cell(i).RightInterface().Phi(k)*grid.Cell(i).Node(k).GetMassFlux();
                    
                    // momentum
                    disc_mom_flux += grid.Cell(i).Node(j).Phi(k)*grid.Cell(i).Node(k).GetMomentumFlux();
                    momentum_corr_L -= grid.Cell(i).LeftInterface().Phi(k)*grid.Cell(i).Node(k).GetMomentumFlux();
                    momentum_corr_R -= grid.Cell(i).RightInterface().Phi(k)*grid.Cell(i).Node(k).GetMomentumFlux();
                    
                    // energy
                    disc_energy_flux += grid.Cell(i).Node(j).Phi(k)*grid.Cell(i).Node(k).GetEnergyFlux();
                    energy_corr_L -= grid.Cell(i).LeftInterface().Phi(k)*grid.Cell(i).Node(k).GetEnergyFlux();
                    energy_corr_R -= grid.Cell(i).RightInterface().Phi(k)*grid.Cell(i).Node(k).GetEnergyFlux();
                }
                
                // temporarily store the continous mass flux for node j
                mass_flux[j] =
                (-2.0/dx)*(disc_mass_flux +
                           mass_corr_L*grid.Cell(i).Node(j).LeftCorrection() +
                           mass_corr_R*grid.Cell(i).Node(j).RightCorrection());
                
                // temporarily store the continous momentum flux for node j
                momentum_flux[j] =
                (-2.0/dx)*(disc_mom_flux +
                           momentum_corr_L*grid.Cell(i).Node(j).LeftCorrection() +
                           momentum_corr_R*grid.Cell(i).Node(j).RightCorrection());
                
                // temporarily store the continous energy flux for node j
                energy_flux[j] =
                (-2.0/dx)*(disc_energy_flux +
                           energy_corr_L*grid.Cell(i).Node(j).LeftCorrection() +
                           energy_corr_R*grid.Cell(i).Node(j).RightCorrection());
            }
            
            // store the local fluxes, coming from the continuous
            // flux funcitons, in each node
            for (int j = 0; j < K; j++)
            {
                grid.Cell(i).Node(j).SetLocalFlux(mass_flux[j],
                                                  momentum_flux[j],
                                                  energy_flux[j]);
            }
        }
        
    }
}

// CHG_325

template <unsigned int C, unsigned int N>
void update_local_fluxes(OneDGrid<C ,N>& grid, int interp_method=0)
{
    /////////////////////////
    //
    // i chooses the method by which to interpolate information
    // to the boundary
    //
    // i = 0: FV (0th-order reconstruction) [order 1]
    // i = 1: FV (1st-order reconstruction) [order 2]
    // i = 2: FV (2nd-order reconstruction) [NOT IMPLEMENTED]
    // i = 3: FV (3rd-order reconstruction) [NOT IMPLEMENTED]
    //
    // i = 12: FRDG (2 interior points) [order 3]
    // i = 13: FRDG (3 interior points) [order 5]
    // i = 14: FRDG (4 interior points) [order 7]
    // i = 15: FRDG (5 interior points) [order 9]
    //
    // i = 22: FRLo (2 interior points) [order 2]
    // i = 23: FRLo (3 interior points) [order 4]
    // i = 24: FRLo (4 interior points) [order 6]
    // i = 25: FRLo (5 interior points) [order 8]
    //
    /////////////////////////
    
    switch (interp_method)
    {
        case 0:
            cell_bd_values_FV_0(grid);
            break;
        
        case 1:
            cell_bd_values_FV_1(grid);
            break;
            
        default:
            cell_bd_values_FR(grid);
            break;
    }
    compute_common_fluxes(grid);
    compute_local_fluxes(grid);
}





template <unsigned int C, unsigned int N>
void update_solution(OneDGrid<C, N> grid,
                     OneDGrid<C, N> prime_grid,
                     OneDGrid<C, N>& out_grid,
                     double step_size)
{
    // variables to hold updated cons values
    double mass, momentum, energy;
    
    // store the updated conserved vars in out_grid
    // using grid and prime_grid to compute the updated values
    for (int i = 0; i < grid.NumberOfCells(); i++)
    {
        for (int j = 0; j < grid.NumberOfNodesPerCell(); j++)
        {
            // CHG_27
            mass = grid.Cell(i).Node(j).GetMass() +
                   step_size*prime_grid.Cell(i).Node(j).GetMassFlux();
            momentum = grid.Cell(i).Node(j).GetMomentumX() +
                       step_size*prime_grid.Cell(i).Node(j).GetMomentumFlux();
            energy = grid.Cell(i).Node(j).GetEnergy() +
                     step_size*prime_grid.Cell(i).Node(j).GetEnergyFlux();
            
            out_grid.Cell(i).Node(j).SetCons(mass, momentum, energy);
        }
    }
    
    // update the prim vars in out_grid using the
    // cons vars that were just set
    out_grid.SetPrimUsingCons();
}

#endif /* FluxFunctions_h */

















