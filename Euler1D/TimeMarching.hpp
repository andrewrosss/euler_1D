//
//  TimeMarching.hpp
//  Euler1D
//
//  Created by Andrew Ross on 2017-04-07.
//  Copyright © 2017 Andrew Ross. All rights reserved.
//

#ifndef TimeMarching_hpp
#define TimeMarching_hpp

#include "OneDNode.hpp"
#include "OneDCell.hpp"
#include "OneDGrid.hpp"
#include "ExactFlux.hpp"
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <cassert>
#include <fstream>
#include <cstdlib>




template <unsigned int C, unsigned int N>
double next_time_step(OneDGrid<C, N>& grid, double CFL=0.2)
{
    // variables to aid readability
    double dx = grid.CellWidth();
    double abs_vel = fabs(grid.Cell(0).Node(0).GetVelocityX());
    double sspd = grid.Cell(0).Node(0).GetSoundSpeed();
    
    // initialize a minimum
    double minimum = CFL*dx/(abs_vel + sspd);
    
    // iterate over the nodes in a grid to find the min
    for (int i = 0; i < grid.NumberOfCells(); i++)
    {
        for (int j = 0; j < grid.NumberOfNodesPerCell(); j++)
        {
            dx = grid.CellWidth();
            
            // CHG_25 ??????????????????????????????????????????????????????
            // CHG_26
            abs_vel = fabs(grid.Cell(0).Node(0).GetVelocityX());
            sspd = grid.Cell(0).Node(0).GetSoundSpeed();
            
            // check and update the min
            if (CFL*dx/(abs_vel + sspd) < minimum)
            {
                minimum = CFL*dx/(abs_vel + sspd);
            }
        }
    }
    
    return minimum;
}


template <unsigned int C, unsigned int N>
void last_RK4_step(OneDGrid<C, N>& grid,
                   OneDGrid<C, N>& grid_hat,
                   OneDGrid<C, N>& grid_tilde,
                   OneDGrid<C, N>& grid_bar,
                   double step_size)
{
    // variables to aid reability
    double mass, momentum, energy;
    
    // visit each node in grid and update the cons vars
    // using the milne corrector step (last step of the RK4 method)
    for (int i = 0; i < grid.NumberOfCells(); i++)
    {
        for (int j = 0; j < grid.NumberOfNodesPerCell(); j++)
        {
            mass = grid.Cell(i).Node(j).GetMass() +
            (step_size/6.0)*(grid.Cell(i).Node(j).GetMassFlux() +
                             2.0*(grid_hat.Cell(i).Node(j).GetMassFlux() +
                                  grid_tilde.Cell(i).Node(j).GetMassFlux()) +
                             grid_bar.Cell(i).Node(j).GetMassFlux());
            
            momentum = grid.Cell(i).Node(j).GetMomentumX() +
            (step_size/6.0)*(grid.Cell(i).Node(j).GetMomentumFlux() +
                             2.0*(grid_hat.Cell(i).Node(j).GetMomentumFlux() +
                                  grid_tilde.Cell(i).Node(j).GetMomentumFlux()) +
                             grid_bar.Cell(i).Node(j).GetMomentumFlux());
            
            // ADD_24
            
            energy = grid.Cell(i).Node(j).GetEnergy() +
            (step_size/6.0)*(grid.Cell(i).Node(j).GetEnergyFlux() +
                             2.0*(grid_hat.Cell(i).Node(j).GetEnergyFlux() +
                                  grid_tilde.Cell(i).Node(j).GetEnergyFlux()) +
                             grid_bar.Cell(i).Node(j).GetEnergyFlux());
            
            // store the updated mass, momentum and energy
            grid.Cell(i).Node(j).SetCons(mass, momentum, energy);
        }
    }
    
    // set the prim vars using the updated cons vars
    grid.SetPrimUsingCons();
}



template <unsigned int C, unsigned int N>
void RK1(OneDGrid<C, N>& grid,
         double step_size,
         int interp_method)
{
    update_local_fluxes(grid, interp_method);
    update_solution(grid, grid, grid, step_size);
}



template <unsigned int C, unsigned int N>
void RK2(OneDGrid<C, N>& grid,
         OneDGrid<C, N>& grid_hat,
         double step_size,
         int interp_method)
{
    // predictor step
    update_local_fluxes(grid, interp_method);
    update_solution(grid, grid, grid_hat, 0.5*step_size);
    
    // corrector step
    update_local_fluxes(grid_hat, interp_method);
    update_solution(grid, grid_hat, grid, step_size);
}


template <unsigned int C, unsigned int N>
void RK4(OneDGrid<C, N>& grid,
         OneDGrid<C, N>& grid_hat,
         OneDGrid<C, N>& grid_tilde,
         OneDGrid<C, N>& grid_bar,
         double step_size,
         int interp_method)
{
    // Euler predictor step
    update_local_fluxes(grid, interp_method);
    update_solution(grid, grid, grid_hat, 0.5*step_size);
    
    // Euler corrector step
    update_local_fluxes(grid_hat, interp_method);
    update_solution(grid, grid_hat, grid_tilde, 0.5*step_size);
    
    // leapfrog predictor step
    update_local_fluxes(grid_tilde, interp_method);
    update_solution(grid, grid_tilde, grid_bar, step_size);
    
    // Milne corrector step
    update_local_fluxes(grid_bar, interp_method);
    last_RK4_step(grid, grid_hat, grid_tilde, grid_bar, step_size);
}




template <unsigned int C, unsigned int N>
void take_step(OneDGrid<C, N>& grid,
               OneDGrid<C, N>& grid_hat,
               OneDGrid<C, N>& grid_tilde,
               OneDGrid<C, N>& grid_bar,
               double CFL,
               int interp_method,
               int time_march_method,
               double& elapsed_time,
               double run_time)
{
    // compute the size of the next time step
    double h = next_time_step(grid, CFL);
    if (h > run_time - elapsed_time)
    {
        h = run_time - elapsed_time;
    }
    
    // choose a time marching algorithm to update the solution
    switch (time_march_method) {
        case 1:
            RK1(grid,
                h,
                interp_method);
            break;
            
        case 2:
                RK2(grid,
                grid_hat,
                h,
                interp_method);
            break;
            
        case 4:
            RK4(grid,
                grid_hat,
                grid_tilde,
                grid_bar,
                h,
                interp_method);
            break;
    }
    elapsed_time += h;
}







#endif /* TimeMarching_hpp */






















