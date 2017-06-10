//
//  main.cpp
//  Euler1D
//
//  Created by Andrew Ross on 2017-03-17.
//  Copyright © 2017 Andrew Ross. All rights reserved.
//

#include <iostream>
#include <iomanip>
#include <vector>
#include "OneDNode.hpp"
#include "OneDCell.hpp"
#include "OneDGrid.hpp"
#include "ExactFlux.hpp"
#include "TimeMarching.hpp"





/******** LOCAL FUNCITONS ***********************************************/

// CASE 1 INITIAL CONDITIONS
double initial_density1(double x);
double initial_velocityX1(double x);
double initial_pressure1(double x);

// CASE 2 INITIAL CONDITIONS
double initial_density2(double x);
double initial_velocityX2(double x);
double initial_pressure2(double x);

// CASE 3 INITIAL CONDITIONS
double initial_density3(double x);
double initial_velocityX3(double x);
double initial_pressure3(double x);

// CASE 4 INITIAL CONDITIONS
double initial_density4(double x);
double initial_velocityX4(double x);
double initial_pressure4(double x);

// CASE 5 INITIAL CONDITIONS
double initial_density5(double x);
double initial_velocityX5(double x);
double initial_pressure5(double x);

// CASE 6 INITIAL CONDITIONS
double initial_density6(double x);
double initial_velocityX6(double x);
double initial_pressure6(double x);

// SIMULATION RUN TIME
double total_run_time(std::string IC);

// CHOOSE INITIAL CONDITIONS GIVEN USER INPUT
typedef double (*fptr)(double x);
fptr initial_density(std::string IC);
fptr initial_velocityX(std::string IC);
fptr initial_pressure(std::string IC);





/************************************************************************/

int main(int argc, const char * argv[]) {
    
    /******** TO RUN: *******
     *                      *
     *      make clean      *
     *      make main       *
     *      ./main  E*      *
     *                      *
     ************************/
    
    
    
    
    
    /******** SPECIFY SIMULATION PARAMETERS: ****************************/
 
    // NUMBER OF CELLS
    const unsigned int interior_cells = 20;
    
    // NUMBER OF SOLUTION POINTS per CELL
    const unsigned int nodes_per_cell = 4;
    
    // CFL NUMBER
    double CFL = 0.01;
    
    // HOW TO INTERPOLATE TO THE BOUNDARY [see note above]
    int interp_method = 24;
    
    // TIME MARCHING METHOD [1: RK1 | 2: RK2 | 4: RK4]
    int time_march_method = 4;
    
    // DOMAIN ENDPOINTS
    double L_endpt = 0.0;
    double R_endpt = 10.0;
    
    
    
    
    
    /******** GENERATE THE MESH *****************************************/
    
    // INITIAL CONDITIONS [set by user from command line]
    std::string IC = argv[1];
    double runtime = total_run_time(IC);
    double (*rho)(double) = initial_density(IC);
    double (*u)(double) = initial_velocityX(IC);
    double (*p)(double) = initial_pressure(IC);
    
    // instantiate the grids
    OneDGrid<interior_cells, nodes_per_cell> grid(L_endpt,
                                                  R_endpt);
    OneDGrid<interior_cells, nodes_per_cell> grid_hat(L_endpt,
                                                      R_endpt);
    OneDGrid<interior_cells, nodes_per_cell> grid_tilde(L_endpt,
                                                        R_endpt);
    OneDGrid<interior_cells, nodes_per_cell> grid_bar(L_endpt,
                                                      R_endpt);
    
    // INITIALIZE THE GRIDS
    grid.Initialize(rho, u, p, interp_method);
    grid_hat.Initialize(rho, u, p, interp_method);
    grid_tilde.Initialize(rho, u, p, interp_method);
    grid_bar.Initialize(rho, u, p, interp_method);
    
    

    
    
    /******** RUN THE SIMULATION ****************************************/
    
    // RUN THE SIMULATION
    double elapsed_time = 0.0;
    while (elapsed_time < runtime)
    {
        printf("%4.8f\n", elapsed_time);
        take_step(grid, grid_hat, grid_tilde, grid_bar, CFL, interp_method, time_march_method, elapsed_time, runtime);
    }
    std::cout << "\n";
    
    
    /******** OUTPUT NUMERICAL SOLUTION TO FILE *************************/
    
    // CREATE/OPEN DENSITY FILE
    std::string density_file = ("density_" + std::to_string(interp_method)
                                + "_(cells=" + std::to_string(interior_cells)
                                + ")_(CFL=" + std::to_string(CFL) + ").txt");
    std::ofstream density_init(density_file);
    std::ofstream density(density_file, std::ios::app);
    assert(density_init.is_open());
    assert(density.is_open());
    
    // CREATE/OPEN X-VELOCITY FILE
    std::string velocity_file = ("velocity_" + std::to_string(interp_method)
                                 + "_(cells=" + std::to_string(interior_cells)
                                 + ")_(CFL=" + std::to_string(CFL) + ").txt");
    std::ofstream velocity_init(velocity_file);
    std::ofstream velocity(velocity_file, std::ios::app);
    assert(velocity_init.is_open());
    assert(velocity.is_open());
    
    // CREATE/OPEN PRESSURE FILE
    std::string pressure_file = ("pressure_" + std::to_string(interp_method)
                                 + "_(cells=" + std::to_string(interior_cells)
                                 + ")_(CFL=" + std::to_string(CFL) + ").txt");
    std::ofstream pressure_init(pressure_file);
    std::ofstream pressure(pressure_file, std::ios::app);
    assert(pressure_init.is_open());
    assert(pressure.is_open());
    
    // CREATE/OPEN ERROR FILE
    std::string error_file = ("error_" + std::to_string(interp_method)
                                + "_(cells=" + std::to_string(interior_cells)
                                + ")_(CFL=" + std::to_string(CFL) + ").txt");
    std::ofstream error_init(error_file);
    std::ofstream error(error_file, std::ios::app);
    assert(error_init.is_open());
    assert(error.is_open());
    
    
    // OUTPUT DATA TO FILES
    double error_buffer = 0;
    double p_ref = 404400.0, r_ref = 4.696;
    double P, R;
    for (int i = 0; i < grid.NumberOfCells(); i++)
    {
        for (int j = 0; j < grid.NumberOfNodesPerCell(); j++)
        {
            density
            << std::setprecision(15)
            << grid.Cell(i).Node(j).GetXCoordinate()
            << "\t"
            << grid.Cell(i).Node(j).GetDensity() << "\n";
            velocity
            << grid.Cell(i).Node(j).GetXCoordinate()
            << "\t"
            << grid.Cell(i).Node(j).GetVelocityX() << "\n";
            pressure
            << grid.Cell(i).Node(j).GetXCoordinate()
            << "\t"
            << grid.Cell(i).Node(j).GetPressure() << "\n";
            
            P = grid.Cell(i).Node(j).GetPressure();
            R = grid.Cell(i).Node(j).GetDensity();
            error_buffer += pow((P/p_ref)*pow(r_ref/R, 1.4) - 1.0, 2.0);
        }
    }
    
    error_buffer /= grid.NumberOfCells()*grid.NumberOfNodesPerCell();
    error << sqrt(error_buffer);

    
    // CLOSE FILES
    density.close();
    density_init.close();
    velocity.close();
    velocity_init.close();
    pressure.close();
    pressure_init.close();
    error.close();
    error_init.close();
    
    
    std::cout << "\nthe simulation ran for : " << elapsed_time << " seconds\n\n";
    
    
    return 0;
}






double initial_density1(double x){return x <= 2.0 ? 2.281 : 1.408;}
double initial_velocityX1(double x){return x <= 2.0 ? 164.83 : 0.0;}
double initial_pressure1(double x){return x <= 2.0 ? 201170.0 : 101100.0;}

double initial_density2(double x){return x <= 2.0 ? 1.045 : 3.483;}
double initial_velocityX2(double x){return x <= 2.0 ? 200.0 : 200.0;}
double initial_pressure2(double x){return x <= 2.0 ? 300000.0 : 300000.0;}

double initial_density3(double x){return x <= 5.0 ? 1.598 : 2.787;}
double initial_velocityX3(double x){return x <= 5.0 ? -383.64 : -216.97;}
double initial_pressure3(double x){return x <= 5.0 ? 91880.0 : 200000.0;}

double initial_density4(double x){return x <= 5.0 ? 4.696 : 1.408;}
double initial_velocityX4(double x){return x <= 5.0 ? 0.0 : 0.0;}
double initial_pressure4(double x){return x <= 5.0 ? 404400.0 : 101100.0;}

double initial_pressure5(double x)
{
    double pL = 404400.0, pR = 101100.0;
    if (x <= 6.0)
    {
        return pL;
    }
    else if (x >= 8.5)
    {
        return pR;
    }
    else
    {
        return (pR + pL)/2.0 - (pR - pL)*cos(2.0*M_PI*(x - 6)/5.0)/2.0;
    }
}
double initial_density5(double x)
{
    if (x <= 6.0)
    {
        return 4.696;
    }
    else if (x >= 8.5)
    {
        return 1.744;
    }
    else
    {
        return 4.696*pow(initial_pressure5(x)/404400.0, 5.0/7.0);
    }
}
double initial_velocityX5(double x)
{
    double uL = 0.0, uR = 311.92, aL = sqrt(7.0*404400.0/5.0/4.696);
    if (x <= 6.0)
    {
        return uL;
    }
    else if (x >= 8.5)
    {
        return uR;
    }
    else
    {
        return uL +
        2.0*aL/(7.0/5.0 - 1.0) -
        2.0/(7.0/5.0 - 1.0)*sqrt(7.0*initial_pressure5(x)/5.0/initial_density5(x));
    }
}

double initial_pressure6(double x)
{
    double pL = 404400.0, pR = 101100.0;
    if (x <= 6.0)
    {
        return pL;
    }
    else if (x >= 8.5)
    {
        return pR;
    }
    else
    {
        return (pR - pL)*(tanh(sqrt(2)*tan(M_PI*(2.0*(x - 6.0)/2.5 - 1.0)/2.0))
                          + 1.0)/2.0 + pL;
        
    }
}
double initial_density6(double x)
{
    if (x <= 6.0)
    {
        return 4.696;
    }
    else if (x >= 8.5)
    {
        return 1.7445572954467772;
    }
    else
    {
        return 4.696*pow(initial_pressure6(x)/404400.0, 5.0/7.0);
    }
}
double initial_velocityX6(double x)
{
    double uL = 0.0, uR = 311.92, aL = sqrt(7.0*404400.0/5.0/4.696);
    if (x <= 6.0)
    {
        return uL;
    }
    else if (x >= 8.5)
    {
        return uR;
    }
    else
    {
        return uL +
        2.0*aL/(7.0/5.0 - 1.0) -
        2.0/(7.0/5.0 - 1.0)*sqrt(7.0*initial_pressure6(x)/5.0/initial_density6(x));
    }
}



fptr initial_density(std::string IC)
{
    // ADD_3
    if (IC == "E1")
    {
        return initial_density1;
    }
    else if (IC == "E2")
    {
        return initial_density2;
    }
    else if (IC == "E3")
    {
        return initial_density3;
    }
    else if (IC == "E4")
    {
        return initial_density4;
    }
    else if (IC == "E5")
    {
        return initial_density5;
    }
    else
    {
        return initial_density6;
    }
}

fptr initial_velocityX(std::string IC)
{
    // ADD_4
    if (IC == "E1")
    {
        return initial_velocityX1;
    }
    else if (IC == "E2")
    {
        return initial_velocityX2;
    }
    else if (IC == "E3")
    {
        return initial_velocityX3;
    }
    else if (IC == "E4")
    {
        return initial_velocityX4;
    }
    else if (IC == "E5")
    {
        return initial_velocityX5;
    }
    else
    {
        return initial_velocityX6;
    }
}

// ADD_6

fptr initial_pressure(std::string IC)
{
    // ADD_5
    if (IC == "E1")
    {
        return initial_pressure1;
    }
    else if (IC == "E2")
    {
        return initial_pressure2;
    }
    else if (IC == "E3")
    {
        return initial_pressure3;
    }
    else if (IC == "E4")
    {
        return initial_pressure4;
    }
    else if (IC == "E5")
    {
        return initial_pressure5;
    }
    else
    {
        return initial_pressure6;
    }
}


double total_run_time(std::string IC)
{
    // ADD_2
    if (IC == "E1")
    {
        return 0.012;
    }
    else if (IC == "E2")
    {
        return 0.025;
    }
    else if (IC == "E3")
    {
        return 0.035;
    }
    else if (IC == "E4")
    {
        return 0.007;
    }
    else
    {
        return 0.014;
    }
}



























