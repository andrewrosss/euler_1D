//
//  OneDCell.hpp
//  Euler1D
//
//  Created by Andrew Ross on 2017-03-18.
//  Copyright © 2017 Andrew Ross. All rights reserved.
//

#ifndef OneDCell_hpp
#define OneDCell_hpp

#include "OneDNode.hpp"
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <cassert>
#include <fstream>
#include <cstdlib>

template <unsigned int numNodes>
class OneDCell {
private:
    double cell_center, cell_width; // ADD_33
    OneDNode mNodes[numNodes];
    OneDNode mLeftInterface, mRightInterface; // ADD_34
    
public:
    /////////////////////////
    //
    // OneDCell  [constructor]
    //
    // SetNodes
    // SetNodePrim
    // SetNodeCons
    // SetPrimUsingCons
    // SetConsUsingPrim
    //
    // GetDensity
    // GetVelocityX
    // GetPressure
    // GetSoundSpeed
    // GetEnthapy
    // GetMass
    // GetMomentumX
    // GetEnergy
    //
    // SetFluxLeft
    // SetFluxRight
    // MassFluxLeft
    // MomentumXFluxLeft
    // EnergyFluxLeft
    // MassFluxRight
    // MomentumXFluxRight
    // EnergyFluxRight
    //
    // SetCellWidth
    // Width
    // SetCellCenter
    // Center
    // LeftEndPoint
    // RightEndPoint
    // len
    //
    /////////////////////////
    
    
    OneDCell()
    {
        /////////////////////////
        //
        // override the default constructor
        //
        /////////////////////////
        
        cell_center = 0.0; cell_width = 2; // CHG_35
    }
    
    void SetNodes(double rho, double u, double p)
    {
        /////////////////////////
        //
        // set primitive and conserved variables of **ALL**
        // nodes in the cell using rho, u and p
        //
        /////////////////////////
        
        for (int i = 0; i < numNodes; i++)
        {
            mNodes[i].SetPrim(rho, u, p);
            mNodes[i].SetConsUsingPrim();
        }
    }
    
    OneDNode& Node(int j)
    {
        /////////////////////////
        //
        // return a reference to the jth node
        //
        /////////////////////////
        
        assert(j < numNodes);
        return mNodes[j];
    }
    
    OneDNode& LeftInterface()
    {
        /////////////////////////
        //
        // return a reference to the left boundary node
        //
        /////////////////////////
        
        return mLeftInterface;
    }
    
    OneDNode& RightInterface()
    {
        /////////////////////////
        //
        // return a reference to the right boundary node
        //
        /////////////////////////
        
        return mRightInterface;
    }
    
    // ADD_36
    
    void SetCellWidth(double dx)
    {
        /////////////////////////
        //
        // set the cell width
        //
        /////////////////////////
        
        cell_width = dx;
    }
    
    // ADD_37
    
    double Width()
    {
        /////////////////////////
        //
        // return the cell width
        //
        /////////////////////////
        
        return cell_width;
    }
    
    // ADD_38
    
    void SetCellCenter(double center)
    {
        /////////////////////////
        //
        // set the cells center of mass (mid-point in this case)
        //
        /////////////////////////
        
        // CHG_39
        cell_center = center;
    }
    
    double Center()
    {
        /////////////////////////
        //
        // return the coordinates of the cell center
        //
        /////////////////////////
        
        return cell_center;
    }
    
    // ADD_40
    
    double LeftEndPoint()
    {
        /////////////////////////
        //
        // return the coordinates of the left cell boundary
        //
        /////////////////////////
        
        return cell_center - cell_width/2.0;
    }
    
    double ReftEndPoint()
    {
        /////////////////////////
        //
        // return the coordinates of the right cell boundary
        //
        /////////////////////////
        
        return cell_center + cell_width/2.0;
    }
    
    // ADD_41
    
    friend int len(const OneDCell cell)
    {
        /////////////////////////
        //
        // return the number of node contained in this cell
        //
        /////////////////////////
        
        return numNodes;
    }
    
    
    
    
    
    
    
    
//    void SetNodes(); // pass functions
    
//    void SetNodePrim(int j, double rho, double u, double p)
//    {
//        /////////////////////////
//        //
//        // set the primitive values at node j to rho, u and p
//        //
//        /////////////////////////
//        
//        assert(j < numNodes);
//        mNodes[j].SetPrim(rho, u, p);
//    }
//    
//    void SetNodeCons(int j, double rho, double momentum, double energy)
//    {
//        /////////////////////////
//        //
//        // set the conservative values at node j to rho, momentum and energy
//        //
//        /////////////////////////
//        
//        assert(j < numNodes);
//        mNodes[j].SetCons(rho, momentum, energy);
//    }
//    
//    void SetConsUsingPrim(int j)
//    {
//        /////////////////////////
//        //
//        // set the conservative values at node j using the
//        // primitive values at node j
//        //
//        /////////////////////////
//        
//        assert(j < numNodes);
//        mNodes[j].SetConsUsingPrim();
//    }
//    
//    void SetPrimUsingCons(int j)
//    {
//        /////////////////////////
//        //
//        // set the primitive values at node j using the
//        // conserved values at node j
//        //
//        /////////////////////////
//        
//        assert(j < numNodes);
//        mNodes[j].SetPrimUsingCons();
//    }
//
//    double GetDensity(int j)
//    {
//        /////////////////////////
//        //
//        // return the density at node j
//        //
//        /////////////////////////
//        
//        assert(j < numNodes);
//        return mNodes[j].GetDensity();
//    }
//    
//    double GetVelocityX(int j)
//    {
//        /////////////////////////
//        //
//        // return the velocity at node j
//        //
//        /////////////////////////
//        
//        assert(j < numNodes);
//        return mNodes[j].GetVelocityX();
//    }
//    
//    double GetPressure(int j)
//    {
//        /////////////////////////
//        //
//        // return the pressure at node j
//        //
//        /////////////////////////
//        
//        assert(j < numNodes);
//        return mNodes[j].GetPressure();
//    }
//    
//    double GetSoundSpeed(int j)
//    {
//        /////////////////////////
//        //
//        // return the sound speed at node j
//        //
//        /////////////////////////
//        
//        assert(j < numNodes);
//        return mNodes[j].GetSoundSpeed();
//    }
//    
//    double GetEnthalpy(int j)
//    {
//        /////////////////////////
//        //
//        // return the enthalpy at node j
//        //
//        /////////////////////////
//        
//        assert(j < numNodes);
//        return mNodes[j].GetEnthalpy();
//    }
//    
//    double GetMass(int j)
//    {
//        /////////////////////////
//        //
//        // return the mass (conserved density - really just the density)
//        // at node j
//        //
//        /////////////////////////
//        
//        assert(j < numNodes);
//        return mNodes[j].GetMass();
//    }
//    
//    double GetMomentumX(int j)
//    {
//        /////////////////////////
//        //
//        // return the momentum at node j
//        //
//        /////////////////////////
//        
//        assert(j < numNodes);
//        return mNodes[j].GetMomentumX();
//    }
//    
//    double GetEnergy(int j)
//    {
//        /////////////////////////
//        //
//        // return the energry at node j
//        //
//        /////////////////////////
//        
//        assert(j < numNodes);
//        return mNodes[j].GetEnergy();
//    }
//    
//    void SetLeftBoundaryPrim(double r, double u, double p)
//    {
//        /////////////////////////
//        //
//        // set primitive values at the left boundary
//        //
//        /////////////////////////
//        
//        mLeftInterface.SetPrim(r, u, p);
//    }
//    
//    void SetRightBoundaryPrim(double r, double u, double p)
//    {
//        /////////////////////////
//        //
//        // set primitive values at the right boundary
//        //
//        /////////////////////////
//        
//        mRightInterface.SetPrim(r, u, p);
//    }
    
    
    
    
    
    
    

//    void SetFluxLeft(double mass_flux, double momentum_flux, double energy_flux)
//    {
//        /////////////////////////
//        //
//        // set the flux values at the left boundary of the cell
//        //
//        /////////////////////////
//        
//        mLeftInterface.SetLocalFlux(mass_flux, momentum_flux, energy_flux);
//    }
//    
//    void SetFluxLeftFromPrim(double rho, double u, double p)
//    {
//        /////////////////////////
//        //
//        // set the flux values at the left boundary
//        // of the cell using prim vars
//        //
//        /////////////////////////
//        
//        mLeftInterface.SetLocalFlux(rho*u,
//                                    p + rho*u*u,
//                                    rho*u*(Gamma*p/rho/(Gamma - 1.0) + u*u/2.0));
//    }
//    
//    void SetFluxRight(double mass_flux, double momentum_flux, double energy_flux)
//    {
//        /////////////////////////
//        //
//        // set the flux values at the right boundary of the cell
//        //
//        /////////////////////////
//        
//        mRightInterface.SetLocalFlux(mass_flux, momentum_flux, energy_flux);
//    }
//    
//    void SetFluxRightFromPrim(double rho, double u, double p)
//    {
//        /////////////////////////
//        //
//        // set the flux values at the right boundary
//        // of the cell using prim vars
//        //
//        /////////////////////////
//        
//        mRightInterface.SetLocalFlux(rho*u,
//                                    p + rho*u*u,
//                                    rho*u*(Gamma*p/rho/(Gamma - 1.0) + u*u/2.0));
//    }
//    
//    double MassFluxLeft()
//    {
//        /////////////////////////
//        //
//        // return the mass flux at the left cell boundary
//        //
//        /////////////////////////
//        
//        return flux_left[0];
//    }
//    
//    double MomentumXFluxLeft()
//    {
//        /////////////////////////
//        //
//        // return the momentum flux at the left cell boundary
//        //
//        /////////////////////////
//        
//        return flux_left[1];
//    }
//    
//    double EnergyFluxLeft()
//    {
//        /////////////////////////
//        //
//        // return the energy flux at the left cell boundary
//        //
//        /////////////////////////
//        
//        return flux_left[2];
//    }
//    
//    double MassFluxRight()
//    {
//        /////////////////////////
//        //
//        // return the mass flux at the right cell boundary
//        //
//        /////////////////////////
//        
//        return flux_right[0];
//    }
//    
//    double MomentumXFluxRight()
//    {
//        /////////////////////////
//        //
//        // return the momentum flux at the right cell boundary
//        //
//        /////////////////////////
//        
//        return flux_right[1];
//    }
//    
//    double EnergyFluxRight()
//    {
//        /////////////////////////
//        //
//        // return the energy flux at the right cell boundary
//        //
//        /////////////////////////
//        
//        return flux_right[2];
//    }
};

#endif /* OneDCell_hpp */












