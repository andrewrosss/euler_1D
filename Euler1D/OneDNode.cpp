//
//  1DNode.cpp
//  Euler1D
//
//  Created by Andrew Ross on 2017-03-17.
//  Copyright © 2017 Andrew Ross. All rights reserved.
//

#include "OneDNode.hpp"
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <cassert>
#include <fstream>
#include <cstdlib>

OneDNode::OneDNode()
{
    /////////////////////////
    //
    // override default constructor and set all quantites to 0
    //
    /////////////////////////
    
    // set gamma
    m_Gamma = 1.4;
    
    // set primitive vars
    m_R = 0.0; m_U = 0.0; m_P = 0.0; m_A = 0.0; m_H = 0.0;
    
    // set conserved vars
    m_Mass = 0.0; m_Mom = 0.0; m_E = 0.0;
    
    // set local flux
    m_Mass_Flux = 0.0; m_MomentumX_Flux = 0.0; m_Energy_Flux = 0.0;
    
    // set the x-coordinate for this node
    node_coords = 0.0;
}

void OneDNode::SetXCoordinate(double x)
{
    /////////////////////////
    //
    // set the x-coordinate of this node
    //
    /////////////////////////
    
    node_coords = x;
}

double OneDNode::GetXCoordinate()
{
    /////////////////////////
    //
    // return the x-coordinate of this node
    //
    /////////////////////////
    
    return node_coords;
}

void OneDNode::SetSoundSpeed()
{
    /////////////////////////
    //
    // set the sound speed given that rho, u and p
    // are set at this node
    //
    /////////////////////////
    
    m_A = sqrt(m_Gamma*m_P/m_R);
}

void OneDNode::SetEnthalpy()
{
    /////////////////////////
    //
    // set the enthalpy given that rho, u and p
    // are set at this node
    //
    /////////////////////////
    
    m_H = m_A*m_A/(m_Gamma - 1.0) + 0.5*m_U*m_U;
}

double OneDNode::GetGamma()
{
    /////////////////////////
    //
    // return the density at this node
    //
    /////////////////////////
    
    return m_Gamma;
}

double OneDNode::GetDensity()
{
    /////////////////////////
    //
    // return the density at this node
    //
    /////////////////////////
    
    return m_R;
}

double OneDNode::GetVelocityX()
{
    /////////////////////////
    //
    // return the velocity at this node
    //
    /////////////////////////
    
    return m_U;
}

double OneDNode::GetPressure()
{
    /////////////////////////
    //
    // return the pressure at this node
    //
    /////////////////////////
    
    return m_P;
}

double OneDNode::GetSoundSpeed()
{
    /////////////////////////
    //
    // return the sound speed at this node
    //
    /////////////////////////
    
    return m_A;
}

double OneDNode::GetEnthalpy()
{
    /////////////////////////
    //
    // return the enthalpy at this node
    //
    /////////////////////////
    
    return m_H;
}

double OneDNode::GetMass()
{
    /////////////////////////
    //
    // return the mass (i.e. "conserved density") at this node
    //
    /////////////////////////
    
    return m_Mass;
}

double OneDNode::GetMomentumX()
{
    /////////////////////////
    //
    // return the momentum at this node
    //
    /////////////////////////
    
    return m_Mom;
}

double OneDNode::GetEnergy()
{
    /////////////////////////
    //
    // return the energy at this node
    //
    /////////////////////////
    
    return m_E;
}

double OneDNode::GetMassFlux()
{
    /////////////////////////
    //
    // return the mass flux at this node
    //
    /////////////////////////
    
    return m_Mass_Flux;
}

double OneDNode::GetMomentumFlux()
{
    /////////////////////////
    //
    // return the momentum flux at this node
    //
    /////////////////////////
    
    return m_MomentumX_Flux;
}

double OneDNode::GetEnergyFlux()
{
    /////////////////////////
    //
    // return the energy flux at this node
    //
    /////////////////////////
    
    return m_Energy_Flux;
}


// ---------------------------------------------------------------------


void OneDNode::SetPrimUsingCons()
{
    /////////////////////////
    //
    // Set the primitive vars at this node using the
    // current conserved vars
    //
    /////////////////////////
    
    m_R = m_Mass;
    m_U = m_Mom/m_Mass;
    m_P = (m_E - m_Mom*m_Mom/m_Mass/2.0)*(m_Gamma - 1.0);
    SetSoundSpeed();
    SetEnthalpy();
}

void OneDNode::SetConsUsingPrim()
{
    /////////////////////////
    //
    // set the conserved vars at this node using the
    // current primitive vars
    //
    /////////////////////////
    
    m_Mass = m_R;
    m_Mom = m_R*m_U;
    m_E = m_P/(m_Gamma - 1.0) + m_R*m_U*m_U/2.0;
}

void OneDNode::SetCons(double rho, double mom, double energy)
{
    /////////////////////////
    //
    // Set the conserved variables at this node
    // from an external source
    //
    /////////////////////////
    
    m_Mass = rho;
    m_Mom = mom;
    m_E = energy;
}

void OneDNode::SetPrim(double rho, double u, double p)
{
    /////////////////////////
    //
    // set the primitive vars at this node
    // from an external source
    //
    /////////////////////////
    
    m_R = rho;
    m_U = u;
    m_P = p;
    SetSoundSpeed();
    SetEnthalpy();
}

void OneDNode::SetGamma(double G)
{
    /////////////////////////
    //
    // set gamma (ratio of specific heats) to G
    //
    /////////////////////////
    
    m_Gamma = G;
}

void OneDNode::SetNode(double rho, double u, double p)
{
    /////////////////////////
    //
    // set both primitive and conserved vars using the
    // given values of rho, u and p
    //
    /////////////////////////
    
    SetPrim(rho, u, p);
    SetConsUsingPrim();
}

void OneDNode::SetLocalFlux(double mass_flux, double momentum_flux,
                          double energy_flux)
{
    /////////////////////////
    //
    // set the net flux at this node
    //
    /////////////////////////
    
    m_Mass_Flux = mass_flux;
    m_MomentumX_Flux = momentum_flux;
    m_Energy_Flux = energy_flux;
}

void OneDNode::SetLocalFluxFromPrim(double rho, double u, double p)
{
    /////////////////////////
    //
    // set the net flux at this node using primitive variables
    // rho, u and p ought to be the Godunov "star states"
    //
    /////////////////////////
    
    m_Mass_Flux = rho*u;
    m_MomentumX_Flux = p + rho*u*u;
    m_Energy_Flux = rho*u*(m_Gamma*p/rho/(m_Gamma - 1.0) + u*u/2.0);
}

void OneDNode::SetPhi(double phi[10], int N)
{
    // add a way to dynamically set the length of m_Phi
    for (int i = 0; i < N; i ++)
    {
        m_Phi[i] = phi[i];
    }
}

void OneDNode::SetGLeft(double g_left)
{
    m_G_Left = g_left;
}

void OneDNode::SetGRight(double g_right)
{
    m_G_Right = g_right;
}

double OneDNode::Phi(int i)
{
    // add a check to make sure the index is in range
    return m_Phi[i];
}

double OneDNode::LeftCorrection()
{
    return m_G_Left;
}

double OneDNode::RightCorrection()
{
    return m_G_Right;
}














