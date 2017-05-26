//
//  NdArray.hpp
//  Euler1D
//
//  Created by Andrew Ross on 2017-03-17.
//  Copyright © 2017 Andrew Ross. All rights reserved.
//

#ifndef NdArray_hpp
#define NdArray_hpp

#include <stdio.h>
#include <iostream>
#include <cmath>
#include <cassert>
#include <fstream>
#include <cstdlib>

class OneDNode
{
private:
    double m_Gamma;                                      // ratio of spec heats
    double m_R, m_U, m_P, m_A, m_H;                      // primitive variables
    double m_Mass, m_Mom, m_E;                           // conserved variables
    double m_Mass_Flux, m_MomentumX_Flux, m_Energy_Flux; // fluxes at this node
    double node_coords;                                  // nodes coordinates
    double m_Phi[10];
    double m_G_Left, m_G_Right;
    
public:
    OneDNode();
    void SetXCoordinate(double x);
    double GetXCoordinate();
    void SetSoundSpeed();
    void SetEnthalpy();
    double GetGamma();
    double GetDensity();
    double GetVelocityX();
    double GetPressure();
    double GetSoundSpeed();
    double GetEnthalpy();
    double GetMass();
    double GetMomentumX();
    double GetEnergy();
    double GetMassFlux();
    double GetMomentumFlux();
    double GetEnergyFlux();
    void SetPrimUsingCons();
    void SetConsUsingPrim();
    void SetCons(double rho, double mom, double energ);
    void SetPrim(double rho, double u, double p);
    void SetGamma(double Gamma);
    void SetNode(double rho, double u, double p);
    void SetLocalFlux(double mass_flux, double momentum_flux, double energy_flux);
    void SetLocalFluxFromPrim(double rho, double u, double p);
    
    void SetPhi(double phi[10], int N);
    void SetGLeft(double g_left);
    void SetGRight(double g_right);
    double Phi(int i);
    double LeftCorrection();
    double RightCorrection();
};

#endif /* NdArray_hpp */
