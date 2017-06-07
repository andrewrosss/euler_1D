//
//  OneDGrid.hpp
//  Euler1D
//
//  Created by Andrew Ross on 2017-04-04.
//  Copyright © 2017 Andrew Ross. All rights reserved.
//


#ifndef OneDGrid_hpp
#define OneDGrid_hpp

#include <stdio.h>

#include "OneDNode.hpp"
#include "OneDCell.hpp"
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <cassert>
#include <fstream>
#include <cstdlib>

template <unsigned int numInteriorCells, unsigned int numInteriorNodesPerCell>
class OneDGrid {
private:
    double mleftEndPt, mrightEndPt; // ADD_11
    OneDCell<numInteriorNodesPerCell> mCell[numInteriorCells]; // CHG_12
    OneDCell<numInteriorNodesPerCell> mLeftGhost, mRightGhost; // CHG_13
    double mStandardNodePositions[numInteriorNodesPerCell];
    
public:
    
    /////////////////////////
    //
    // OneDGrid()
    // Ghost
    // Cell
    // NumberOfCells
    // NumberOfNodesPerCell
    // SetLeftEndPt
    // SetRightEndPt
    // LeftEndPt
    // RightEndPt
    // CellWidth
    //
    /////////////////////////
    
    OneDGrid(double x_left, double x_right)
    {
        /////////////////////////
        //
        // override the defaut constructor
        //
        /////////////////////////
        
        mleftEndPt = x_left;    // set the left end pt
        mrightEndPt = x_right;  // set the right end pt
        // ADD_14
//        for (int i = 0; i < numInteriorNodesPerCell; i++)
//        {
//            mStandardNodePositions[i] = std_node_positions[i];
//        }
        
        // define the cell wdith
        // ADD_15
        double dx = (mrightEndPt - mleftEndPt)/((double)(numInteriorCells));
        double x0; // cell center
        
        // set ghost cell widths
        // ADD_16
        mLeftGhost.SetCellWidth(dx);
        mRightGhost.SetCellWidth(dx);
        
        // set ghost cell centers
        // ADD_17
        mLeftGhost.SetCellCenter(mleftEndPt - dx/2.0);
        mRightGhost.SetCellCenter(mrightEndPt + dx/2.0);
        
//        // set node positions in ghost cells
//        for (int i = 0; i < numInteriorNodesPerCell; i++)
//        {
//            x0 = mLeftGhost.Center();
//            mLeftGhost.Node(i).SetXCoordinate(dx*mStandardNodePositions[i]/2.0 + x0);
//            x0 = mRightGhost.Center();
//            mRightGhost.Node(i).SetXCoordinate(dx*mStandardNodePositions[i]/2.0 + x0);
//        }
        
        
        // ADD_18
        for (int i = 0; i < numInteriorCells; i++)
        {
            // set cell width
            mCell[i].SetCellWidth(dx);
            
            // set cell center
            x0 = mleftEndPt + i*dx + dx/2.0;
            mCell[i].SetCellCenter(x0);
            
//            for (int j = 0; j < numInteriorNodesPerCell; j++)
//            {
//                mCell[i].Node(j).SetXCoordinate(dx*mStandardNodePositions[j]/2.0 + x0);
//            }
            
        }
    }
    
    OneDCell<numInteriorNodesPerCell>& Ghost(int i)
    {
        /////////////////////////
        //
        // return a reference to the ghost cells
        // -1, -2, ... : going left
        // 1, 2, ... : going right
        //
        /////////////////////////
        
        if (i == -1)
        {
            return mLeftGhost;
        }
        else
        {
            return mRightGhost;
        }
    }
    
    OneDCell<numInteriorNodesPerCell>& Cell(int i)
    {
        /////////////////////////
        //
        // return a reference to the ith cell
        //
        /////////////////////////
        
        // CHG_19
        assert(i > -2);
        assert(i < ((int)(numInteriorCells)) + 1);
        if (i == -1)
        {
            return mLeftGhost;
        }
        else if (i == numInteriorCells)
        {
            return mRightGhost;
        }
        else
        {
            return mCell[i];
        }
    }
    
    int NumberOfCells()
    {
        /////////////////////////
        //
        // return the number of cells in the grid
        // -- use for looping over the cells --
        //
        /////////////////////////
        
        return numInteriorCells;
    }
    
    int NumberOfNodesPerCell()
    {
        /////////////////////////
        //
        // return the number of nodes per cell
        // --  use for looping over the nodes in a given cell --
        //
        /////////////////////////
        
        return numInteriorNodesPerCell;
    }
    
    void SetPrimUsingCons()
    {
        // CHG_20
        // set the prim vars using the (updated) cons vars
        for (int i = 0; i < numInteriorCells; i++)
        {
            for (int j = 0; j < numInteriorNodesPerCell; j++)
            {
                mCell[i].Node(j).SetPrimUsingCons();
            }
        }
        
        // also update the prim vars in the ghost cells
        for (int j = 0; j < numInteriorNodesPerCell; j++)
        {
            mLeftGhost.Node(j).SetPrimUsingCons();
            mRightGhost.Node(j).SetPrimUsingCons();
        }
    }
    
    double LeftEndPt()
    {
        /////////////////////////
        //
        // return the left grid end point
        //
        /////////////////////////
        
        return mleftEndPt;
    }
    
    double RightEndPt()
    {
        /////////////////////////
        //
        // return the right grid end point
        //
        /////////////////////////
        
        return mrightEndPt;
    }
    
    // ADD_21
    
    double CellWidth()
    {
        /////////////////////////
        //
        // return the cell widths (assuming a uniform mesh)
        //
        /////////////////////////
        
        return (mrightEndPt - mleftEndPt)/((double)(numInteriorCells));
    }
    
    // ADD_22
    
    void Initialize(double (*rho)(double),
                    double (*u)(double),
                    double (*p)(double),
                    int interp_method=0)
    {
        /////////////////////////
        //
        // use the functions rho, u and p (functions of x)
        // to initialize the given state variables at
        // each of the nodes
        //
        // interp_method only has to be given when using
        // a FR scheme to run a simulation
        //
        /////////////////////////
        
        // cell-center and cell width
        double x0;
        double dx = (mrightEndPt - mleftEndPt)/((double)(numInteriorCells));
        
        
        /***********************************
         *                                 *
         *          FINITE VOLUME          *
         *                                 *
         ***********************************/
        if (interp_method < 10)
        {
            mStandardNodePositions[0] = 0.0;
            
            // set node positions in ghost cells
            for (int i = 0; i < numInteriorNodesPerCell; i++)
            {
                x0 = mLeftGhost.Center();
                mLeftGhost.Node(i).SetXCoordinate(dx*mStandardNodePositions[i]/2.0 + x0);
                x0 = mRightGhost.Center();
                mRightGhost.Node(i).SetXCoordinate(dx*mStandardNodePositions[i]/2.0 + x0);
            }
            // set node positions in interior cells
            for (int i = 0; i < numInteriorCells; i++)
            {
                x0 = mleftEndPt + i*dx + dx/2.0;
                for (int j = 0; j < numInteriorNodesPerCell; j++)
                {
                    mCell[i].Node(j).SetXCoordinate(dx*mStandardNodePositions[j]/2.0 + x0);
                }
            }
        }
        
        
        
        
        /***********************************
         *                                 *
         *    FLUX-RECONSTRUCTION (DG)     *
         *                                 *
         ***********************************/
        else if (interp_method == 12)
        {
            // FLUX-RECONSTRUCTION
            // (DG 2 interior points)
            
            // set the local node coordinates
            mStandardNodePositions[0] = -0.5;
            mStandardNodePositions[1] = 0.5;
            
            // derivative matrix
            double D[2][2] =
            { {-1.0, 1.0},
              {-1.0, 1.0} };
            
            // interpolating polys evaluated at the left/right interface
            double poly_left[2] = {1.5, -0.5};
            double poly_right[2] = {-0.5, 1.5};
            
            double g_Left[2] = {-1.25, 0.25};
            double g_Right[2] = {-0.25, 1.25};
            
            // set node positions in ghost cells
            for (int i = 0; i < numInteriorNodesPerCell; i++)
            {
                x0 = mLeftGhost.Center();
                mLeftGhost.Node(i).SetXCoordinate(dx*mStandardNodePositions[i]/2.0 + x0);
                x0 = mRightGhost.Center();
                mRightGhost.Node(i).SetXCoordinate(dx*mStandardNodePositions[i]/2.0 + x0);
            }
            // set node positions in interior cells
            for (int i = 0; i < numInteriorCells; i++)
            {
                x0 = mleftEndPt + i*dx + dx/2.0;
                for (int j = 0; j < numInteriorNodesPerCell; j++)
                {
                    mCell[i].Node(j).SetXCoordinate(dx*mStandardNodePositions[j]/2.0 + x0);
                }
            }
            
            // ---------------------------------------
            for (int i = 0; i < numInteriorCells; i++)
            {
                mCell[i].LeftInterface().SetPhi(poly_left, numInteriorNodesPerCell);
                for (int j = 0; j < numInteriorNodesPerCell; j++)
                {
                    mCell[i].Node(j).SetGLeft(g_Left[j]);
                    mCell[i].Node(j).SetPhi(D[j], numInteriorNodesPerCell);
                    mCell[i].Node(j).SetGRight(g_Right[j]);
                }
                mCell[i].RightInterface().SetPhi(poly_right, numInteriorNodesPerCell);
            }
        }
        
        else if (interp_method == 13)
        {
            // FLUX-RECONSTRUCTION
            // (DG 3 interior points)
            
            
            // set the local node coordinates
            mStandardNodePositions[0] = -2.0/3.0;
            mStandardNodePositions[1] = 0.0;
            mStandardNodePositions[2] = 2.0/3.0;
            
            // derivative matrix
            double D[3][3] =
            { {-2.25, 3.0, -0.75},
              {-0.75, 0.0, 0.75},
              {0.75, -3.0, 2.25} };
            
            // interpolating polys evaluated at the left/right interface
            double poly_left[3] = {1.875, -1.25, 0.375};
            double poly_right[3] = {0.375, -1.25, 1.875};
            
            double g_Left[3] = {-23.0/12.0, 3.0/4.0, 1.0/12.0};
            double g_Right[3] = {-1.0/12.0, -3.0/4.0, 23.0/12.0};
            
            // set node positions in ghost cells
            for (int i = 0; i < numInteriorNodesPerCell; i++)
            {
                x0 = mLeftGhost.Center();
                mLeftGhost.Node(i).SetXCoordinate(dx*mStandardNodePositions[i]/2.0 + x0);
                x0 = mRightGhost.Center();
                mRightGhost.Node(i).SetXCoordinate(dx*mStandardNodePositions[i]/2.0 + x0);
            }
            // set node positions in interior cells
            for (int i = 0; i < numInteriorCells; i++)
            {
                x0 = mleftEndPt + i*dx + dx/2.0;
                for (int j = 0; j < numInteriorNodesPerCell; j++)
                {
                    mCell[i].Node(j).SetXCoordinate(dx*mStandardNodePositions[j]/2.0 + x0);
                }
            }
            
            // ---------------------------------------
            for (int i = 0; i < numInteriorCells; i++)
            {
                mCell[i].LeftInterface().SetPhi(poly_left, numInteriorNodesPerCell);
                for (int j = 0; j < numInteriorNodesPerCell; j++)
                {
                    mCell[i].Node(j).SetGLeft(g_Left[j]);
                    mCell[i].Node(j).SetPhi(D[j], numInteriorNodesPerCell);
                    mCell[i].Node(j).SetGRight(g_Right[j]);
                }
                mCell[i].RightInterface().SetPhi(poly_right, numInteriorNodesPerCell);
            }
        }
        
        else if (interp_method == 14)
        {
            // FLUX-RECONSTRUCTION
            // (DG 4 interior points)
            
            // set the local node coordinates
            mStandardNodePositions[0] = -0.75;
            mStandardNodePositions[1] = -0.25;
            mStandardNodePositions[2] = 0.25;
            mStandardNodePositions[3] = 0.75;
            
            // derivative matrix
            double D[4][4] =
            { {-11.0/3.0, 6.0, -3.0, 2.0/3.0},
              {-2.0/3.0, -1.0, 2.0, -1.0/3.0},
              {1.0/3.0, -2.0, 1.0, 2.0/3.0},
              {-2.0/3.0, 3.0, -6.0, 11.0/3.0} };
            
            // interpolating polys evaluated at the left/right interface
            double poly_left[4] = {35.0/16.0, -35.0/16.0, 21.0/16.0, -5.0/16.0};
            double poly_right[4] = {-5.0/16.0, 21.0/16.0, -35.0/16.0, 35.0/16.0};
            
            double g_Left[4] = {-573.0/256.0, 337.0/256.0, -73.0/256.0, -123.0/256.0};
            double g_Right[4] = {123.0/256.0, 73.0/256.0, -337.0/256.0, 573.0/256.0};
            
            // set node positions in ghost cells
            for (int i = 0; i < numInteriorNodesPerCell; i++)
            {
                x0 = mLeftGhost.Center();
                mLeftGhost.Node(i).SetXCoordinate(dx*mStandardNodePositions[i]/2.0 + x0);
                x0 = mRightGhost.Center();
                mRightGhost.Node(i).SetXCoordinate(dx*mStandardNodePositions[i]/2.0 + x0);
            }
            // set node positions in interior cells
            for (int i = 0; i < numInteriorCells; i++)
            {
                x0 = mleftEndPt + i*dx + dx/2.0;
                for (int j = 0; j < numInteriorNodesPerCell; j++)
                {
                    mCell[i].Node(j).SetXCoordinate(dx*mStandardNodePositions[j]/2.0 + x0);
                }
            }
            
            // ---------------------------------------
            for (int i = 0; i < numInteriorCells; i++)
            {
                mCell[i].LeftInterface().SetPhi(poly_left, numInteriorNodesPerCell);
                for (int j = 0; j < numInteriorNodesPerCell; j++)
                {
                    mCell[i].Node(j).SetGLeft(g_Left[j]);
                    mCell[i].Node(j).SetPhi(D[j], numInteriorNodesPerCell);
                    mCell[i].Node(j).SetGRight(g_Right[j]);
                }
                mCell[i].RightInterface().SetPhi(poly_right, numInteriorNodesPerCell);
            }
        }
        
        else if (interp_method == 15)
        {
            // FLUX-RECONSTRUCTION
            // (DG 5 interior points)
            
            // set the local node coordinates
            mStandardNodePositions[0] = -0.8;
            mStandardNodePositions[1] = -0.4;
            mStandardNodePositions[2] = 0.0;
            mStandardNodePositions[3] = 0.4;
            mStandardNodePositions[4] = 0.8;
            
            // derivative matrix
            double D[5][5] =
            {   {-5.20833, 10., -7.5, 3.33333, -0.625},
                {-0.625, -2.08333, 3.75, -1.25, 0.208333},
                {0.208333, -1.66667, 0, 1.66667, -0.208333},
                {-0.208333, 1.25, -3.75, 2.08333, 0.625},
                {0.625, -3.33333, 7.5, -10., 5.20833}   };
            
            // interpolating polys evaluated at the left/right interface
            double poly_left[5] =
            {2.46094, -3.28125, 2.95312, -1.40625, 0.273438};
            double poly_right[5] =
            {0.273437, -1.40625, 2.95312, -3.28125, 2.46094};
            
            double g_Left[5] =
            {-2.0815, 1.5985, -0.9375, -0.2815, 0.8785};
            double g_Right[5] =
            {-0.8785, 0.2815, 0.9375, -1.5985, 2.0815};
            
            // set node positions in ghost cells
            for (int i = 0; i < numInteriorNodesPerCell; i++)
            {
                x0 = mLeftGhost.Center();
                mLeftGhost.Node(i).SetXCoordinate(dx*mStandardNodePositions[i]/2.0 + x0);
                x0 = mRightGhost.Center();
                mRightGhost.Node(i).SetXCoordinate(dx*mStandardNodePositions[i]/2.0 + x0);
            }
            // set node positions in interior cells
            for (int i = 0; i < numInteriorCells; i++)
            {
                x0 = mleftEndPt + i*dx + dx/2.0;
                for (int j = 0; j < numInteriorNodesPerCell; j++)
                {
                    mCell[i].Node(j).SetXCoordinate(dx*mStandardNodePositions[j]/2.0 + x0);
                }
            }
            
            // ---------------------------------------
            for (int i = 0; i < numInteriorCells; i++)
            {
                mCell[i].LeftInterface().SetPhi(poly_left, numInteriorNodesPerCell);
                for (int j = 0; j < numInteriorNodesPerCell; j++)
                {
                    mCell[i].Node(j).SetGLeft(g_Left[j]);
                    mCell[i].Node(j).SetPhi(D[j], numInteriorNodesPerCell);
                    mCell[i].Node(j).SetGRight(g_Right[j]);
                }
                mCell[i].RightInterface().SetPhi(poly_right, numInteriorNodesPerCell);
            }
        }
        
        
        
        
        /***********************************
         *                                 *
         *    FLUX-RECONSTRUCTION (Lo)     *
         *                                 *
         ***********************************/
        else if (interp_method == 22)
        {
            // FLUX-RECONSTRUCTION
            // (LumpLo 2 interior points)
            
            // set the local node coordinates
            mStandardNodePositions[0] = -1.0;
            mStandardNodePositions[1] = 1.0;
            
            // derivative matrix
            double D[2][2] =
            { {-0.5, 0.5},
              {-0.5, 0.5} };
            
            // interpolating polys evaluated at the left/right interface
            double poly_left[2] = {1.0, 0.0};
            double poly_right[2] = {0.0, 1.0};
            
            double g_Left[2] = {-1.0, 0.0};
            double g_Right[2] = {0.0, 1.0};
            
            // set node positions in ghost cells
            for (int i = 0; i < numInteriorNodesPerCell; i++)
            {
                x0 = mLeftGhost.Center();
                mLeftGhost.Node(i).SetXCoordinate(dx*mStandardNodePositions[i]/2.0 + x0);
                x0 = mRightGhost.Center();
                mRightGhost.Node(i).SetXCoordinate(dx*mStandardNodePositions[i]/2.0 + x0);
            }
            // set node positions in interior cells
            for (int i = 0; i < numInteriorCells; i++)
            {
                x0 = mleftEndPt + i*dx + dx/2.0;
                for (int j = 0; j < numInteriorNodesPerCell; j++)
                {
                    mCell[i].Node(j).SetXCoordinate(dx*mStandardNodePositions[j]/2.0 + x0);
                }
            }
            
            // ---------------------------------------
            for (int i = 0; i < numInteriorCells; i++)
            {
                mCell[i].LeftInterface().SetPhi(poly_left, numInteriorNodesPerCell);
                for (int j = 0; j < numInteriorNodesPerCell; j++)
                {
                    mCell[i].Node(j).SetGLeft(g_Left[j]);
                    mCell[i].Node(j).SetPhi(D[j], numInteriorNodesPerCell);
                    mCell[i].Node(j).SetGRight(g_Right[j]);
                }
                mCell[i].RightInterface().SetPhi(poly_right, numInteriorNodesPerCell);
            }
        }
        
        
        
        else if (interp_method == 23)
        {
            // FLUX-RECONSTRUCTION
            // (LumpLo 3 interior points)
            
            // set the local node coordinates
            mStandardNodePositions[0] = -1.0;
            mStandardNodePositions[1] = 0.0;
            mStandardNodePositions[2] = 1.0;
            
            // derivative matrix
            double D[3][3] =
            { {-1.5, 2.0, -0.5},
              {-0.5, 0.0, 0.5},
              {0.5, -2.0, 1.5} };
            
            // interpolating polys evaluated at the left/right interface
            double poly_left[3] = {1.0, 0.0, 0.0};
            double poly_right[3] = {0.0, 0.0, 1.0};
            
            double g_Left[3] = {-3.0, 0.0, 0.0};
            double g_Right[3] = {0.0, 0.0, 3.0};
            
            // set node positions in ghost cells
            for (int i = 0; i < numInteriorNodesPerCell; i++)
            {
                x0 = mLeftGhost.Center();
                mLeftGhost.Node(i).SetXCoordinate(dx*mStandardNodePositions[i]/2.0 + x0);
                x0 = mRightGhost.Center();
                mRightGhost.Node(i).SetXCoordinate(dx*mStandardNodePositions[i]/2.0 + x0);
            }
            // set node positions in interior cells
            for (int i = 0; i < numInteriorCells; i++)
            {
                x0 = mleftEndPt + i*dx + dx/2.0;
                for (int j = 0; j < numInteriorNodesPerCell; j++)
                {
                    mCell[i].Node(j).SetXCoordinate(dx*mStandardNodePositions[j]/2.0 + x0);
                }
            }
            
            // ---------------------------------------
            for (int i = 0; i < numInteriorCells; i++)
            {
                mCell[i].LeftInterface().SetPhi(poly_left, numInteriorNodesPerCell);
                for (int j = 0; j < numInteriorNodesPerCell; j++)
                {
                    mCell[i].Node(j).SetGLeft(g_Left[j]);
                    mCell[i].Node(j).SetPhi(D[j], numInteriorNodesPerCell);
                    mCell[i].Node(j).SetGRight(g_Right[j]);
                }
                mCell[i].RightInterface().SetPhi(poly_right, numInteriorNodesPerCell);
            }
        }
        
        
        
        
        else if (interp_method == 24)
        {
            // FLUX-RECONSTRUCTION
            // (LumpLo 4 interior points)
            
            // set the local node coordinates
            mStandardNodePositions[0] = -1.0;
            mStandardNodePositions[1] = -0.447214;
            mStandardNodePositions[2] = 0.447214;
            mStandardNodePositions[3] = 1.0;
            
            // derivative matrix
            double D[4][4] =
            { {-3.0, 4.04508, -1.54508, 0.5},
              {-0.809017, 0.0, 1.11803, -0.309017},
              {0.309017, -1.11803, 0.0, 0.809017},
              {-0.5, 1.54508, -4.04508, 3.0} };
            
            // interpolating polys evaluated at the left/right interface
            double poly_left[4] = {1.0, 0.0, 0.0, 0.0};
            double poly_right[4] = {0.0, 0.0, 0.0, 1.0};
            
            double g_Left[4] = {-6.0, 0.0, 0.0, 0.0};
            double g_Right[4] = {0.0, 0.0, 0.0, 6.0};
            
            // set node positions in ghost cells
            for (int i = 0; i < numInteriorNodesPerCell; i++)
            {
                x0 = mLeftGhost.Center();
                mLeftGhost.Node(i).SetXCoordinate(dx*mStandardNodePositions[i]/2.0 + x0);
                x0 = mRightGhost.Center();
                mRightGhost.Node(i).SetXCoordinate(dx*mStandardNodePositions[i]/2.0 + x0);
            }
            // set node positions in interior cells
            for (int i = 0; i < numInteriorCells; i++)
            {
                x0 = mleftEndPt + i*dx + dx/2.0;
                for (int j = 0; j < numInteriorNodesPerCell; j++)
                {
                    mCell[i].Node(j).SetXCoordinate(dx*mStandardNodePositions[j]/2.0 + x0);
                }
            }
            
            // ---------------------------------------
            for (int i = 0; i < numInteriorCells; i++)
            {
                mCell[i].LeftInterface().SetPhi(poly_left, numInteriorNodesPerCell);
                for (int j = 0; j < numInteriorNodesPerCell; j++)
                {
                    mCell[i].Node(j).SetGLeft(g_Left[j]);
                    mCell[i].Node(j).SetPhi(D[j], numInteriorNodesPerCell);
                    mCell[i].Node(j).SetGRight(g_Right[j]);
                }
                mCell[i].RightInterface().SetPhi(poly_right, numInteriorNodesPerCell);
            }
        }
        
    
        
        else if (interp_method == 25)
        {
            // FLUX-RECONSTRUCTION
            // (LumpLo 5 interior points)
            
            // set the local node coordinates
            mStandardNodePositions[0] = -1.0;
            mStandardNodePositions[1] = -0.6546536707079771;
            mStandardNodePositions[2] = 0.0;
            mStandardNodePositions[3] = 0.6546536707079771;
            mStandardNodePositions[4] = 1.0;
            
            // derivative matrix
            double D[5][5] =
            { {-5.0, 6.75650248872424, -8.0/3.0, 1.4101641779424268, -1.0/2.0},
                {-1.2409902530309824, 0.0, 1.7457431218879393, -0.7637626158259733, 0.2590097469690171},
              {3.0/8.0, -1.3365845776954535, 0.0, 1.336584577695453, -3.0/8.0},
                {-0.25900974696901713, 0.7637626158259735, -1.7457431218879391, 0.0, 1.2409902530309824},
              {1.0/2.0, -1.410164177942426, 8.0/3.0, -6.756502488724239, 5.0} };
            
            // interpolating polys evaluated at the left/right interface
            double poly_left[5] = {1.0, 0.0, 0.0, 0.0, 0.0};
            double poly_right[5] = {0.0, 0.0, 0.0, 0.0,1.0};
            
            double g_Left[5] = {-10.0, 0.0, 0.0, 0.0, 0.0};
            double g_Right[5] = {0.0, 0.0, 0.0, 0.0, 10.0};
            
            // set node positions in ghost cells
            for (int i = 0; i < numInteriorNodesPerCell; i++)
            {
                x0 = mLeftGhost.Center();
                mLeftGhost.Node(i).SetXCoordinate(dx*mStandardNodePositions[i]/2.0 + x0);
                x0 = mRightGhost.Center();
                mRightGhost.Node(i).SetXCoordinate(dx*mStandardNodePositions[i]/2.0 + x0);
            }
            // set node positions in interior cells
            for (int i = 0; i < numInteriorCells; i++)
            {
                x0 = mleftEndPt + i*dx + dx/2.0;
                for (int j = 0; j < numInteriorNodesPerCell; j++)
                {
                    mCell[i].Node(j).SetXCoordinate(dx*mStandardNodePositions[j]/2.0 + x0);
                }
            }
            
            // ---------------------------------------
            for (int i = 0; i < numInteriorCells; i++)
            {
                mCell[i].LeftInterface().SetPhi(poly_left, numInteriorNodesPerCell);
                for (int j = 0; j < numInteriorNodesPerCell; j++)
                {
                    mCell[i].Node(j).SetGLeft(g_Left[j]);
                    mCell[i].Node(j).SetPhi(D[j], numInteriorNodesPerCell);
                    mCell[i].Node(j).SetGRight(g_Right[j]);
                }
                mCell[i].RightInterface().SetPhi(poly_right, numInteriorNodesPerCell);
            }
        }

        
        
        // load initial values into the grid
        double x;
        for (int j = 0; j < numInteriorNodesPerCell; j++)
        {
            x = mLeftGhost.Node(j).GetXCoordinate();
            mLeftGhost.Node(j).SetNode(rho(x), u(x), p(x));
            x = mRightGhost.Node(j).GetXCoordinate();
            mRightGhost.Node(j).SetNode(rho(x), u(x), p(x));
        }
        for (int i = 0; i < numInteriorCells; i++)
        {
            for (int j = 0; j < numInteriorNodesPerCell; j++)
            {
                x = mCell[i].Node(j).GetXCoordinate();
                mCell[i].Node(j).SetNode(rho(x), u(x), p(x));
            }
        }
    }
    
    // ADD_23
};

#endif /* OneDGrid_hpp */




















