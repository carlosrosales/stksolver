/******************************************************************************
* File      : pressure_tria3.c                                                *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Calculates the pressure at a given point Xin[]. Assumes the solution vector *
* vB corresponds to a single layer potential solution of the Stokes problem.  *
* Works for linear interpolation in triangular elements (3-noded triangles).  *
******************************************************************************/

/******************************************************************************
* COPYRIGHT & LICENSE INFORMATION                                             *
*                                                                             *
* Copyright 2006 Carlos Rosales Fernandez and The Institute of High           *
* Performance Computing (A*STAR)                                              *
*                                                                             *
* This file is part of stkSolver.                                             *
*                                                                             *
* stkSolver is free software; you can redistribute it and/or modify           *
* it under the terms of the GNU General Public License as published by        *
* the Free Software Foundation; either version 2 of the License, or           *
* (at your option) any later version.                                         *
*                                                                             *
* stkSolver is distributed in the hope that it will be useful,                *
* but WITHOUT ANY WARRANTY; without even the implied warranty of              *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               *
* GNU General Public License for more details.                                *
*                                                                             *
* You should have received a copy of the GNU General Public License           *
* along with stkSolver; if not, write to the Free Software                    *
* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA  *
******************************************************************************/

#include "constants.h"
#include "pressure_tria3.h"

double pressure_tria3(double *Xin, double **mNodes, unsigned int **mElems, 
                     double *vB)
{
    const unsigned int ELEMS = nElems, NODES_IN_ELEM = 3;
    unsigned int currentNode, i, j, SinNode, test, xNode, yNode, zNode; 
    double dx, dy, dz, P;
    double X[3][3], Int[3][3];

    /* Initialize */
    P = 0.0;
 
    for(i = 0; i < ELEMS; i++){
        for(j = 0; j < NODES_IN_ELEM; j++){
            currentNode = mElems[i][j] - 1;
            X[j][0] = mNodes[currentNode][0];
            X[j][1] = mNodes[currentNode][1];
            X[j][2] = mNodes[currentNode][2];
        }

        /* Check for singular case */
        test = 0;
        for(j = 0; j < NODES_IN_ELEM; j++){
            dx = X[j][0] - Xin[0];
            dy = X[j][1] - Xin[1];
            dz = X[j][2] - Xin[2];
            if(dx == 0.0 && dy == 0.0 && dz == 0.0){
                test = 1;
                SinNode = j+1;
                break;
            }
        }

        if(test == 0) intHStk_tria3(X,Xin,Int);
        else intSingularHStk_tria3(SinNode,X,Xin,Int);

        /* Add contribution to pressure from all j nodes in element i */
        for(j = 0; j < NODES_IN_ELEM; j++){
            xNode = mElems[i][j] - 1;
            yNode = xNode + nNodes;
            zNode = yNode + nNodes;
            P += Int[0][j]*vB[xNode] + Int[1][j]*vB[yNode] + Int[2][j]*vB[zNode];
        }
    }

    return -0.5*P/pi;
}
