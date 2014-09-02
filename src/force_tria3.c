/******************************************************************************
* Name      : force_tria3.c                                                   *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Computes the integral of the single layer potential density over a particle *
* surface. Since the single layer potential actually corresponds to the local *
* stress density, this actually yields the force on the particle. Analytical  *
* integration is used (the integral of N[j] for linear triangles is always    *
* 1/6 and J is constant).                                                     *
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

#include "force_tria3.h"

int force_tria3(double **mNodes, unsigned int **mElems, double *vB, double *F)
{
    const unsigned int  NODES_IN_ELEM = 3;
    unsigned int    currentNode, i, j, SinNode, test, xNode, yNode, zNode; 
    double   A, dx1, dx2, dy1, dy2, dz1, dz2, j1, j2, j3, X[3][3];

    /* Initialize */
    F[0] = F[1] = F[2] = 0.0;
    
    for(i = nElemsWall; i < nElems; i++){
        for(j = 0; j < NODES_IN_ELEM; j++){
            currentNode = mElems[i][j] - 1;
            X[j][0] = mNodes[currentNode][0];
            X[j][1] = mNodes[currentNode][1];
            X[j][2] = mNodes[currentNode][2];
        }
        
        /* Derivatives with respect to L1 and L2 */
        dx1 = X[0][0] - X[2][0];
        dx2 = X[1][0] - X[2][0];
        dy1 = X[0][1] - X[2][1];
        dy2 = X[1][1] - X[2][1];
        dz1 = X[0][2] - X[2][2];
        dz2 = X[1][2] - X[2][2];

        /* Get Jacobian */
        j1 = dy1*dz2 - dy2*dz1;
        j2 = dx2*dz1 - dx1*dz2;
        j3 = dx1*dy2 - dx2*dy1;
        A = sqrt(j1*j1 + j2*j2 + j3*j3)/6.0;
	
        /* Add contribution to force from element i */
        for(j = 0; j < NODES_IN_ELEM; j++){
            xNode = mElems[i][j] - 1;
            yNode = xNode + nNodes;
            zNode = yNode + nNodes;
	    
            F[0] += A*vB[xNode];
            F[1] += A*vB[yNode];
            F[2] += A*vB[zNode];
        }
    }
	
	return 0;
}
