/******************************************************************************
* Name      : torque_tria6.c                                                  *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Computes the integral of the cross product of the single layer potential    *
* density with the relative position with respect to the centre of mass of    *
* particle over the surface a particle. Since the single layer potential      *
* actually corresponds to the local stress density, this yields the torque on *
* the particle around its centre of mass. Gaussian quadrature integration is  *
* used.                                                                       *
* Works for quadratic interpolation in triangular elements (6-noded triangles)*
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

#include "integral_tria6.h"
#include "torque_tria6.h"

int torque_tria6(double **mNodes, unsigned int **mElems, double *vB, 
                 double *Xp, double *T)
{
	const unsigned int NODES_IN_ELEM = 6;
	unsigned int currentNode, i, j, k, xNode, yNode, zNode; 
	double A, B, dx, dy, dz, J, W, L[2], N[6], XL[3], X[6][3];

    /* Initialize */
    T[0] = T[1] = T[2] = 0.0; 

    for(i = nElemsWall; i < nElems; i++){
        for(j = 0; j < NODES_IN_ELEM; j++){
            currentNode = mElems[i][j] - 1;
            X[j][0] = mNodes[currentNode][0];
            X[j][1] = mNodes[currentNode][1];
            X[j][2] = mNodes[currentNode][2];
        }
        
        /* Integrate using gaussian quadrature */
        for(j = 0; j < TNGAUSS; j++){
            L[0] = TGauss[j][0];
            L[1] = TGauss[j][1];
            W = TGauss[j][2];

            /* Get shape function N for these (L1,L2) values */
            shape_tria6(L,N);		
            J = X2L_tria6(X,XL,L,N);
            dx = XL[0] - Xp[0];
            dy = XL[1] - Xp[1];
            dz = XL[2] - Xp[2];
            A = W*J;

            /* Add contribution from each node in the element i */
            for(k = 0; k < NODES_IN_ELEM; k++){
                xNode = mElems[i][k] - 1;
                yNode = xNode + nNodes;
                zNode = yNode + nNodes;
                B = A*N[k];
	    
                T[0] += B*(vB[yNode]*dz - vB[zNode]*dy);
                T[1] += B*(vB[zNode]*dx - vB[xNode]*dz);
                T[2] += B*(vB[xNode]*dy - vB[yNode]*dx);
            }
        }
    }

    return 0;
}
