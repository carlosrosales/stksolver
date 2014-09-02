/******************************************************************************
* Name      : intHStk_tria6.c                                                 *
* Author    : Carlos Rosales Fernandez (carlos@ihpc.a-star.edu.sg)            *
* Date      : 01-09-2006                                                      *
* Revision  : 1.0                                                             *
*******************************************************************************
* DESCRIPTION                                                                 *
* Computes the integral of the Green's function of the pressure from the      *
* Stokes equation for one element and returns the integrals in Int[][] for    *
* all shape functions. Uses standard gaussian quadrature in a triangle with   *
* TNGAUSS points.                                                             *
* Works for quadratic interpolation in triangular elements (6-noded triangles)*
*                                                                             *
* Int[1][2]    : contribution to P from Phi_2 and node 3.                     *
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

int intHStk_tria6(double X[6][3], double *Xeval, double Int[3][6])
{
    const unsigned int  NODES_IN_ELEM = 6;
    unsigned int    i, j, k;
    double   A, J, r, W;
    double   B[3], dx, dy, dz, L[2], N[6], XL[3];

    /* Initialize */
    for(i = 0; i < NODES_IN_ELEM; i++){
        Int[0][i] = 0.0;
        Int[1][i] = 0.0;
        Int[2][i] = 0.0;
    }

    /* Integrate using gauss quadrature */
    for(i = 0; i < TNGAUSS; i++){
        L[0] = TGauss[i][0];
        L[1] = TGauss[i][1];
        W = TGauss[i][2];

        /* Get shape functions N for these (L1,L2) values */
        shape_tria6(L,N);
        J = X2L_tria6(X,XL,L,N);

        /* Get auxiliar quantities */
        dx = Xeval[0] - XL[0];
        dy = Xeval[1] - XL[1];
        dz = Xeval[2] - XL[2];
        r = dx*dx + dy*dy + dz*dz;
        A = J*W/(r*sqrt(r));

	    /* Get contribution from each node in the element */
	    B[0] = A*dx;
	    B[1] = A*dy;
	    B[2] = A*dz;
        for(j = 0; j < NODES_IN_ELEM; j++){
            Int[0][j] += B[0]*N[j];
            Int[1][j] += B[1]*N[j];
            Int[2][j] += B[2]*N[j];
        }
	}

	return 0;
}
